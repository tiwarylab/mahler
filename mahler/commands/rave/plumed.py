import sys
from os import PathLike
from pathlib import Path
from typing import IO, TextIO

import mdtraj as md
import numpy as np
from af2rave.spib import SPIBResult
from mahler.utils.contact import NativeContact
from mahler.utils.topology import chain_idx_from_chain_id, check_chains

import logging
LOGGER = logging.getLogger("mahler.rave")

def sprint(*args: str, indent: int = 4, file: TextIO | None = None, **kwargs) -> None:
    """Formatted print:
    - Strips ragged leading/trailing whitespace from each line
    - Applies indentation between ...\n and \n... markers
    """
    def process_indent(text: str, indent: int = 4) -> str:
        lines = text.splitlines()
        out_lines = []
        level = 0
        for line in lines:
            line = line.strip()  # clean ragged stuff
            if line.startswith("..."):  # closing marker
                level = max(0, level - 1)
                out_lines.append(" " * (level * indent) + line)
            elif line.endswith("..."):  # opening marker
                out_lines.append(" " * (level * indent) + line)
                level += 1
            else:
                out_lines.append(" " * (level * indent) + line)
        return "\n".join(out_lines)

    pargs = [process_indent(arg, indent) for arg in args]
    if file is not None:
        kwargs["file"] = file
    print(*pargs, **kwargs)


def print_header(
        structure_file: PathLike[str] = "structure.pdb",
        file: TextIO | None = None,
) -> None:
    sprint(
        f"""UNITS LENGTH=A TIME=ps ENERGY=kcal/mol
              FLUSH STRIDE=500
              MOLINFO STRUCTURE={structure_file} MOLTYPE=protein""",
        file=file,
    )


def print_separator(file: TextIO | None = None) -> None:
    sprint("# " + "=" * 46, file=file)


def print_comment(comment: str, file: TextIO | None = None) -> None:
    sprint("# " + comment, file=file)

def print_wholemolecule(
        topology: md.Topology,
        file: TextIO | None = None,
) -> None:
    '''Print WHOLEMOLECULES directive for plumed input.'''
    protein_top = topology.subset(topology.select("protein"))
    wm_line = ["WHOLEMOLECULES"]
    for i, c in enumerate(protein_top.chains):
        chain_atoms = [atom.index for atom in c.atoms]
        start = chain_atoms[0] + 1  # 1-based
        end = chain_atoms[-1] + 1   # 1-based
        chain_range = f"ENTITY{i}={start}-{end}"
        wm_line.append(chain_range)
    sprint(" ".join(wm_line), file=file)

def print_com_dist(
        native_contact: NativeContact,
        ag_chains: list[str] = ["A"],
        ab_chains: list[str] = ["B", "C"],
        file: TextIO | None = None,
) -> None:
    ca_idx = native_contact.contacts_ca
    ag = []
    ab = []
    for cag in ag_chains:
        if cag not in ca_idx:
            raise ValueError(f"Chain '{cag}' not found in Topology.")
        ag += ca_idx[cag]
    ag = (np.array(ag) + 1).astype(str).tolist()  # 1-based
    for cab in ab_chains:
        if cab not in ca_idx:
            raise ValueError(f"Chain '{cab}' not found in Topology.")
        ab += ca_idx[cab]
    ab = (np.array(ab) + 1).astype(str).tolist()  # 1-based
    sprint(f"com_ag: CENTER ATOMS={','.join(ag)}", file=file)
    sprint(f"com_ab: CENTER ATOMS={','.join(ab)}", file=file)
    sprint("com_dist: DISTANCE ATOMS=com_ab,com_ag NOPBC", file=file)

def print_latent_variables(
        model: SPIBResult,
        index: np.ndarray,
        file: TextIO | None = None,
) -> None:

    wt = model.apparent_weight
    bias = model.apparent_bias

    n_cvs = index.shape[0]

    for i, (x, y) in enumerate(index):
        sprint(f"d_{i}: DISTANCE ATOMS={x+1},{y+1} NOPBC", file=file)
    for i, b in enumerate(bias):
        sprint(f"b_{i}: CONSTANT VALUE={b[0]:.6f}", file=file)

    fmt = ("COMBINE ...\n  LABEL=sig_{}\n  ARG=b_{}," +
        ",".join([f"d_{i}" for i in range(n_cvs)]) +
        "\n  COEFFICIENTS=1.000000," +
        ",".join(["{:.6f}" for _ in range(n_cvs)]) +
        "\n  PERIODIC=NO\n... COMBINE")
    # Assuming you want to format the coefficients with wt and include 1 in the label
    sprint(fmt.format(1, 0, *wt[0]), file=file)
    sprint(fmt.format(2, 1, *wt[1]), file=file)

def print_tail(nc: NativeContact, file: TextIO | None = None) -> None:

    dist = nc.ref_ifd[0] + 4.0
    LOGGER.debug(f"Using contacting CA: {nc.contacts_ca}")
    LOGGER.info(f"Setting committor max distance to {dist:.2f} Å (IFD + 4 Å).")

    sprint('''METAD ...
                LABEL=metad
                ARG=sig_1,sig_2
                PACE=5000
                HEIGHT=0.30
                BIASFACTOR=12
                SIGMA=0.05,0.05
                file=HILLS
                GRID_MIN=-3,-3
                GRID_MAX=3,3
                GRID_BIN=600,600
                TEMP=310.0
                CALC_RCT
                RCT_USTRIDE=1
                ACCELERATION
            ... METAD
            # =======================================================
            COMMITTOR ...
                ARG=cmap,com_dist
                STRIDE=50
                BASIN_LL1=-1,<max_dist>
                BASIN_UL1=0.5,100
            ... COMMITTOR
            # =======================================================
            # Book-keeping
            PRINT ARG=com_dist,cmap,sig_1,sig_2,metad.* FILE=COLVAR STRIDE=50
            '''.replace("<max_dist>", f"{dist:.2f}"),
            file=file,
    )

def _resolve_existing_path(path: PathLike[str]) -> Path:
    path_obj = Path(path)
    if not path_obj.exists():
        raise FileNotFoundError(path_obj)
    return path_obj


def _resolve_output(out: PathLike[str] | IO[str] | None) -> tuple[TextIO, bool]:
    """Return a writable text stream and whether it should be closed afterwards."""
    if out is None:
        return sys.stdout, False
    if hasattr(out, "write"):
        return out, False
    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    return open(out_path, "w"), True

def print_plumed(
        input_file: PathLike[str],
        ag_chains: list[str],
        ab_chains: list[str],
        spib_model: PathLike[str] | SPIBResult,
        index_file: PathLike[str] | np.ndarray,
        out: PathLike[str] | IO[str] | None = "plumed.dat",
) -> int:

    input_path = _resolve_existing_path(input_file)
    if not isinstance(spib_model, SPIBResult):
        spib_model = SPIBResult.from_file(str(spib_model))
    if not isinstance(index_file, np.ndarray):
        index = np.load(str(index_file))

    traj = md.load(str(input_path))

    top = traj.topology
    _ = check_chains(traj)          # will throw warnings if missing chain IDs

    # check that all protein chains are accounted for
    chain_map = chain_idx_from_chain_id(top, protein_only=True)
    try:
        ag_chain_idx = [chain_map[c] for c in ag_chains]
        LOGGER.debug(f"Antigen chains resolved to chain {ag_chain_idx}")
    except KeyError:
        missing = set(ag_chains).difference(set(chain_map.keys()))
        raise ValueError(f"Antigen chains not found in topology: {missing}")
    try:
        ab_chain_idx = [chain_map[c] for c in ab_chains]
        LOGGER.debug(f"Antibody chains resolved to chain {ab_chain_idx}")
    except KeyError:
        missing = set(ab_chains).difference(set(chain_map.keys()))
        raise ValueError(f"Antibody chains not found in topology: {missing}")

    nc = NativeContact(traj, [
        f"chainid {" ".join(map(str, ag_chain_idx))} and mass > 1.1",
        f"chainid {" ".join(map(str, ab_chain_idx))} and mass > 1.1"
    ])

    stream, should_close = _resolve_output(out)

    try:
        print_header(input_file, file=stream)
        print_wholemolecule(top, file=stream)
        print_separator(file=stream)
        print_comment("Interface residue distance", file=stream)
        print_com_dist(nc, ag_chains, ab_chains, file=stream)
        print_separator(file=stream)
        print_comment("Latent variables", file=stream)
        print_latent_variables(spib_model, index, file=stream)
        print_separator(file=stream)
        print_comment("Native contacts", file=stream)
        nc.print_plumed(header=False, filename=stream)
        print_separator(file=stream)
        print_tail(nc, file=stream)
    finally:
        if should_close:
            stream.close()

    return 0
