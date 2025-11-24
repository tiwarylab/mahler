# MAHLER

MAHLER (Metadynamics-Anchored Hybrid Learning for Engineering off-Rates) provides a command-line interface intended to orchestrate future simulation workflows. This repository currently offers a lightweight backbone that can be installed with `pip` and extended with new commands.

## Installation

```bash
pip install .
```

During development you can use an editable install instead:

```bash
pip install -e .
```

## Usage

After installation the `mahler` CLI becomes available on your `$PATH`.

```bash
mahler --help
```

The scaffold ships with the following example commands:

- `mahler info [--json]` prints diagnostic information about the CLI runtime.
- `mahler fold` will eventually run AlphaFold to acquire representative structures.
- `mahler rave` will run the RAVE protocol and acquire a latent space.
- `mahler ifmetad` represents infrequent metadynamics workflows.
- `mahler reweight` will reweight collected trajectories.
- `mahler mdprep` will cluster MD trajectories via `--n_clusters` and optionally align to a reference PDB via `--reference`.
- `mahler mdrun` will run MD production trajectories using inputs such as `--pdb`, `--index`, and outputs `--colvar`, `--traj`, `--log`.
