# Tempo

TEMPO (Transferable Estimation via Metadynamics of Perturbations in Off-rates) provides a command-line interface intended to orchestrate future simulation workflows. This repository currently offers a lightweight backbone that can be installed with `pip` and extended with new commands.

## Installation

```bash
pip install .
```

During development you can use an editable install instead:

```bash
pip install -e .
```

## Usage

After installation the `tempo` CLI becomes available on your `$PATH`.

```bash
tempo --help
```

The scaffold ships with the following example commands:

- `tempo info [--json]` prints diagnostic information about the CLI runtime.
- `tempo fold` will eventually run AlphaFold to acquire representative structures.
- `tempo rave` will run the RAVE protocol and acquire a latent space.
- `tempo ifmetad` represents infrequent metadynamics workflows.
- `tempo reweight` will reweight collected trajectories.
- `tempo mdprep` will cluster MD trajectories via `--n_clusters` and optionally align to a reference PDB via `--reference`.
