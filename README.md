# MAHLER

MAHLER (Metadynamics-Anchored Hybrid Learning for Engineering off-Rates) provides a command-line interface for our method of reweighting IfMetaD trajectories for predicting antibody off-rates. 

## Installation

```bash
pip install .
```

## Usage

After installation the `mahler` CLI becomes available on your `$PATH`.

```bash
mahler --help
```

The package ships with the following example commands:

- `mahler info [--json]` prints diagnostic information about the CLI runtime.
- `mahler fold` will eventually run AlphaFold to acquire representative structures.
- `mahler mdprep` will cluster folded stuctures and create MD setups.
- `mahler mdrun` will run short MD for representation learning with RAVE.
- `mahler rave` will run the RAVE protocol and acquire a latent space.
- `mahler reweight` will reweight collected trajectories.


