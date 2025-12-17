# VQE Repository

This repo holds notebooks and slide decks for VQE-related experiments and supporting quantum workflows.

## Repository layout
- `PPT/` slide decks used for presentations.
- `QAOA/` notebooks for QAOA experiments.
- `SQD/` notebooks for SQD workflows (see `SQD/README.md` for environment setup).
- `VQE_/` notebooks for VQE runs and visualizations (see `VQE_/README.md` for environment setup).

## Install Conda (Python 3.9 compatible)
1) Download and install Miniconda (recommended) from: https://docs.conda.io/en/latest/miniconda.html
2) During install on Windows, allow it to add the `conda` command to your shell (or later run `conda init powershell`).
3) Close and reopen PowerShell, then verify: `conda --version`.

## Set up environments
- **SQD notebooks**: follow `SQD/README.md` to create `sqd-env` (Python 3.9) and install the pinned packages.
- **VQE_ notebooks**: follow `VQE_/README.md` to create `vqe-env` (Python 3.9) and install the pinned packages.
- **QAOA notebooks**: create a dedicated env as needed; match the Python 3.9 baseline used in other folders for best compatibility.

## Using the notebooks
1) Clone the repo: `git clone https://github.com/sirenecdac/VQE.git`
2) Change into the repo: `cd VQE`
3) Create/activate the env for the folder you plan to run (see above).
4) Open in VS Code or Jupyter and select the matching kernel (e.g., `sqd-env` or `vqe-env`).
5) Run the notebooks from their respective folder.

## Notes

- Keep Python at 3.9 for these environments unless you update dependencies accordingly.
