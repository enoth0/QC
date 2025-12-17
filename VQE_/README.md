# VQE_ Notebook Setup

Follow these steps to run the notebooks in this folder with the expected dependencies.

## 1) Create and activate the environment (Python 3.11)
```bash
conda create -n vqe-env python=3.9 -y
conda activate vqe-env
```

## 2) Install required packages (pinned versions)
```bash
pip install numpy==2.3.1 matplotlib==3.10.3 qiskit==1.2.4 qiskit-aer==0.17.1 qiskit-algorithms==0.3.0 qiskit-nature==0.7.2 pyscf==2.9.0
```

## 3) Optional: verify versions
```bash
python -V
pip show qiskit qiskit-aer qiskit-algorithms qiskit-nature numpy matplotlib pyscf
```

## 4) Run the notebooks
Use VS Code, Jupyter, or `jupyter notebook` from this directory and select the `vqe-env` kernel.
