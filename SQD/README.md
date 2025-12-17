# SQD Notebook Setup

Follow these steps to create the environment and install all required dependencies.

## 1) Create and activate the environment (Python 3.11)
```bash
conda create -n sqd-env python=3.11 -y
conda activate sqd-env
```

## 2) Install required packages (pinned versions)
```bash
pip install \
  qiskit==2.2.3 \
  qiskit-addon-sqd==0.12.0 \
  qiskit-nature==0.7.2 \
  qiskit-aer==0.17.1 \
  qiskit-ibm-runtime==0.43.1 \
  qiskit-algorithms==0.3.0 \
  pyscf==2.9.0 \
  rdkit==2025.9.3 \
  ffsim==0.0.63 \
  numpy==2.3.1 \
  scipy==1.16.3 \
  pandas==2.3.3 \
  matplotlib==3.10.3 \
  sympy==1.14.0 \
  networkx==3.6
```

## 3) Optional: verify versions
```bash
python -V
pip show qiskit qiskit-addon-sqd qiskit-nature qiskit-aer qiskit-ibm-runtime qiskit-algorithms pyscf rdkit ffsim numpy scipy pandas matplotlib sympy networkx
```

## 4) Run the notebooks
Launch VS Code, Jupyter, or `jupyter notebook` in this directory and select the `sqd-env` kernel.
