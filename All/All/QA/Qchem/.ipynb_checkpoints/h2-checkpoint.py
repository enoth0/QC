from qiskit.utils import QuantumInstance
qi = QuantumInstance(backend=backend, shots=4096, seed_simulator=42, seed_transpiler=42)