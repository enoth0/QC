import networkx as nx
from qiskit_algorithms import QAOA, SamplingVQE
from qiskit_algorithms.optimizers import COBYLA
from qiskit_optimization.applications import Maxcut
from qiskit.primitives import Sampler  # This will work now!
from qiskit_aer.primitives import Sampler as AerSampler # For simulation

# --- Step 1: Define the Problem ---

# Create a sample graph (a simple 4-node diamond graph)
G = nx.Graph()
G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])
n_nodes = 4

# Create the Max-Cut problem from the graph
max_cut = Maxcut(G)
qubit_op = max_cut.to_quadratic_program().to_ising()[0]

# --- Step 2 & 3: Set up QAOA and the Hybrid Loop ---

# 1. The Classical Optimizer
optimizer = COBYLA(maxiter=100)

# 2. The Quantum "Machine"
# Use the high-performance AerSampler for simulation
# Use 'Sampler()' for a basic statevector simulation
sampler = AerSampler(run_options={"shots": 1024}, 
                     default_shots=1024)

# 3. The QAOA Algorithm
# We combine the optimizer, sampler, and define p (reps)
qaoa = QAOA(sampler=sampler, optimizer=optimizer, reps=1)

# We wrap this in SamplingVQE to run the optimization loop
vqe = SamplingVQE(sampler=sampler, ansatz=qaoa.ansatz, optimizer=optimizer)

# Set the Hamiltonian for the VQE to optimize
vqe.problem = max_cut

# Run the algorithm!
result = vqe.compute_minimum_eigenvalue(qubit_op)

# --- Step 4: Get the Solution ---

# The 'result' object contains the findings.
print("--- QAOA Result ---")
print(f"Optimal parameters: {result.optimal_parameters}")
print(f"Optimal cost: {result.optimal_value}")

# To get the solution bitstring, interpret the result
solution = max_cut.sample_most_likely(result.eigenstate)
print(f"\nSolution bitstring: {solution}")
print(f"Solution partitioning: {max_cut.interpret(solution)}")