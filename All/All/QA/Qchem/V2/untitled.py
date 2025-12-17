import numpy as np
import math

# Import Qiskit Nature and PySCF components
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
from qiskit.primitives import Estimator
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import COBYLA

# Import PySCF components for classical calculations
from pyscf import gto, scf, fci


def get_ch4_coords(d):
    """
    Returns CH4 geometry with tetrahedral symmetry for a given C-H bond distance d (Å).
    The Carbon atom is placed at the origin.
    """
    # Calculate coordinates based on tetrahedral geometry
    a = d / math.sqrt(3.0)
    coords = [
        ["C", [0.0, 0.0, 0.0]],
        ["H", [a, a, a]],
        ["H", [-a, -a, a]],
        ["H", [-a, a, -a]],
        ["H", [a, -a, -a]],
    ]
    return coords


def get_pyscf_mol(coords, charge, spin, basis):
    """Builds a PySCF Mole object."""
    mol = gto.Mole()
    mol.atom = coords
    mol.basis = basis
    mol.charge = charge
    mol.spin = spin
    mol.build()
    return mol


def get_qiskit_problem(coords, charge, spin, basis):
    """Generates a Qiskit Nature ElectronicStructureProblem."""
    # Format atom string for PySCFDriver
    atom_str = "; ".join([f"{atom[0]} {atom[1][0]} {atom[1][1]} {atom[1][2]}" for atom in coords])
    driver = PySCFDriver(atom=atom_str, unit=DistanceUnit.ANGSTROM,
                         charge=charge, spin=spin, basis=basis)
    problem = driver.run()
    
    print(f"Qiskit Problem Summary:")
    print(f"  Number of spatial orbitals: {problem.num_spatial_orbitals}")
    print(f"  Number of particles (alpha, beta): {problem.num_particles}")
    
    return problem


def run_uhf_pyscf(mol):
    """Runs a classical UHF calculation with PySCF."""
    mf = scf.UHF(mol)
    ehf = mf.kernel()
    return ehf


def run_fci_pyscf(mol):
    """Runs a classical Full CI calculation with PySCF (exact solution for the basis)."""
    # Run an RHF calculation first to get molecular orbitals
    mf = scf.RHF(mol).run()
    # Create an FCI solver
    cisolver = fci.FCI(mol, mf.mo_coeff)
    # kernel() returns energy and wavefunction
    efci, _ = cisolver.kernel()
    return efci


def run_vqe(problem, mapper):
    """Runs the VQE algorithm to find the ground state energy."""
    estimator = Estimator()
    # Using COBYLA as in your example, with a reasonable maxiter for a single point
    optimizer = COBYLA(maxiter=2000)

    # Setup the initial state and ansatz
    init_state = HartreeFock(problem.num_spatial_orbitals, problem.num_particles, mapper)
    ansatz = UCCSD(problem.num_spatial_orbitals, problem.num_particles, mapper, initial_state=init_state)

    vqe = VQE(estimator, ansatz, optimizer, initial_point=[0] * ansatz.num_parameters)
    
    # Map the Hamiltonian to qubits
    qubit_op = mapper.map(problem.hamiltonian.second_q_op())
    
    print(f"  Number of qubits after mapping: {qubit_op.num_qubits}")
    
    # Run VQE
    result = vqe.compute_minimum_eigenvalue(qubit_op)

    # Interpret the result
    return problem.interpret(result).total_energies[0].real


def main():
    # === Main Setup ===
    # Methane is a neutral, closed-shell molecule.
    basis = 'sto3g'
    charge = 0
    spin = 0
    # Standard equilibrium C-H bond distance for methane
    bond_distance_angstrom = 1.087

    print("Starting calculation for CH4 molecule...")
    print(f"  Basis: {basis}")
    print(f"  Charge: {charge}, Spin: {spin}")
    print(f"  C-H Bond Distance: {bond_distance_angstrom} Å")
    print("-" * 60)

    # === Calculation ===
    try:
        coords = get_ch4_coords(bond_distance_angstrom)

        # PySCF classical calculations
        mol = get_pyscf_mol(coords, charge, spin, basis)
        
        print("Running classical UHF (PySCF)...")
        uhf_energy = run_uhf_pyscf(mol)
        
        print("Running classical FCI (PySCF)...")
        fci_energy = run_fci_pyscf(mol)

        # Qiskit VQE calculation
        print("Setting up Qiskit VQE problem...")
        problem = get_qiskit_problem(coords, charge, spin, basis)
        
        # Use ParityMapper as in your example.
        # It reduces the qubit count by 2 for this problem.
        mapper = ParityMapper(num_particles=problem.num_particles)
        
        print("Running VQE (Qiskit)...")
        vqe_energy = run_vqe(problem, mapper)

        print("-" * 60)
        print("Calculation complete. Results:")
        print(f"  Classical UHF Energy:     {uhf_energy: .8f} Ha")
        print(f"  VQE (UCCSD) Energy:       {vqe_energy: .8f} Ha")
        print(f"  Classical FCI (Exact):    {fci_energy: .8f} Ha")
        print("-" * 60)
        
        # Check VQE accuracy
        vqe_error = abs(vqe_energy - fci_energy)
        print(f"  VQE Error vs FCI: {vqe_error: .2e} Ha")

    except Exception as e:
        print(f"\nAn error occurred during the calculation: {e}")


if __name__ == "__main__":
    main()
