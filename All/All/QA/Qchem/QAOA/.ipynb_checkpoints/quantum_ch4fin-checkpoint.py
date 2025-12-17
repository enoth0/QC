import numpy as np
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import COBYLA
from qiskit.circuit.library import TwoLocal
from qiskit.primitives import Estimator
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.drivers import PySCFDriver
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def quantum_energy(r):
    """Compute quantum electronic energy for CH₄ at a given C–H bond length (Å)."""
    mol = f"""C 0.000000 0.000000 0.000000
H {r:.4f} 0.000000 0.000000
H {-r/3.0:.4f} {r*(2*np.sqrt(2)/3):.4f} 0.000000
H {-r/3.0:.4f} {-r*(np.sqrt(2)/3):.4f} {r*(np.sqrt(6)/3):.4f}
H {-r/3.0:.4f} {-r*(np.sqrt(2)/3):.4f} {-r*(np.sqrt(6)/3):.4f}
"""

    driver = PySCFDriver(atom=mol, basis="sto3g", spin=0, charge=0)
    problem = driver.run()

    mapper = ParityMapper()
    ansatz = TwoLocal(rotation_blocks="ry", entanglement_blocks="cz", reps=1)
    estimator = Estimator()

    vqe = VQE(estimator, ansatz, COBYLA(maxiter=100))
    solver = GroundStateEigensolver(mapper, vqe)
    result = solver.solve(problem)

    return result.total_energies[0].real


def morse_potential(r, De, re, a, E0):
    """Morse potential: E = E0 + De*(1 - exp(-a*(r - re)))^2"""
    return E0 + De * (1 - np.exp(-a * (r - re))) ** 2


def main():
    print("=== Quantum Vibrational Frequency Calculation (CH₄ Symmetric Stretch) ===")

    bond_lengths = np.linspace(0.9, 1.3, 9)
    energies = []

    for r in bond_lengths:
        print(f"Computing energy for r = {r:.4f} Å ...")
        try:
            E = quantum_energy(r)
            print(f"  ✅ Energy = {E:.8f} Hartree")
            energies.append(E)
        except Exception as e:
            print(f"  ❌ Error for r={r:.3f} Å: {e}")
            energies.append(np.nan)

    bond_lengths = np.array(bond_lengths)
    energies = np.array(energies)

    mask = ~np.isnan(energies)
    bond_lengths, energies = bond_lengths[mask], energies[mask]

    if len(energies) < 4:
        raise ValueError("Not enough valid energy points for Morse fitting.")

    
    print("\nPerforming Morse potential fit...")
    try:
        p0 = [0.2, bond_lengths[np.argmin(energies)], 1.5, min(energies)]  # initial guess
        popt, _ = curve_fit(morse_potential, bond_lengths, energies, p0=p0, maxfev=20000)
        De, re, a, E0 = popt
        print(f"✅ Morse fit succeeded:")
        print(f"   De = {De:.6f} Hartree")
        print(f"   re = {re:.6f} Å")
        print(f"   a  = {a:.6f} Å⁻¹")
        print(f"   E0 = {E0:.6f} Hartree")
    except Exception as e:
        print(f"❌ Morse fit failed: {e}")
        return

    
    hartree_to_wavenumber = 219474.63137  # cm⁻¹
    mu = (1.00784 * 12.000) / (1.00784 + 12.000) * 1.66054e-27  # reduced mass (kg)
    Eh = 4.359744e-18  # Hartree to Joules
    a_m = a * 1e10  # convert Å⁻¹ → m⁻¹

    omega_e = a_m * np.sqrt(2 * De * Eh / mu) / (2 * np.pi * 3e10)  # cm⁻¹
    omega_exe = omega_e**2 / (4 * De * hartree_to_wavenumber)
    v01 = omega_e - 2 * omega_exe

    print("\n=== Vibrational Results (Anharmonic) ===")
    print(f"ωₑ      = {omega_e:10.2f} cm⁻¹  (harmonic frequency)")
    print(f"ωₑxₑ    = {omega_exe:10.2f} cm⁻¹  (anharmonic correction)")
    print(f"ν₀→₁    = {v01:10.2f} cm⁻¹  (fundamental transition)")

    
    plt.figure(figsize=(7, 5))
    plt.plot(bond_lengths, (energies - E0) * hartree_to_wavenumber, "o", label="Quantum Energies")
    r_fit = np.linspace(min(bond_lengths), max(bond_lengths), 300)
    E_fit = morse_potential(r_fit, *popt)
    plt.plot(r_fit, (E_fit - E0) * hartree_to_wavenumber, "-", label="Morse Fit")
    plt.xlabel("C–H bond length (Å)")
    plt.ylabel("Relative Energy (cm⁻¹)")
    plt.title("CH₄ Symmetric Stretch – Morse Potential Fit")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
