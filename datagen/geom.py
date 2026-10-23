"""Geometry helpers shared by the datagen scripts."""
import math
from pathlib import Path

# Pyykko single-bond covalent radii (Angstrom)
RCOV = {"H": 0.32, "C": 0.75, "N": 0.71, "O": 0.63, "F": 0.64, "P": 1.11, "S": 1.03, "Cl": 0.99, "Br": 1.14,
        "Cr": 1.22, "Mn": 1.19, "Fe": 1.16, "Co": 1.11, "Ni": 1.10, "Ru": 1.25, "Rh": 1.25}
METALS = {"Cr", "Mn", "Fe", "Co", "Ni", "Ru", "Rh"}


def read_xyz(path):
    """(elements, coordinates, comment line) of an xyz file."""
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    atoms = [line.split() for line in lines[2:2 + n]]
    return [a[0] for a in atoms], [tuple(map(float, a[1:4])) for a in atoms], lines[1]


def graph(el, xyz):
    """Bonded index pairs: d < 1.3 (r_A + r_B); 1.35 (r_A + r_B) for metal bonds, which keeps long high-spin/JT
    bonds but not the second-shell atoms of five-membered chelate rings (M...C ~2.7 A). H-H pairs are ignored:
    there is no H2 ligand here, so such contacts are only steric clashes of a start structure."""
    return {(i, j) for i in range(len(el)) for j in range(i) if not el[i] == el[j] == "H"
            and math.dist(xyz[i], xyz[j]) < (1.35 if {el[i], el[j]} & METALS else 1.3) * (RCOV[el[i]] + RCOV[el[j]])}


def same_bonding(xyz_a, xyz_b):
    """True if two xyz files of the same molecule have the same bond graph."""
    return graph(*read_xyz(xyz_a)[:2]) == graph(*read_xyz(xyz_b)[:2])
