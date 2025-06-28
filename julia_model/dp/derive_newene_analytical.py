#!/usr/bin/env python3
"""
Derive analytical solution for NEWENE from CES production function
"""

import sympy as sp

# Define symbols
YN, KN, L = sp.symbols("YN KN L", positive=True)
ELEC, NELE = sp.symbols("ELEC NELE", positive=True)
a, b = sp.symbols("a b", positive=True)
rho = sp.symbols("rho", real=True)
kpvs, elvs = sp.symbols("kpvs elvs", positive=True)


# Capital-labor aggregate
KL_term = a * KN ** (rho * kpvs) * L ** (rho * (1 - kpvs))

# Energy aggregate
E_term = b * ELEC ** (rho * elvs) * NELE ** (rho * (1 - elvs))

# Full production function
YN_expr = (KL_term + E_term) ** (1 / rho)

print("CES Production Function:")
print(f"YN = {YN_expr}")


YN_rho = YN**rho
KL_term_value = a * KN ** (rho * kpvs) * L ** (rho * (1 - kpvs))

# Energy term must equal:
E_term_needed = YN_rho - KL_term_value
print(f"E_term_needed = {E_term_needed}")


r = sp.symbols("r", positive=True)  # ELEC/NELE ratio
NELE_solution = ((YN_rho - KL_term_value) / (b * r ** (rho * elvs))) ** (1 / rho)
ELEC_solution = r * NELE_solution

print("\nAnalytical solution:")
print(f"NELE = {NELE_solution}")
print(f"ELEC = {ELEC_solution}")

# Total NEWENE (sum of ELEC and NELE in production function terms)
NEWENE_total = ELEC_solution + NELE_solution

print(f"NEWENE_total = {NEWENE_total}")
