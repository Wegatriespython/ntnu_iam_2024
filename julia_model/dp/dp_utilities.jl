# dp_utilities.jl - Utility functions for DP macro model
# Includes CES production, calibration, and helper functions

using Parameters
using Roots
using ForwardDiff

# Utility functions
u(c::Float64) = log(max(c, 1e-10))
u_prime(c::Float64) = 1.0 / max(c, 1e-10)
u_prime_inv(up::Float64) = 1.0 / up

# CES production function structure
@with_kw struct CESProduction
  ρ::Float64 = -2.33      # CES exponent (esub = 0.3)
  kpvs::Float64 = 0.28     # Capital value share
  elvs::Float64 = 0.42     # Electricity value share  
  a::Float64 = 0.0         # Capital-labor coefficient
  b::Float64 = 0.0         # Energy coefficient
end

# Nested CES production function (MACRO formulation)
function ces_production(K, L, ELEC, NELE, ces::CESProduction)
  # Y = [a*K^(ρκ)*L^(ρ(1-κ)) + b*ELEC^(ρε)*NELE^(ρ(1-ε))]^(1/ρ)

  # Promote all inputs to common type
  K, L, ELEC, NELE = promote(K, L, ELEC, NELE)

  # Ensure positive inputs
  K = max(K, 0.1)
  L = max(L, 0.1)
  ELEC = max(ELEC, 0.1)
  NELE = max(NELE, 0.1)

  # Capital-labor term
  kl_term = ces.a * K^(ces.ρ * ces.kpvs) * L^(ces.ρ * (1 - ces.kpvs))

  # Energy term
  energy_term = ces.b * ELEC^(ces.ρ * ces.elvs) * NELE^(ces.ρ * (1 - ces.elvs))

  # Total production
  Y = (kl_term + energy_term)^(1 / ces.ρ)

  return Y
end

# Marginal products
function marginal_products(K::Float64, L::Float64, ELEC::Float64, NELE::Float64, ces::CESProduction)
  # Use automatic differentiation for accuracy

  # Marginal product of capital
  MPK = ForwardDiff.derivative(k -> ces_production(k, L, ELEC, NELE, ces), K)

  # Marginal product of labor  
  MPL = ForwardDiff.derivative(l -> ces_production(K, l, ELEC, NELE, ces), L)

  # Marginal product of electricity
  MPE = ForwardDiff.derivative(e -> ces_production(K, L, e, NELE, ces), ELEC)

  # Marginal product of non-electricity
  MPN = ForwardDiff.derivative(n -> ces_production(K, L, ELEC, n, ces), NELE)

  return (MPK=MPK, MPL=MPL, MPE=MPE, MPN=MPN)
end

# Calibrate CES parameters to base year data
function calibrate_ces(;
  Y0::Float64=71.0,      # Base year GDP (trillion)
  K0::Float64=198.8,     # Base year capital (trillion)
  L0::Float64=1.0,       # Normalized labor
  ELEC0::Float64=22.6,   # Base year electricity (PWh)
  NELE0::Float64=87.3,   # Base year non-electricity (PWh)
  p_elec::Float64=0.0567,# Electricity price
  p_nele::Float64=0.020, # Non-electricity price
  w::Float64=35.0,       # Wage (implied)
  r::Float64=0.10,       # Return to capital (implied)
  ρ::Float64=-2.33,
  kpvs::Float64=0.28,
  elvs::Float64=0.42
)
  # From FOCs in base year:
  # p_nele = ∂Y/∂NELE
  # p_elec = ∂Y/∂ELEC

  # Calculate coefficient b for energy using MACRO formulation
  # From FOC: p_nele = ∂Y/∂NELE
  # For Y = [a*K^(ρκ)*L^(ρ(1-κ)) + b*ELEC^(ρε)*NELE^(ρ(1-ε))]^(1/ρ)
  # We have: ∂Y/∂NELE = Y^(1-ρ) * (1-elvs) * b * ELEC^(ρ*elvs) * NELE^(ρ*(1-elvs)-1)

  # From MACRO calibration:
  b = p_nele / ((1 - elvs) * Y0^(1 - ρ) * ELEC0^(ρ * elvs) * NELE0^(ρ * (1 - elvs) - 1))

  # Calculate coefficient a for capital-labor
  # From production function: Y0 = [a*K0^(ρκ)*L0^(ρ(1-κ)) + b*ELEC0^(ρε)*NELE0^(ρ(1-ε))]^(1/ρ)
  # So: Y0^ρ = a*K0^(ρκ)*L0^(ρ(1-κ)) + b*ELEC0^(ρε)*NELE0^(ρ(1-ε))

  energy_term = b * ELEC0^(ρ * elvs) * NELE0^(ρ * (1 - elvs))
  kl_term = Y0^ρ - energy_term
  a = kl_term / (K0^(ρ * kpvs) * L0^(ρ * (1 - kpvs)))

  # Verify calibration
  ces = CESProduction(ρ=ρ, kpvs=kpvs, elvs=elvs, a=a, b=b)
  Y_check = ces_production(K0, L0, ELEC0, NELE0, ces)

  println("Calibration results:")
  println("  a = $a")
  println("  b = $b")
  println("  Y check: $Y_check (should be $Y0)")

  # Check FOCs
  mps = marginal_products(K0, L0, ELEC0, NELE0, ces)
  println("  MPE: $(mps.MPE) (should be ~$p_elec)")
  println("  MPN: $(mps.MPN) (should be ~$p_nele)")

  return ces
end

# Solve for optimal energy given prices and production target
function solve_optimal_energy(K::Float64, L::Float64, Y_target::Float64,
  p_elec::Float64, p_nele::Float64, ces::CESProduction)
  # From FOCs: MPE/MPN = p_elec/p_nele
  # This gives us the ratio ELEC/NELE

  # Also need: Y = production(K, L, ELEC, NELE)

  # Solve for energy ratio from price ratio
  price_ratio = p_elec / p_nele

  # For CES: MPE/MPN = (elvs/(1-elvs)) * (NELE/ELEC)^(1-ρ)
  # So: ELEC/NELE = (elvs/(1-elvs) / price_ratio)^(1/(1-ρ))
  energy_ratio = ((ces.elvs / (1 - ces.elvs)) / price_ratio)^(1 / (1 - ces.ρ))

  # Now solve for scale to match Y_target
  # Let NELE = x, then ELEC = energy_ratio * x
  # Solve: Y_target = production(K, L, energy_ratio*x, x)

  function production_error(nele)
    elec = energy_ratio * nele
    Y = ces_production(K, L, elec, nele, ces)
    return Y - Y_target
  end

  # Find NELE that gives target production
  nele_opt = find_zero(production_error, (1.0, 500.0))
  elec_opt = energy_ratio * nele_opt

  return elec_opt, nele_opt
end

# Calculate steady state values
function calculate_steady_state(;
  g::Float64=0.02,       # Long-run growth rate
  δ::Float64=0.05,       # Annual depreciation
  drate::Float64=0.05,   # Discount rate
  duration::Float64=10.0,
  ces::CESProduction,
  L_2080::Float64=2.0,   # Terminal labor
  p_elec::Float64=0.0567,
  p_nele::Float64=0.020
)
  # In steady state with growth g:
  # K/Y = I/Y / (g + δ)
  # C/Y = 1 - I/Y - EC/Y
  # I/Y = (g + δ) * K/Y

  # From Euler equation in steady state:
  # r = MPK - δ = g + drate

  r_ss = g + drate  # Steady state return

  # For Cobb-Douglas approximation: MPK = α * Y/K
  # So: K/Y = α / (r_ss + δ)

  # With CES, this is more complex, but approximate:
  KY_ratio = ces.kpvs / (r_ss + δ)

  # Investment ratio
  IY_ratio = (g + δ) * KY_ratio

  # Energy cost ratio (approximate from data)
  ECY_ratio = 0.07  # About 7% of GDP

  # Consumption ratio
  CY_ratio = 1 - IY_ratio - ECY_ratio

  return (
    KY_ratio=KY_ratio,
    IY_ratio=IY_ratio,
    CY_ratio=CY_ratio,
    ECY_ratio=ECY_ratio,
    r_ss=r_ss
  )
end

# NOTE: load_macro_parameters function removed - now using standard MACRO parameter loading
# from macro_data_load.jl to ensure consistency between DP and standard models
