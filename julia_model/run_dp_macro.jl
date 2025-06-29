# run_dp_macro.jl - Main runner for Dynamic Programming MACRO model
# Solves the model using backward induction with vintage capital structure

include("shared.jl")  # Load common sets and parameters
include("macro/macro_data_load.jl")  # Use standard MACRO parameter loading
include("dp/dp_macro_model.jl")
include("dp/energy_cost_surrogate.jl")
include("dp/dp_utilities.jl")

using Printf
using Plots
using DataFrames

function run_dp_macro(;
  use_egm::Bool=true,          # Use endogenous grid method
  use_surrogate::Bool=false,   # Use Python surrogate (vs simple approximation)
  n_grid::Int=100,             # Grid size
  plot_results::Bool=true,
  use_linear_interpolation::Bool=true  # Use linear interpolation (vs cubic spline)
)

  println("="^60)
  println("Dynamic Programming MACRO Model")
  println("="^60)

  # Use standard MACRO parameters (already loaded from macro_data_load.jl)
  println("\nUsing standard MACRO parameters...")
  
  # Create parameter struct for DP model compatibility
  params = (
    years = year_all,
    gdp_base = gdp_base,
    k0 = k0,
    c0 = c0,
    i0 = i0,
    depr = depr,
    drate = drate,
    labor = labor,
    grow = grow,
    growth_factor = growth_factor,
    aeei_factor = aeei_factor
  )

  # Calibrate CES production function using standard MACRO coefficients
  println("\nUsing standard MACRO calibration...")
  ces = (
    a = a,  # From macro_data_load.jl
    b = b,  # From macro_data_load.jl
    ρ = rho,  # From macro_data_load.jl
    kpvs = kpvs,
    elvs = elvs
  )

  # Initialize energy cost function
  println("\nInitializing energy cost surrogate...")
  if use_surrogate
    # Use AR2 surrogate model
    energy_cost_func = (elec, nele, year) -> energy_cost_surrogate(elec, nele, year, model_type=:ar2)
  else
    energy_cost_func = simple_energy_cost
  end


  # Use standard MACRO calculated values
  # udf, newlab, finite_time_corr already calculated in macro_data_load.jl

  # Create DP model with vintage structure using standard MACRO parameters
  model = DPMacroModel(
    years=params.years,
    n_K=n_grid,
    a=ces.a,
    b=ces.b,
    ρ=ces.ρ,
    kpvs=ces.kpvs,
    elvs=ces.elvs,
    depr=params.depr,
    drate=params.drate,
    labor=params.labor,
    newlab=newlab,  # From macro_data_load.jl
    aeei_factor=params.aeei_factor,
    growth=params.grow,
    growth_factor=params.growth_factor,
    udf=udf,  # From macro_data_load.jl
    finite_time_corr=finite_time_corr,  # From macro_data_load.jl
    y0=params.gdp_base,
    k0=params.k0,
    c0=params.c0,
    i0=params.i0,
    use_linear_interpolation=use_linear_interpolation
  )

  # Create Bellman problem
  bp = BellmanProblem(model, energy_cost_func)

  # Solve using backward induction
  println("\nSolving model using backward induction...")
  solve_dp!(bp)

  # Simulate trajectory
  println("\nSimulating optimal trajectory...")
  # Initial PRODENE approximation
  PRODENE0 = model.PRODENE_base
  trajectory = simulate_trajectory(model, params.k0, params.gdp_base, PRODENE0, energy_cost_func)


# Display results
println("\n" * "="^60)
println("RESULTS")
println("="^60)

results_df = DataFrame(
  Year=trajectory.years,
  Capital=round.(trajectory.K, digits=1),
  NewCapital=round.(trajectory.KN, digits=1),
  Consumption=round.(trajectory.C, digits=1),
  Investment=round.(trajectory.I, digits=1),
  GDP=round.(trajectory.GDP_MACRO, digits=1),
  NewProduction=round.(trajectory.YN, digits=1),
  PRODENE=round.(trajectory.PRODENE_total, digits=1),
  NEWENE=round.(trajectory.NEWENE_total, digits=1),
  ElecPhys=round.(trajectory.PHYSENE_ELEC, digits=1),
  NelePhys=round.(trajectory.PHYSENE_NELE, digits=1),
  EnergyCost=round.(trajectory.EC, digits=2)
)

println(results_df)

# Calculate growth rates
println("\nGrowth Rates:")
for i in 2:length(trajectory.years)
  gdp_growth = (trajectory.GDP_MACRO[i] / trajectory.GDP_MACRO[i-1])^(1 / 10) - 1
  c_growth = (trajectory.C[i] / trajectory.C[i-1])^(1 / 10) - 1
  k_growth = (trajectory.K[i] / trajectory.K[i-1])^(1 / 10) - 1

  @printf("  %d-%d: GDP=%.1f%%, C=%.1f%%, K=%.1f%%\n",
    trajectory.years[i-1], trajectory.years[i],
    gdp_growth * 100, c_growth * 100, k_growth * 100)
end

# Key ratios
println("\nKey Ratios:")
for i in 1:length(trajectory.years)
  saving_rate = trajectory.I[i] / trajectory.GDP_MACRO[i]
  capital_output = trajectory.K[i] / trajectory.GDP_MACRO[i]
  energy_cost_share = trajectory.EC[i] / trajectory.GDP_MACRO[i]

  @printf("  %d: S/Y=%.1f%%, K/Y=%.1f, EC/Y=%.1f%%\n",
    trajectory.years[i], saving_rate * 100, capital_output, energy_cost_share * 100)
end

# Display utility
println("\nUtility: ", round(trajectory.UTILITY, digits=2))

return trajectory, results_df
end

if abspath(PROGRAM_FILE) == @__FILE__
  trajectory, results = run_dp_macro(use_egm=true, use_surrogate=false, n_grid=20, use_linear_interpolation=true)
end
