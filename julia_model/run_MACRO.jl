# run_MACRO.jl - Equivalent to run_MACRO.gms

using JuMP
using GAMS

# ------------------------------------------------------------------------------
# include model specification
# ------------------------------------------------------------------------------

# include set, parameter and variable definitions that are shared across the models
include("shared.jl")

include("macro/macro_data_load.jl")
include("macro/macro_core.jl")
include("macro/macro_presolve.jl")

# ------------------------------------------------------------------------------
# solve model
# ------------------------------------------------------------------------------

# Create model and set up MACRO
model = Model(GAMS.Optimizer)
set_optimizer_attribute(model, "NLP", "CONOPT")

# Create shared variables
@variable(model, PHYSENE[sector, year_all] >= 0)
# Don't create COST_ANNUAL for standalone MACRO - this triggers the wrong constraint

# For standalone MACRO, fix energy variables
for s in sector, y in year_all
    fix(PHYSENE[s,y], enestart[(s,y)], force=true)
end

# Build MACRO model
create_macro_model!(model)
set_macro_bounds_and_initial_values!(model)
additional_macro_preprocessing!(model)

# Set objective - MAXIMIZING UTILITY
@objective(model, Max, model[:UTILITY])

# SOLVE MACRO MAXIMIZING UTILITY USING NLP
optimize!(model)

# ------------------------------------------------------------------------------
# reporting
# ------------------------------------------------------------------------------

# GDP_MACRO.L(year) = (I.L(year) + C.L(year) + EC.L(year))
GDP_MACRO = Dict(y => value(model[:I][y]) + value(model[:C][y]) + value(model[:EC][y]) for y in year_all)

# Display results in table format like dp_run_macro.jl
println("\n" * "="^60)
println("RESULTS")
println("="^60)

using CSV, DataFrames
results_df = DataFrame(
    Year = collect(year_all),
    Capital = round.([value(model[:K][y]) for y in year_all], digits=1),
    NewCapital = round.([value(model[:KN][y]) for y in year_all], digits=1),
    Consumption = round.([value(model[:C][y]) for y in year_all], digits=1),
    Investment = round.([value(model[:I][y]) for y in year_all], digits=1),
    GDP = round.([GDP_MACRO[y] for y in year_all], digits=1),
    NewProduction = round.([value(model[:YN][y]) for y in year_all], digits=1),
    PRODENE = round.([sum(value(model[:PRODENE][s,y]) for s in sector) for y in year_all], digits=1),
    NEWENE = round.([sum(value(model[:NEWENE][s,y]) for s in sector) for y in year_all], digits=1),
    ElecPhys = round.([value(PHYSENE["ELEC",y]) for y in year_all], digits=1),
    NelePhys = round.([value(PHYSENE["NELE",y]) for y in year_all], digits=1),
    EnergyCost = round.([value(model[:EC][y]) for y in year_all], digits=2)
)

println(results_df)

# Calculate growth rates
println("\nGrowth Rates:")
gdp_values = [GDP_MACRO[y] for y in year_all]
c_values = [value(model[:C][y]) for y in year_all]
k_values = [value(model[:K][y]) for y in year_all]

for i in 2:length(year_all)
    gdp_growth = (gdp_values[i] / gdp_values[i-1])^(1 / 10) - 1
    c_growth = (c_values[i] / c_values[i-1])^(1 / 10) - 1
    k_growth = (k_values[i] / k_values[i-1])^(1 / 10) - 1
    
    println("  $(year_all[i-1])-$(year_all[i]): GDP=$(round(gdp_growth * 100, digits=1))%, C=$(round(c_growth * 100, digits=1))%, K=$(round(k_growth * 100, digits=1))%")
end

# Key ratios
println("\nKey Ratios:")
for (i, y) in enumerate(year_all)
    saving_rate = value(model[:I][y]) / gdp_values[i]
    capital_output = k_values[i] / gdp_values[i]
    energy_cost_share = value(model[:EC][y]) / gdp_values[i]
    
    println("  $y: S/Y=$(round(saving_rate * 100, digits=1))%, K/Y=$(round(capital_output, digits=1)), EC/Y=$(round(energy_cost_share * 100, digits=1))%")
end

# Display objective value
println("\nUtility: ", round(objective_value(model), digits=2))

# Also save summary with objective value
results_summary = Dict(
    "objective_UTILITY" => objective_value(model),
    "solver_status" => string(termination_status(model))
)

# Save both files (use the same DataFrame for CSV export)
CSV.write("macro_results.csv", results_df)
using JSON
open("macro_results_summary.json", "w") do f
    JSON.print(f, results_summary, 4)
end
