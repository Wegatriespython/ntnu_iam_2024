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

# DISPLAY GDP_MACRO.L
println("\nGDP_MACRO.L:")
for y in year_all
    println("$y: ", round(GDP_MACRO[y], digits=4))
end

# Display objective value
println("\nUTILITY objective: ", objective_value(model))

# execute_unload "macro_results.gdx"
using CSV, DataFrames
results_df = DataFrame(
    year = collect(year_all),
    GDP_MACRO = [GDP_MACRO[y] for y in year_all],
    I = [value(model[:I][y]) for y in year_all],
    C = [value(model[:C][y]) for y in year_all],
    EC = [value(model[:EC][y]) for y in year_all]
)

# Also save summary with objective value
results_summary = Dict(
    "objective_UTILITY" => objective_value(model),
    "solver_status" => string(termination_status(model))
)

# Save both files
CSV.write("macro_results.csv", results_df)
using JSON
open("macro_results_summary.json", "w") do f
    JSON.print(f, results_summary, 4)
end
