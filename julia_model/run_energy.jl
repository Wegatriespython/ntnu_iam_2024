# run_energy.jl - Standalone energy model using corrected GAMS-matching parameters

using JuMP, GAMS

println("Starting standalone energy model (GAMS-equivalent)...")

# Load all definitions from energy_model_world.jl
include("shared.jl")
include("energy/energy_model_world.jl")

# ------------------------------------------------------------------------------
# Create and solve model
# ------------------------------------------------------------------------------

model = Model(GAMS.Optimizer)

set_optimizer_attribute(model, "LP", "CPLEX")

# Create the energy model with all constraints and objective
create_energy_model!(model)

# ------------------------------------------------------------------------------
# Solve model
# ------------------------------------------------------------------------------

println("Model variables: $(num_variables(model))")
println("Model constraints: $(num_constraints(model; count_variable_in_set_constraints=false))")

optimize!(model)

# ------------------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------------------

println("\nTermination status: $(termination_status(model))")

if termination_status(model) ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
  total_cost = value(model[:TOTAL_COST])
  println("✅ Optimal solution found")
  println("Total system cost: $(round(total_cost, digits=2)) billion USD")

  println("\nKey results for 2020:")
  for t in ["coal_ppl", "gas_ppl", "oil_ppl", "bio_ppl", "hydro_ppl", "wind_ppl", "solar_PV_ppl", "nuclear_ppl"]
    if t in technology
      act = value(model[:ACT][t, 2020])
      if act > 0.01
        println("  ACT[$t] = $(round(act, digits=2)) PWh")
      end
    end
  end

  println("\nNon-electric technologies 2020:")
  for t in ["coal_nele", "gas_nele", "oil_nele", "bio_nele", "solar_nele", "other_nele"]
    if t in technology
      act = value(model[:ACT][t, 2020])
      if act > 0.01
        println("  ACT[$t] = $(round(act, digits=2)) PWh")
      end
    end
  end

  println("\nAnnual costs (billion USD):")
  for y in year_all
    annual_cost = value(model[:COST_ANNUAL][y])
    println("  $y: $(round(annual_cost, digits=1)) billion USD")
  end

  println("\nEmissions (Mt CO2):")
  for y in year_all
    emiss = value(model[:EMISS][y])
    println("  $y: $(round(emiss, digits=1)) Mt CO2")
  end


else
  println("❌ Model failed: $(termination_status(model))")
end

println("\n✓ Energy model completed")
