# run_energy_reduced.jl - Standalone runner for reduced form energy model

using JuMP
using HiGHS
using DataFrames
using CSV
using JSON

# Include the reduced energy model
include("energy/energy_model_reduced.jl")

function run_reduced_energy_model()
  println("=== Running Reduced Form Energy Model ===")

  # Create model
  model = Model(HiGHS.Optimizer)
  set_optimizer_attribute(model, "log_to_console", true)
  set_optimizer_attribute(model, "time_limit", 300.0)

  # Set up model extensions for shared parameters
  model.ext[:year_all] = [2020, 2030, 2040, 2050, 2060, 2070, 2080]
  model.ext[:period_length] = 10
  model.ext[:discount_rate] = 0.05

  # Create the reduced energy model
  create_reduced_energy_model!(model)

  # Solve
  println("\nSolving reduced energy model...")
  optimize!(model)

  # Check solution status
  status = termination_status(model)
  println("\nSolution status: ", status)

  if status == MOI.LOCALLY_SOLVED || status == MOI.OPTIMAL
    # Extract results
    years = model.ext[:year_all]

    # Activity levels by technology group
    println("\n=== Activity Levels (GWh) ===")
    activity_df = DataFrame(Year=years)
    groups = collect(keys(TECH_GROUPS))
    for group in groups
      activity_df[!, Symbol(group)] = [value(model[:ACT_R][group, y]) for y in years]
    end
    println(activity_df)

    # Capacity by technology group
    println("\n=== Capacity (MW) ===")
    capacity_df = DataFrame(Year=years)
    for group in groups
      capacity_df[!, Symbol(group)] = [value(model[:CAP_NEW_R][group, y]) for y in years]
    end
    println(capacity_df)

    # Cost analysis
    println("\n=== Cost Analysis ===")
    annual_costs = [value(model[:COST_ANNUAL_R][y]) for y in years]
    cost_df = DataFrame(
      Year=years,
      Annual_Cost_MUSD=annual_costs,
      Discounted_Cost=[annual_costs[i] * 10 * (1.05)^(-10 * (i - 1)) for i in 1:length(years)]
    )
    println(cost_df)
    println("\nTotal Discounted Cost: ", round(value(model[:TOTAL_COST_R]), digits=2), " Million USD")

    # Technology shares
    println("\n=== Electricity Generation Shares (%) ===")
    shares_df = DataFrame(Year=years)

    elec_groups = ["fossil_power", "renewable_power", "nuclear_power", "biomass_power"]

    for group in elec_groups
      if group in groups
        tech_shares = Float64[]
        for y in years
          total_elec = sum(value(model[:ACT_R][g, y]) *
                           get(group_params[g]["output_coeffs"], ("electricity", "secondary"), 0.0)
                           for g in elec_groups if g in groups)

          if total_elec > 0
            share = value(model[:ACT_R][group, y]) *
                    get(group_params[group]["output_coeffs"], ("electricity", "secondary"), 0.0) / total_elec * 100
          else
            share = 0.0
          end
          push!(tech_shares, share)
        end
        shares_df[!, Symbol(group)] = tech_shares
      end
    end
    println(shares_df)

    # Save results
    results = Dict(
      "status" => string(status),
      "total_cost" => value(model[:TOTAL_COST_R]),
      "activity" => Dict(string(group) => Dict(string(y) => value(model[:ACT_R][group, y])
                                               for y in years) for group in groups),
      "capacity" => Dict(string(group) => Dict(string(y) => value(model[:CAP_NEW_R][group, y])
                                               for y in years) for group in groups),
      "annual_costs" => Dict(string(y) => value(model[:COST_ANNUAL_R][y]) for y in years)
    )

    open("reduced_energy_results.json", "w") do f
      JSON.print(f, results, 4)
    end
    println("\nResults saved to reduced_energy_results.json")

    # Model statistics
    println("\n=== Model Statistics ===")
    println("Number of variables: ", num_variables(model))
    println("Number of constraints: ", num_constraints(model; count_variable_in_set_constraints=false))
    println("Solve time: ", round(solve_time(model), digits=3), " seconds")

  else
    println("Model did not solve successfully!")
    println("Termination status: ", status)
    println("Primal status: ", primal_status(model))
    println("Dual status: ", dual_status(model))
  end

  return model
end

# Run the model
if abspath(PROGRAM_FILE) == @__FILE__
  model = run_reduced_energy_model()
end
