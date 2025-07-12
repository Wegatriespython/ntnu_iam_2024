using Distributed

# Add workers for parallel computation
nworkers = min(Sys.CPU_THREADS - 1, 8)
if nprocs() < nworkers + 1
  addprocs(nworkers - nprocs() + 1)
end

# Load necessary modules on main process
using DataFrames
using CSV
using YAML
using JSON
using Dates

# Load necessary modules on all workers
@everywhere begin
  using Pkg
  Pkg.activate(".")
  using JuMP
  using HiGHS
  using DataFrames
  using YAML
  using JSON
  using Statistics
  using Dates
end

# Load model code on workers
@everywhere begin
  include("shared.jl")
  include("energy/energy_model_world.jl")
  include("energy/energy_model_reduced.jl")
end

# Function to extract all decision variables from full model
@everywhere function extract_full_model_results(model)
  results = Dict{String,Any}()

  # Get model dimensions from globals
  years = year_all
  techs = technology

  # Store objective
  results["objective"] = objective_value(model)
  results["solve_time"] = solve_time(model)
  results["status"] = string(termination_status(model))

  # Extract all decision variables
  # ACT - Activity levels [tech, year]
  results["ACT"] = Dict(
    "$(t)_$(y)" => value(model[:ACT][t, y])
    for t in techs, y in years
  )

  # CAP_NEW - New capacity [tech, year]
  results["CAP_NEW"] = Dict(
    "$(t)_$(y)" => value(model[:CAP_NEW][t, y])
    for t in techs, y in years
  )


  # COST_ANNUAL - Annual costs [year]
  results["COST_ANNUAL"] = Dict(
    string(y) => value(model[:COST_ANNUAL][y])
    for y in years
  )

  # EMISS - Emissions [year]
  results["EMISS"] = Dict(
    string(y) => value(model[:EMISS][y])
    for y in years
  )

  # CUM_EMISS - Cumulative emissions
  results["CUM_EMISS"] = value(model[:CUM_EMISS])

  # TOTAL_COST
  results["TOTAL_COST"] = value(model[:TOTAL_COST])

  return results
end

# Function to extract all decision variables from reduced model
@everywhere function extract_reduced_model_results(model)
  results = Dict{String,Any}()

  # Get model dimensions
  years = [2020, 2030, 2040, 2050, 2060, 2070, 2080]
  groups = collect(keys(TECH_GROUPS))

  # Store objective
  results["objective"] = objective_value(model)
  results["solve_time"] = solve_time(model)
  results["status"] = string(termination_status(model))

  # Extract all decision variables
  # ACT_R - Activity levels [group, year]
  results["ACT_R"] = Dict(
    "$(g)_$(y)" => value(model[:ACT_R][g, y])
    for g in groups, y in years
  )

  # CAP_NEW_R - New capacity [group, year]
  results["CAP_NEW_R"] = Dict(
    "$(g)_$(y)" => value(model[:CAP_NEW_R][g, y])
    for g in groups, y in years
  )

  # COST_ANNUAL_R - Annual costs [year]
  results["COST_ANNUAL_R"] = Dict(
    string(y) => value(model[:COST_ANNUAL_R][y])
    for y in years
  )


  # TOTAL_COST_R
  results["TOTAL_COST_R"] = value(model[:TOTAL_COST_R])

  return results
end

# Helper function to run full model
@everywhere function run_full_model(config_path::String)
  try
    full_model = Model(HiGHS.Optimizer)
    set_optimizer_attribute(full_model, "log_to_console", false)
    set_optimizer_attribute(full_model, "time_limit", 60.0)

    Base.invokelatest(create_energy_model!, full_model, config_file=config_path)
    optimize!(full_model)

    if termination_status(full_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
      return extract_full_model_results(full_model)
    else
      return Dict("status" => string(termination_status(full_model)))
    end
  catch e
    return Dict("status" => "error", "error" => string(e))
  end
end

# Helper function to run reduced model
@everywhere function run_reduced_model(config_path::String)
  try
    reduced_model = Model(HiGHS.Optimizer)
    set_optimizer_attribute(reduced_model, "log_to_console", false)
    set_optimizer_attribute(reduced_model, "time_limit", 60.0)

    # Set up model extensions
    reduced_model.ext[:year_all] = [2020, 2030, 2040, 2050, 2060, 2070, 2080]
    reduced_model.ext[:period_length] = 10
    reduced_model.ext[:discount_rate] = 0.05

    # Create reduced model with custom parameters (unscaled)
    Base.invokelatest(create_reduced_energy_model!, reduced_model, config_file=config_path)

    optimize!(reduced_model)

    if termination_status(reduced_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
      return extract_reduced_model_results(reduced_model)
    else
      return Dict("status" => string(termination_status(reduced_model)))
    end
  catch e
    return Dict("status" => "error", "error" => string(e))
  end
end

# Function to run both models on a single configuration
@everywhere function run_single_config(config_file::String, sample_id::Int, demand_scenario::AbstractString)
  results = Dict{String,Any}()

  # Metadata
  results["sample_id"] = sample_id
  results["demand_scenario"] = demand_scenario
  results["config_file"] = config_file
  results["timestamp"] = string(now())

  # Load configuration
  config_path = joinpath("config_samples", config_file)
  custom_params = YAML.load_file(config_path)

  # Store key parameter values
  results["beta"] = custom_params["model"]["beta"]

  println("Running models concurrently for sample $sample_id, demand $demand_scenario...")

  # Run both models concurrently
  full_future = @spawn run_full_model(config_path)
  reduced_future = @spawn run_reduced_model(config_path)

  # Collect results
  results["full"] = fetch(full_future)
  results["reduced"] = fetch(reduced_future)

  return results
end

# Main execution
function main()
  # Get list of config files
  config_files = readdir("config_samples")
  filter!(f -> endswith(f, ".yaml") && startswith(f, "config_sample_"), config_files)

  println("Found $(length(config_files)) configuration files")
  println("Using $nworkers workers for parallel execution")

  # Prepare tasks
  tasks = []
  for file in config_files
    # Parse filename
    m = match(r"config_sample_(\d+)_demand_(\w+)\.yaml", file)
    if !isnothing(m)
      sample_id = parse(Int, m.captures[1])
      demand_scenario = m.captures[2]
      push!(tasks, (file, sample_id, demand_scenario))
    end
  end

  # Run in parallel
  println("\nGenerating bias dataset...")
  println("Processing $(length(tasks)) configurations...")

  all_results = @distributed (vcat) for (file, sample_id, demand) in tasks
    result = run_single_config(file, sample_id, demand)
    println("Completed: sample $sample_id, demand $demand")
    [result]
  end

  # Save full results as JSON
  timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
  json_file = "bias_dataset_full_$(timestamp).json"
  open(json_file, "w") do f
    JSON.print(f, all_results, 4)
  end

  # Create summary CSV for quick analysis
  summary_rows = []
  for result in all_results
    row = Dict{String,Any}()
    row["sample_id"] = result["sample_id"]
    row["demand_scenario"] = result["demand_scenario"]
    row["beta"] = result["beta"]

    # Full model summary
    if haskey(result["full"], "objective")
      row["full_total_cost"] = result["full"]["TOTAL_COST"]
      row["full_cum_emissions"] = result["full"]["CUM_EMISS"]
      row["full_solve_time"] = result["full"]["solve_time"]
      row["full_status"] = result["full"]["status"]

      # Sample annual costs by year
      for y in [2020, 2050, 2080]
        row["full_cost_$y"] = result["full"]["COST_ANNUAL"][string(y)]
      end
    else
      row["full_status"] = get(result["full"], "status", "failed")
    end

    # Reduced model summary
    if haskey(result["reduced"], "objective")
      row["reduced_total_cost"] = result["reduced"]["TOTAL_COST_R"]
      row["reduced_solve_time"] = result["reduced"]["solve_time"]
      row["reduced_status"] = result["reduced"]["status"]

      # Sample annual costs by year
      for y in [2020, 2050, 2080]
        row["reduced_cost_$y"] = result["reduced"]["COST_ANNUAL_R"][string(y)]
      end
    else
      row["reduced_status"] = get(result["reduced"], "status", "failed")
    end

    push!(summary_rows, row)
  end

  summary_df = DataFrame(summary_rows)
  csv_file = "bias_dataset_summary_$(timestamp).csv"
  CSV.write(csv_file, summary_df)

  println("\n✓ Completed $(length(all_results)) runs")
  println("✓ Full results saved to: $json_file")
  println("✓ Summary saved to: $csv_file")

  # Quick statistics
  successful = count(r ->
      haskey(r["full"], "objective") && haskey(r["reduced"], "objective"),
    all_results
  )
  println("\nSuccessful runs: $successful / $(length(all_results))")
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
