# run_energy_refactored.jl
# A robust and modular script to run the energy model with optional perturbations.

using JuMP, GAMS
using CSV, DataFrames
using ArgParse         # For robust command-line parsing
using Random           # For perturbations
using Printf           # For formatted output

# --- Function Definitions ---

"""
    parse_commandline_args() -> Dict

Sets up and parses command-line arguments using ArgParse.jl.
Returns a dictionary of the parsed arguments.
"""
function parse_commandline_args()
  s = ArgParseSettings(
    description="Run the standalone energy model.",
    prog="run_energy.jl"
  )
  @add_arg_table! s begin
    "--perturb"
    help = "Apply a random cost perturbation. Specify the percentage bound."
    arg_type = Float64
    default = 0.0
    "--output"
    help = "Filename to save the results to a CSV file."
    arg_type = String
    default = ""
    "--seed"
    help = "Random seed for reproducible perturbations."
    arg_type = Int
    default = -1 # A value of -1 means no seed will be set
  end
  return parse_args(s)
end

"""
    apply_perturbations!(cost_data::Vector{<:Dict}, bound::Float64)

Applies a random perturbation within [-bound, +bound] to all values in the
provided cost dictionaries. The dictionaries are modified in-place.
"""
function apply_perturbations!(cost_data::Vector{<:Dict}, percent_bound::Float64)
  if percent_bound <= 0.0
    return
  end

  bound = percent_bound / 100.0
  println("\n[ Info: Applying random perturbations within ±$(percent_bound)% bounds")

  for cost_dict in cost_data
    for (key, value) in cost_dict
      random_factor = 1 + (2 * rand() - 1) * bound
      # Apply perturbation and ensure value remains positive
      cost_dict[key] = max(0.001, value * random_factor)
    end
  end
  println("[ Success: Random cost perturbations applied")
end

"""
    save_perturbed_costs(cost_data::Vector{<:Dict}, filename::String)

Save the perturbed cost data to a CSV file for debugging purposes.
"""
function save_perturbed_costs(cost_data::Vector{<:Dict}, filename::String)
  println("[ Info: Saving perturbed costs to '$filename'")
  
  cost_names = ["inv_cost", "fom", "vom"]
  all_costs = DataFrame()
  
  for (i, cost_dict) in enumerate(cost_data)
    for (tech, cost) in cost_dict
      temp_df = DataFrame(
        cost_type = cost_names[i],
        technology = tech,
        cost = cost
      )
      append!(all_costs, temp_df)
    end
  end
  
  CSV.write(filename, all_costs)
  println("[ Success: Perturbed costs saved ($(nrow(all_costs)) records)")
end

"""
    extract_results(model::Model) -> DataFrame

Extracts all key variables (ACT, CAP_NEW, COST_ANNUAL, etc.) from a solved JuMP model
and returns them as a single, tidy DataFrame.
"""
function extract_results(model::Model)
  results_df = DataFrame()

  # Helper function to extract a 2D variable and append it to the main DataFrame
  function append_2d_variable!(df, var_name, unit)
    var = model[var_name]
    
    # Get the axes
    techs = var.axes[1]
    years = var.axes[2]
    
    # Create result vectors
    tech_vec = String[]
    year_vec = Int[]
    value_vec = Float64[]
    
    # Extract non-zero values
    for tech in techs
      for year in years
        val = JuMP.value(var[tech, year])
        if val > 1e-6
          push!(tech_vec, tech)
          push!(year_vec, year)
          push!(value_vec, val)
        end
      end
    end
    
    # Create temporary DataFrame
    temp_df = DataFrame(
      variable=string(var_name),
      technology=tech_vec,
      year=year_vec,
      value=value_vec,
      unit=unit
    )
    
    # Append to main DataFrame
    append!(df, temp_df)
  end

  # Extract multi-dimensional variables
  append_2d_variable!(results_df, :ACT, "PWh")
  append_2d_variable!(results_df, :CAP_NEW, "GW")

  # Extract annual variables (single dimension)
  for var_info in [
    (name=:COST_ANNUAL, unit="billion_USD"),
    (name=:EMISS, unit="Mt_CO2")
  ]
    var = model[var_info.name]
    years = var.axes[1]
    
    # Manual extraction for 1D arrays
    year_vec = Int[]
    value_vec = Float64[]
    
    for year in years
      push!(year_vec, year)
      push!(value_vec, JuMP.value(var[year]))
    end
    
    temp_df = DataFrame(
      variable=string(var_info.name),
      technology="ALL",
      year=year_vec,
      value=value_vec,
      unit=var_info.unit
    )
    append!(results_df, temp_df)
  end

  # Extract scalar variables
  for var_info in [
    (name=:TOTAL_COST, unit="billion_USD"),
    (name=:CUM_EMISS, unit="Mt_CO2")
  ]
    temp_df = DataFrame(
      variable=string(var_info.name),
      technology="ALL",
      year=0,
      value=value(model[var_info.name]),
      unit=var_info.unit
    )
    append!(results_df, temp_df)
  end

  return results_df
end

"""
    main()

Main function to orchestrate the entire model run.
"""
function main()
  # 1. Parse Arguments
  args = parse_commandline_args()
  perturb_percent = args["perturb"]
  output_filename = args["output"]
  seed = args["seed"]

  println("--- Starting Energy Model Run ---")
  if perturb_percent > 0
    println("[ Config: Perturbation enabled at ±$(perturb_percent)%")
    if seed != -1
      Random.seed!(seed)
      println("[ Config: Random seed set to $seed")
    end
  end
  if !isempty(output_filename)
    println("[ Config: Results will be saved to '$output_filename'")
  end

  # 2. Load Model Data
  include("shared.jl")
  include("energy/energy_model_world.jl")

  # 3. Apply Perturbations (if requested)
  apply_perturbations!([inv_cost, fom, vom], perturb_percent)
  
  # Save perturbed input costs if output file is specified
  if !isempty(output_filename) && perturb_percent > 0.0
    costs_filename = replace(output_filename, ".csv" => "_costs.csv")
    save_perturbed_costs([inv_cost, fom, vom], costs_filename)
  end

  # After potential perturbation, must recalculate derived cost parameters
  cost_capacity, cost_activity, cost = Base.invokelatest(calculate_costs)

  # 4. Create and Solve the JuMP Model
  println("\n[ Info: Building the optimization model...")
  model = Model(GAMS.Optimizer)
  set_optimizer_attribute(model, "LP", "CPLEX")
  Base.invokelatest(create_energy_model!, model)

  @printf "[ Info: Model has %d variables and %d constraints\n" num_variables(model) num_constraints(model; count_variable_in_set_constraints=false)
  println("[ Info: Starting solver...")
  optimize!(model)

  # 5. Report and Save Results
  println("\n--- Results ---")
  status = termination_status(model)
  println("Termination status: $status")

  if status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
    println("✅ Optimal solution found")
    total_cost = value(model[:TOTAL_COST])
    @printf "Total system cost: %.2f billion USD\n" total_cost

    if !isempty(output_filename)
      println("\n[ Info: Extracting results for saving...")
      results_df = extract_results(model)
      CSV.write(output_filename, results_df)
      println("[ Success: Results saved to '$output_filename' ($(nrow(results_df)) records)")
    else
      println("\n[ Info: Skipping result saving (no --output specified).")
    end
  else
    println("❌ Model solve failed.")
  end

  println("\n--- Energy Model Completed ---")
end


# --- Script Execution ---
if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
