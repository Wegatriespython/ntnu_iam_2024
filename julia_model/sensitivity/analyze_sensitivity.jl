# analyze_sensitivity_results.jl

using CSV
using DataFrames
using Statistics
using Printf

# --- Configuration ---
const RESULTS_DIR = "sensitivity/data"
const OUTPUT_DIR = "."
const MIN_BASELINE_THRESHOLD = 0.1  # PWh - filter out low baseline activities

# Technologies that have zero costs in YAML (should be excluded from elasticity analysis)
const ZERO_COST_TECHNOLOGIES = Set([
  "other_nele", "coal_nele", "oil_nele", "gas_nele", "bio_nele", "solar_nele",  # _nele techs
  "hydro_pot", "wind_pot", "solar_pot", "bio_pot",                              # _pot techs  
  "coal_extr", "gas_extr", "oil_extr", "nuclear_fuel",                         # _extr techs
  "electricity_grid", "appliances"                                             # other zero-cost techs
])

# --- Utility Functions ---

"""
    calculate_elasticity(baseline_value, perturbed_value, perturbation_percent)

Calculate elasticity: (% change in quantity) / (% change in cost)
"""
function calculate_elasticity(baseline_value, perturbed_value, perturbation_percent)
  if perturbation_percent == 0.0 || baseline_value <= 1e-6
    return 0.0
  end

  percent_change_quantity = (perturbed_value - baseline_value) / baseline_value
  percent_change_cost = perturbation_percent / 100.0

  return percent_change_quantity / percent_change_cost
end

"""
    safe_cv(mean_val, std_val)

Calculate coefficient of variation safely, handling zero means
"""
function safe_cv(mean_val, std_val)
  return mean_val > 1e-6 ? std_val / mean_val : 0.0
end

"""
    add_elasticity_metrics!(df, baseline_col, perturbed_col, perturbation_col)

Add elasticity, percent change, and CV columns to a DataFrame
"""
function add_elasticity_metrics!(df, baseline_col, perturbed_col, perturbation_col)
  df.percent_change = (df[!, perturbed_col] .- df[!, baseline_col]) ./ df[!, baseline_col] * 100
  df.elasticity = [calculate_elasticity(b, p, pert) for (b, p, pert) in
                   zip(df[!, baseline_col], df[!, perturbed_col], df[!, perturbation_col])]
  if :std in names(df)
    df.cv = [safe_cv(m, s) for (m, s) in zip(df[!, perturbed_col], df.std)]
  end
  return df
end

# --- Core Analysis Functions ---

"""
    find_and_load_data(dir::String) -> DataFrame

Scans a directory for sensitivity result CSVs and loads them into a unified DataFrame.
"""
function find_and_load_data(dir::String)
  println("[ Info: Starting comprehensive sensitivity analysis")
  all_dfs = DataFrame[]

  filename_regex = r"^(baseline|perturb)_?([\d\.]*)?_run_(\d+)\.csv$"

  files = readdir(dir)
  for filename in files
    m = match(filename_regex, filename)
    if m !== nothing
      run_type = m[1]
      perturb_str = m[2]
      run_num_str = m[3]

      perturb_level = run_type == "baseline" ? 0.0 :
                      isempty(perturb_str) ? error("Missing perturbation level in $filename") :
                      parse(Float64, perturb_str)
      run_number = parse(Int, run_num_str)

      try
        df = CSV.read(joinpath(dir, filename), DataFrame)
        df.run_type .= run_type
        df.perturbation_level .= perturb_level
        df.run_number .= run_number
        push!(all_dfs, df)
      catch e
        @warn "Could not load file: $filename. Error: $e"
      end
    end
  end

  if isempty(all_dfs)
    error("No valid result files found in '$dir'. Run sensitivity simulations first.")
  end

  println("[ Info: Loading $(length(all_dfs)) sensitivity result files")
  combined_df = vcat(all_dfs...)
  println("[ Info: Combined dataset contains $(nrow(combined_df)) records")
  return combined_df
end

"""
    analyze_total_cost_elasticity(df::DataFrame) -> DataFrame

Analyzes elasticity of total system cost with respect to perturbations.
"""
function analyze_total_cost_elasticity(df::DataFrame)
  println("[ Info: Analyzing total cost elasticity")
  total_costs = filter(:variable => ==("TOTAL_COST"), df)

  # Separate baseline and perturbed data
  baseline_data = filter(:run_type => ==("baseline"), total_costs)
  perturbed_data = filter(row -> startswith(row.run_type, "perturb"), total_costs)

  baseline_mean = mean(baseline_data.value)

  # Aggregate perturbed results by perturbation level
  cost_summary = combine(groupby(perturbed_data, :perturbation_level),
    :value => mean => :mean_cost,
    :value => std => :std_cost,
    :value => minimum => :min_cost,
    :value => maximum => :max_cost,
    nrow => :n_runs
  )

  # Add baseline row
  baseline_row = DataFrame(
    perturbation_level=0.0, mean_cost=baseline_mean, std_cost=0.0,
    min_cost=baseline_mean, max_cost=baseline_mean, n_runs=1
  )
  cost_summary = vcat(baseline_row, cost_summary)

  # Add elasticity metrics
  cost_summary.cv_cost = [safe_cv(m, s) for (m, s) in zip(cost_summary.mean_cost, cost_summary.std_cost)]
  cost_summary.percent_change = (cost_summary.mean_cost .- baseline_mean) ./ baseline_mean * 100
  cost_summary.elasticity = [p == 0.0 ? 0.0 : calculate_elasticity(baseline_mean, c, p) for (c, p) in
                             zip(cost_summary.mean_cost, cost_summary.perturbation_level)]

  return sort(cost_summary, :perturbation_level)
end

"""
    analyze_technology_elasticity(df::DataFrame) -> (DataFrame, DataFrame)

Analyzes technology activity elasticity, filtering out low-baseline effects.
Returns detailed analysis and summary DataFrames.
"""
function analyze_technology_elasticity(df::DataFrame)
  println("[ Info: Analyzing technology elasticity")
  act_data = filter(:variable => ==("ACT"), df)
  println("[ Info: Found $(length(unique(act_data.technology))) technologies across $(length(unique(act_data.year))) years")

  # Get baseline data and filter out low activities
  baseline_data = filter(:run_type => ==("baseline"), act_data)
  baseline_summary = combine(groupby(baseline_data, [:technology, :year]),
    :value => mean => :baseline_mean)

  # Filter for meaningful baseline activities
  baseline_summary = filter(:baseline_mean => >(MIN_BASELINE_THRESHOLD), baseline_summary)

  # Filter out technologies with zero costs in YAML
  baseline_summary = filter(:technology => tech -> !(tech in ZERO_COST_TECHNOLOGIES), baseline_summary)
  println("[ Info: Filtering to $(nrow(baseline_summary)) technology-year combinations with activity > $(MIN_BASELINE_THRESHOLD) PWh and non-zero costs")

  # Get perturbed data
  perturbed_data = filter(row -> startswith(row.run_type, "perturb"), act_data)
  perturbed_summary = combine(groupby(perturbed_data, [:technology, :year, :perturbation_level]),
    :value => mean => :perturbed_mean,
    :value => std => :perturbed_std,
    nrow => :n_runs
  )

  # Join and calculate elasticity
  tech_analysis = innerjoin(perturbed_summary, baseline_summary, on=[:technology, :year])

  # Rename column for helper function
  rename!(tech_analysis, :perturbed_std => :std)

  # Add elasticity metrics using helper function
  tech_analysis = add_elasticity_metrics!(tech_analysis, :baseline_mean, :perturbed_mean, :perturbation_level)

  # Create summary by technology and perturbation level
  tech_summary = combine(groupby(tech_analysis, [:technology, :perturbation_level]),
    :elasticity => (x -> maximum(abs.(x))) => :max_abs_elasticity,
    :elasticity => (x -> mean(abs.(x))) => :mean_abs_elasticity,
    :percent_change => (x -> maximum(abs.(x))) => :max_abs_percent_change,
    (:cv in names(tech_analysis) ? :cv : :std) => maximum => :max_cv,
    :baseline_mean => mean => :avg_baseline_activity,
    :year => (x -> length(unique(x))) => :n_years
  )

  return tech_analysis, sort(tech_summary, :max_abs_elasticity, rev=true)
end

"""
    analyze_timeseries_elasticity(df::DataFrame) -> DataFrame

Analyzes elasticity of time-series variables (costs, emissions).
"""
function analyze_timeseries_elasticity(df::DataFrame)
  println("[ Info: Analyzing time series elasticity")
  ts_data = filter(:variable => v -> v in ["COST_ANNUAL", "EMISS"], df)

  # Separate baseline and perturbed data
  baseline_data = filter(:run_type => ==("baseline"), ts_data)
  baseline_summary = combine(groupby(baseline_data, [:variable, :year]),
    :value => mean => :baseline_mean)

  perturbed_data = filter(row -> startswith(row.run_type, "perturb"), ts_data)
  perturbed_summary = combine(groupby(perturbed_data, [:variable, :year, :perturbation_level]),
    :value => mean => :perturbed_mean,
    :value => std => :perturbed_std,
    nrow => :n_runs
  )

  # Join and add elasticity metrics
  ts_analysis = innerjoin(perturbed_summary, baseline_summary, on=[:variable, :year])
  rename!(ts_analysis, :perturbed_std => :std)
  ts_analysis = add_elasticity_metrics!(ts_analysis, :baseline_mean, :perturbed_mean, :perturbation_level)

  return ts_analysis
end

"""
    generate_elasticity_report(results::Dict)

Generates a comprehensive elasticity-focused report.
"""
function generate_elasticity_report(results::Dict)

  header("TOTAL COST ELASTICITY")
  cost_summary = results[:cost_summary]
  display(cost_summary)

  header("TECHNOLOGY ELASTICITY SUMMARY")
  tech_summary = results[:tech_summary]
  println("Top 10 most elastic technologies:")
  display(first(tech_summary, min(10, nrow(tech_summary))))

  header("ELASTICITY INTERPRETATION")

  # Extract key metrics
  baseline_cost = filter(:perturbation_level => ==(0.0), cost_summary).mean_cost[1]
  perturbed_costs = filter(:perturbation_level => >(0.0), cost_summary)

  max_cost_elasticity = isempty(perturbed_costs) ? 0.0 : maximum(abs.(perturbed_costs.elasticity))
  max_tech_elasticity = isempty(tech_summary) ? 0.0 : maximum(tech_summary.max_abs_elasticity)
  max_cost_cv = isempty(perturbed_costs) ? 0.0 : maximum(perturbed_costs.cv_cost)

  # Report key findings
  println("Baseline total cost: $(round(baseline_cost, digits=2)) billion USD")
  @printf "Maximum cost elasticity: %.3f\n" max_cost_elasticity
  @printf "Maximum technology elasticity: %.3f\n" max_tech_elasticity
  @printf "Maximum cost CV: %.2f%%\n" max_cost_cv * 100


  # Highlight highly elastic technologies
  highly_elastic = filter(:max_abs_elasticity => >(1.0), tech_summary)
  if !isempty(highly_elastic)
    println("\n⚡ HIGHLY ELASTIC TECHNOLOGIES (|elasticity| > 1.0):")
    for row in eachrow(first(highly_elastic, 5))
      @printf "   %-15s: elasticity %.2f (baseline: %.1f PWh)\n" row.technology row.max_abs_elasticity row.avg_baseline_activity
    end
  end
end

"""
    save_analysis_results(results::Dict, output_dir::String)

Saves all analysis results to CSV files.
"""
function save_analysis_results(results::Dict, output_dir::String)
  println("\n[ Info: Saving analysis results to $output_dir")

  # Ensure output directory exists
  ispath(output_dir) || mkpath(output_dir)

  for (name, df) in results
    filename = "elasticity_analysis_$(name).csv"
    filepath = joinpath(output_dir, filename)
    CSV.write(filepath, df)
    println("  ✓ Saved $(filename) ($(nrow(df)) rows)")
  end
end

# --- Main Execution ---

"""
    main()

Main function orchestrating the complete elasticity analysis workflow.
"""
function main()
  # Load data
  all_data = find_and_load_data(RESULTS_DIR)

  # Perform analyses
  cost_summary = analyze_total_cost_elasticity(all_data)
  tech_analysis, tech_summary = analyze_technology_elasticity(all_data)
  ts_analysis = analyze_timeseries_elasticity(all_data)

  # Compile results
  results = Dict(
    :cost_summary => cost_summary,
    :technology_detailed => tech_analysis,
    :tech_summary => tech_summary,
    :timeseries_analysis => ts_analysis
  )

  # Generate report and save results
  generate_elasticity_report(results)
  save_analysis_results(results, OUTPUT_DIR)
end

# Execute if run as script
if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
