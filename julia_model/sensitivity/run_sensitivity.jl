# run_sensitivity.jl
# Script to generate sensitivity analysis data by running multiple perturbed energy model runs

using ArgParse
using Printf
using Random

"""
    parse_commandline_args() -> Dict

Sets up and parses command-line arguments for the sensitivity analysis.
Returns a dictionary of the parsed arguments.
"""
function parse_commandline_args()
    s = ArgParseSettings(
        description="Run sensitivity analysis by generating multiple perturbed energy model runs.",
        prog="run_sensitivity.jl"
    )
    @add_arg_table! s begin
        "--num-runs", "-n"
        help = "Number of runs per perturbation level"
        arg_type = Int
        default = 10
        "--perturb-levels", "-p"
        help = "Comma-separated list of perturbation levels (percentages)"
        arg_type = String
        default = "1.0,5.0"
        "--data-dir", "-d"
        help = "Directory to save the generated data files"
        arg_type = String
        default = "sensitivity/data"
        "--baseline-runs", "-b"
        help = "Number of baseline runs (no perturbation)"
        arg_type = Int
        default = 1
        "--seed"
        help = "Base random seed for reproducible results"
        arg_type = Int
        default = 42
    end
    return parse_args(s)
end

"""
    run_energy_model(perturb_level::Float64, run_number::Int, output_file::String, seed::Int)

Runs the energy model with specified perturbation level and saves results to output file.
"""
function run_energy_model(perturb_level::Float64, run_number::Int, output_file::String, seed::Int)
    cmd_args = []
    
    # Add perturbation if not baseline
    if perturb_level > 0.0
        push!(cmd_args, "--perturb", string(perturb_level))
    end
    
    # Add output file
    push!(cmd_args, "--output", output_file)
    
    # Add seed
    push!(cmd_args, "--seed", string(seed))
    
    # Run the energy model
    julia_cmd = `julia --project=@. run_energy.jl $cmd_args`
    
    println("[ Info: Running: $julia_cmd")
    
    try
        run(julia_cmd)
        println("[ Success: Completed run $(run_number) with $(perturb_level)% perturbation")
        return true
    catch e
        println("[ Error: Failed run $(run_number) with $(perturb_level)% perturbation: $e")
        return false
    end
end

"""
    generate_filename(run_type::String, perturb_level::Float64, run_number::Int, data_dir::String) -> String

Generates standardized filename for sensitivity analysis results.
"""
function generate_filename(run_type::String, perturb_level::Float64, run_number::Int, data_dir::String)
    if run_type == "baseline"
        filename = "baseline_run_$(run_number).csv"
    else
        filename = "perturb_$(perturb_level)_run_$(run_number).csv"
    end
    return joinpath(data_dir, filename)
end

"""
    main()

Main function to orchestrate the sensitivity analysis data generation.
"""
function main()
    # Parse arguments
    args = parse_commandline_args()
    num_runs = args["num-runs"]
    perturb_levels_str = args["perturb-levels"]
    data_dir = args["data-dir"]
    baseline_runs = args["baseline-runs"]
    base_seed = args["seed"]
    
    # Parse perturbation levels
    perturb_levels = [parse(Float64, strip(level)) for level in split(perturb_levels_str, ',')]
    
    println("--- Starting Sensitivity Analysis Data Generation ---")
    println("[ Config: Number of runs per perturbation level: $num_runs")
    println("[ Config: Perturbation levels: $perturb_levels")
    println("[ Config: Baseline runs: $baseline_runs")
    println("[ Config: Data directory: $data_dir")
    println("[ Config: Base seed: $base_seed")
    
    # Ensure data directory exists
    mkpath(data_dir)
    
    # Track results
    total_runs = 0
    successful_runs = 0
    failed_runs = 0
    
    # Run baseline cases
    println("\n=== Running Baseline Cases ===")
    for run_num in 1:baseline_runs
        seed = base_seed + run_num
        output_file = generate_filename("baseline", 0.0, run_num, data_dir)
        
        if run_energy_model(0.0, run_num, output_file, seed)
            successful_runs += 1
        else
            failed_runs += 1
        end
        total_runs += 1
    end
    
    # Run perturbed cases
    println("\n=== Running Perturbed Cases ===")
    for perturb_level in perturb_levels
        println("\n--- Perturbation Level: $(perturb_level)% ---")
        
        for run_num in 1:num_runs
            # Ensure unique seeds for each perturbation level and run
            seed = base_seed + Int(perturb_level * 1000) + run_num + 1000
            output_file = generate_filename("perturbed", perturb_level, run_num, data_dir)
            
            if run_energy_model(perturb_level, run_num, output_file, seed)
                successful_runs += 1
            else
                failed_runs += 1
            end
            total_runs += 1
        end
    end
    
    # Summary
    println("\n=== Sensitivity Analysis Summary ===")
    println("Total runs attempted: $total_runs")
    println("Successful runs: $successful_runs")
    println("Failed runs: $failed_runs")
    @printf "Success rate: %.1f%%\n" (successful_runs / total_runs * 100)
    
    if successful_runs > 0
        println("\n[ Info: Data files saved to: $data_dir")
        println("[ Info: You can now run the analysis using:")
        println("         julia --project=@. sensitivity/analyze_sensitivity.jl")
    else
        println("\n[ Error: No successful runs completed. Check model configuration.")
    end
    
    println("\n--- Sensitivity Analysis Data Generation Completed ---")
end

# Script execution
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end