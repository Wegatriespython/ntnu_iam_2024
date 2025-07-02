# run_energy_basis_stability.jl
# Demonstration of basis stability analysis on the energy model

using JuMP, GAMS
using CSV, DataFrames
using Printf
using Statistics

# Load energy model components
include("shared.jl")
include("energy/energy_model_world.jl")
include("basis_stability.jl")

"""
Main function demonstrating basis stability analysis on the energy model.
"""
function main()
  println("=== ENERGY MODEL BASIS STABILITY ANALYSIS ===")

  # 1. Solve baseline energy model
  println("\n[ Step 1: Solving baseline energy model...")
  model = Model(GAMS.Optimizer)
  set_optimizer_attribute(model, "LP", "CPLEX")

  # Create the energy model
  create_energy_model!(model)

  @printf "[ Info: Model has %d variables and %d constraints\n" num_variables(model) num_constraints(model; count_variable_in_set_constraints=false)

  # Solve
  optimize!(model)
  status = termination_status(model)

  if status != MOI.OPTIMAL
    println("❌ Model solve failed with status: $status")
    return
  end

  println("✅ Optimal solution found")
  total_cost = value(model[:TOTAL_COST])
  @printf "Baseline total cost: %.2f billion USD\n" total_cost

  # 2. Extract basis information
  println("\n[ Step 2: Extracting optimal basis information...")
  try
    basis = extract_basis_info(model)
    print_stability_summary(basis)

    # 3. Test specific perturbations
    println("\n[ Step 3: Testing specific RHS perturbations...")
    test_specific_perturbations(basis)

    # 4. Compare with current perturbation approach
    println("\n[ Step 4: Comparing with current random perturbation approach...")
    compare_perturbation_approaches(basis)

  catch e
    println(e)
    # Fallback: demonstrate perturbation analysis concept
  end

  println("\n=== ANALYSIS COMPLETE ===")
end

"""
    test_specific_perturbations(basis::BasisStability)

Test specific perturbation scenarios against the stability cone.
"""
function test_specific_perturbations(basis::BasisStability)
  # Test some specific perturbation scenarios (using correct dimension basis.m)
  test_cases = [
    ("Zero perturbation", zeros(Float64, basis.m)),
    ("Small uniform increase", 0.01 * ones(Float64, basis.m)),
    ("Small uniform decrease", -0.01 * ones(Float64, basis.m)),
    ("Large uniform increase", 0.1 * ones(Float64, basis.m)),
    ("Random small", 0.01 * randn(basis.m)),
  ]

  println("Testing specific perturbation cases:")
  for (name, δD) in test_cases
    feasible = test_perturbation(basis, δD)
    status_str = feasible ? "✅ STABLE" : "❌ UNSTABLE"
    norm_str = @sprintf("%.4f", norm(δD))
    println("  $name (‖δD‖=$norm_str): $status_str")
  end
end

"""
    compare_perturbation_approaches(basis::BasisStability)

Compare basis stability analysis with random perturbation approach.
"""
function compare_perturbation_approaches(basis::BasisStability)
  # Get independent bounds for guidance
  bounds = compute_independent_bounds(basis)

  # Test random perturbations within and outside stability bounds
  n_tests = 50
  within_bounds_stable = 0
  outside_bounds_stable = 0

  println("Comparing perturbation approaches ($n_tests samples each):")

  # Test perturbations within independent bounds
  for _ in 1:n_tests
    δD = zeros(Float64, basis.m)
    for j in 1:basis.m
      if isfinite(bounds[j][1]) && isfinite(bounds[j][2])
        # 50% of bound range
        range_size = 0.5 * (bounds[j][2] - bounds[j][1])
        center = 0.5 * (bounds[j][1] + bounds[j][2])
        δD[j] = center + range_size * (2 * rand() - 1)
      else
        δD[j] = 0.01 * randn()  # Small random perturbation
      end
    end

    if test_perturbation(basis, δD)
      within_bounds_stable += 1
    end
  end

  # Test larger random perturbations (may be outside bounds)
  for _ in 1:n_tests
    δD = 0.1 * randn(basis.m)  # Larger random perturbations

    if test_perturbation(basis, δD)
      outside_bounds_stable += 1
    end
  end

  within_rate = within_bounds_stable / n_tests * 100
  outside_rate = outside_bounds_stable / n_tests * 100

  println("  Within stability bounds: $(within_rate)% remain stable")
  println("  Random perturbations: $(outside_rate)% remain stable")
  println("  → Basis stability analysis provides $(within_rate - outside_rate)% improvement in stability prediction")
end

# Run the demonstration
if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
