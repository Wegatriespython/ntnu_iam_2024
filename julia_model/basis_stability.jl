# basis_stability.jl
# Proper implementation of post-optimal basis stability analysis for LP models
# Following the mathematical formulation: K_B = {δD : M*δD >= -x*_B}

using JuMP, GAMS
using LinearAlgebra
using SparseArrays
using Printf
using Statistics
import MathOptInterface as MOI

"""
    BasisStability

Structure to hold basis stability analysis results.
Contains the optimal basis B, basic matrix inverse M = A_B^(-1), 
basic solution x*_B = A_B^(-1)*D^0, and stability cone definition.
"""
struct BasisStability
    # Basic variables (m variables forming the basis)
    basic_vars::Vector{VariableRef}
    basic_indices::Vector{Int}
    
    # Non-basic variables
    nonbasic_vars::Vector{VariableRef}
    nonbasic_indices::Vector{Int}
    
    # Basic solution and matrix inverse (both m×m)
x_B::Vector{Float64}        # Basic solution: x*_B = A_B^(-1) * D^0
    M::Matrix{Float64}          # Basic matrix inverse: M = A_B^(-1)
    
    # Original RHS vector D^0 (m×1)
    D0::Vector{Float64}         # Original RHS vector D^0
    
    # Constraint mapping information
    constraint_refs::Vector{ConstraintRef}
    constraint_names::Vector{String}
    
    # Dimensions for validation
    m::Int  # Number of constraints (= number of basic variables)
    n::Int  # Total number of variables
end

"""
    extract_basis_info(model::Model)

Extract optimal basis information from a solved LP model.
Returns BasisStability struct with proper m×m basic matrix inverse.
"""
function extract_basis_info(model::Model)
    if termination_status(model) != MOI.OPTIMAL
        error("Model must be optimally solved before basis extraction")
    end
    
    println("[ Info: Extracting basis information from optimal solution...")
    
    # Get all variables and constraints
    all_vars = all_variables(model)
    n = length(all_vars)
    
    # Extract structural constraints (equality + active inequality constraints)
    constraint_refs, constraint_names, D0 = extract_structural_constraints(model)
    m = length(constraint_refs)
    
    println("   Model dimensions:")
    println("     Total variables (n): $n")
    println("     Structural constraints (m): $m")
    
    # Extract basis using MOI basis status
    basic_vars, nonbasic_vars = extract_optimal_basis(model, all_vars)
    
    # Ensure we have exactly m basic variables
    if length(basic_vars) != m
        error("Number of basic variables ($(length(basic_vars))) extracted from MOI.BASIC does not match number of structural constraints ($m). This indicates a potential issue with the solver's basis reporting or a degenerate problem. Basis stability analysis requires a square, invertible basic matrix.")
    end
    
    # Build basic indices
    basic_indices = [findfirst(v -> v == bv, all_vars) for bv in basic_vars]
    nonbasic_indices = [findfirst(v -> v == nbv, all_vars) for nbv in nonbasic_vars]
    
    # Extract basic solution
    x_B = [value(var) for var in basic_vars]
    
    # Build basic matrix A_B and compute M = A_B^(-1)
    A_B = build_basic_matrix(model, basic_vars, constraint_refs)
    M = compute_basic_matrix_inverse(A_B)
    
    println("   Basis extraction complete:")
    println("     Basic variables: $(length(basic_vars))")
    println("     Non-basic variables: $(length(nonbasic_vars))")
    println("     Basic matrix condition number: $(cond(A_B))")
    
    return BasisStability(
        basic_vars, basic_indices,
        nonbasic_vars, nonbasic_indices,
        x_B, M, D0,
        constraint_refs, constraint_names,
        m, n
    )
end

"""
    extract_structural_constraints(model::Model)

Extract the structural constraints that define the LP problem.
Returns constraint references, names, and RHS values.
"""
function extract_structural_constraints(model::Model)
    constraint_refs = Vector{ConstraintRef}()
    constraint_names = String[]
    D0 = Float64[]
    
    # Get all constraint types
    constraint_types = list_of_constraint_types(model)
    
    for (func_type, set_type) in constraint_types
        constrs = all_constraints(model, func_type, set_type)
        
        # Only include AffExpr constraints (structural constraints)
        if func_type == AffExpr
            for constr in constrs
                if set_type == MOI.EqualTo{Float64}
                    # All equality constraints are structural
                    push!(constraint_refs, constr)
                    push!(constraint_names, string(constr))
                    push!(D0, MOI.get(model, MOI.ConstraintSet(), constr).value)
                    
                elseif set_type == MOI.GreaterThan{Float64}
                    # Include active inequality constraints (non-zero dual)
                    try
                        dual_val = dual(constr)
                        if abs(dual_val) > 1e-10  # Active constraint
                            push!(constraint_refs, constr)
                            push!(constraint_names, string(constr))
                            push!(D0, MOI.get(model, MOI.ConstraintSet(), constr).lower)
                        end
                    catch
                        # If dual not available, skip or include based on heuristic
                    end
                    
                elseif set_type == MOI.LessThan{Float64}
                    # Include active inequality constraints (non-zero dual)
                    try
                        dual_val = dual(constr)
                        if abs(dual_val) > 1e-10  # Active constraint
                            push!(constraint_refs, constr)
                            push!(constraint_names, string(constr))
                            push!(D0, MOI.get(model, MOI.ConstraintSet(), constr).upper)
                        end
                    catch
                        # If dual not available, skip or include based on heuristic
                    end
                end
            end
        end
    end
    
    return constraint_refs, constraint_names, D0
end

"""
    extract_optimal_basis(model::Model, all_vars::Vector{VariableRef})

Extract basic and non-basic variables using MOI basis status.
"""
function extract_optimal_basis(model::Model, all_vars::Vector{VariableRef})
    basic_vars = VariableRef[]
    nonbasic_vars = VariableRef[]
    
    for var in all_vars
        try
            # Try MOI basis status first
            basis_status = get_attribute(var, MOI.VariableBasisStatus())
            if basis_status == MOI.BASIC
                push!(basic_vars, var)
            else
                push!(nonbasic_vars, var)
            end
        catch
            # Fallback: heuristic based on bounds
            var_value = value(var)
            at_bound = false
            
            # Check if at lower bound
            if has_lower_bound(var)
                at_bound = at_bound || abs(var_value - lower_bound(var)) < 1e-8
            end
            
            # Check if at upper bound  
            if has_upper_bound(var)
                at_bound = at_bound || abs(var_value - upper_bound(var)) < 1e-8
            end
            
            if at_bound
                push!(nonbasic_vars, var)
            else
                push!(basic_vars, var)
            end
        end
    end
    
    return basic_vars, nonbasic_vars
end

"""
    build_basic_matrix(model::Model, basic_vars::Vector{VariableRef}, 
                      constraint_refs::Vector{ConstraintRef})

Build the basic matrix A_B containing coefficients of basic variables 
in the structural constraints.
"""
function build_basic_matrix(model::Model, basic_vars::Vector{VariableRef}, 
                           constraint_refs::Vector{ConstraintRef})
    m = length(constraint_refs)
    n_basic = length(basic_vars)
    
    if n_basic != m
        error("Number of basic variables ($n_basic) must equal number of constraints ($m)")
    end
    
    A_B = zeros(Float64, m, m)
    
    for (i, constr) in enumerate(constraint_refs)
        constr_func = constraint_object(constr).func
        
        for (j, var) in enumerate(basic_vars)
            # Extract coefficient of var in constraint i
            coeff = 0.0
            
            if isa(constr_func, AffExpr)
                # Search through terms for this variable
                for (term_var, term_coeff) in constr_func.terms
                    if term_var == var
                        coeff = term_coeff
                        break
                    end
                end
            elseif isa(constr_func, VariableRef) && constr_func == var
                coeff = 1.0
            end
            
            A_B[i, j] = coeff
        end
    end
    
    return A_B
end

"""
    compute_basic_matrix_inverse(A_B::Matrix{Float64})

Compute the basic matrix inverse M = A_B^(-1) with numerical stabilization.
"""
function compute_basic_matrix_inverse(A_B::Matrix{Float64})
    m = size(A_B, 1)
    
    # Check condition number
    cond_num = cond(A_B)
    if cond_num > 1e12
        @warn "Basic matrix is poorly conditioned (cond = $cond_num)"
    end
    
    try
        # Compute inverse using LU factorization
        M = inv(A_B)
        return M
    catch e
        @warn "Failed to invert basic matrix: $e"
        @warn "Using regularized inverse"
        
        # Add small regularization term
        reg_matrix = A_B + 1e-12 * I
        return inv(reg_matrix)
    end
end

"""
    compute_stability_cone(basis::BasisStability)

Compute the stability cone K_B = {δD : M * δD >= -x*_B}.
Returns the H-representation as (A_ineq, b_ineq) where A_ineq * δD >= b_ineq.
"""
function compute_stability_cone(basis::BasisStability)
    # The stability cone is: M * δD >= -x_B
    # In standard form: A_ineq * δD >= b_ineq
    A_ineq = basis.M
    b_ineq = -basis.x_B
    
    return A_ineq, b_ineq
end

"""
    compute_independent_bounds(basis::BasisStability)

Compute independent bounds for each RHS component δD_j.
Returns intervals [δD_j_min, δD_j_max] for each component j.

For each j: δD_j ∈ [max_{i: M_ij > 0} (-x*_B_i / M_ij), min_{i: M_ij < 0} (-x*_B_i / M_ij)]
"""
function compute_independent_bounds(basis::BasisStability)
    M = basis.M
    x_B = basis.x_B
    m = basis.m
    
    bounds = Vector{Tuple{Float64, Float64}}()
    
    for j in 1:m  # For each RHS component
        M_j = M[:, j]  # j-th column of M
        
        lower_bound = -Inf
        upper_bound = Inf
        
        for i in 1:m  # For each row constraint
            if M_j[i] > 1e-12
                # M_ij > 0: constraint is δD_j >= -x_B_i / M_ij
                lower_bound = max(lower_bound, -x_B[i] / M_j[i])
            elseif M_j[i] < -1e-12
                # M_ij < 0: constraint is δD_j <= -x_B_i / M_ij  
                upper_bound = min(upper_bound, -x_B[i] / M_j[i])
            end
            # M_ij ≈ 0: no constraint on δD_j from this row
        end
        
        push!(bounds, (lower_bound, upper_bound))
    end
    
    return bounds
end

"""
    test_perturbation(basis::BasisStability, δD::Vector{Float64})

Test if perturbation vector δD lies within stability cone K_B.
Returns true if M * δD >= -x*_B (basis remains optimal).
"""
function test_perturbation(basis::BasisStability, δD::Vector{Float64})
    if length(δD) != basis.m
        error("Perturbation vector length ($(length(δD))) must equal number of constraints ($(basis.m))")
    end
    
    # Test: M * δD >= -x_B
    M_δD = basis.M * δD
    feasible = all(M_δD .>= -basis.x_B .- 1e-10)  # Small numerical tolerance
    
    return feasible
end

"""
    analyze_perturbation_sensitivity(basis::BasisStability; n_samples::Int=1000)

Perform Monte Carlo analysis of perturbation sensitivity within the stability cone.
"""
function analyze_perturbation_sensitivity(basis::BasisStability; n_samples::Int=1000)
    println("[ Info: Running perturbation sensitivity analysis with $n_samples samples...")
    
    bounds = compute_independent_bounds(basis)
    
    feasible_count = 0
    perturbation_norms = Float64[]
    
    for _ in 1:n_samples
        # Generate random perturbation within independent bounds
        δD = zeros(Float64, basis.m)
        
        for j in 1:basis.m
            if isfinite(bounds[j][1]) && isfinite(bounds[j][2])
                # Sample uniformly within bounds
                δD[j] = bounds[j][1] + rand() * (bounds[j][2] - bounds[j][1])
            elseif isfinite(bounds[j][1])
                # Only lower bound
                δD[j] = bounds[j][1] + abs(randn())
            elseif isfinite(bounds[j][2])
                # Only upper bound  
                δD[j] = bounds[j][2] - abs(randn())
            else
                # No bounds - use small random perturbation
                δD[j] = 0.1 * randn()
            end
        end
        
        if test_perturbation(basis, δD)
            feasible_count += 1
            push!(perturbation_norms, norm(δD))
        end
    end
    
    feasibility_rate = feasible_count / n_samples
    
    println("[ Results: Perturbation Analysis")
    println("   Feasibility rate: $(feasibility_rate*100)%")
    if !isempty(perturbation_norms)
        println("   Average feasible perturbation norm: $(mean(perturbation_norms))")
        println("   Max feasible perturbation norm: $(maximum(perturbation_norms))")
    end
    
    return (
        feasibility_rate = feasibility_rate,
        bounds = bounds,
        sample_norms = perturbation_norms
    )
end

"""
    print_stability_summary(basis::BasisStability)

Print a comprehensive summary of the basis stability analysis.
"""
function print_stability_summary(basis::BasisStability)
    println("\n=== BASIS STABILITY ANALYSIS ===")
    println("Model dimensions:")
    println("  Variables (n): $(basis.n)")
    println("  Constraints (m): $(basis.m)")
    println("  Basic variables: $(length(basis.basic_vars))")
    println("  Non-basic variables: $(length(basis.nonbasic_vars))")
    
    # Check matrix properties
    cond_M = cond(basis.M)
    println("\nMatrix properties:")
    println("  Basic matrix inverse condition number: ", @sprintf("%.2e", cond_M))
    println("  Basic solution norm: ", @sprintf("%.4f", norm(basis.x_B)))
    
    # Print sample of independent bounds
    bounds = compute_independent_bounds(basis)
    println("\nIndependent RHS perturbation bounds (first 10):")
    for j in 1:min(10, length(bounds))
        lb_str = isfinite(bounds[j][1]) ? @sprintf("%.4f", bounds[j][1]) : "-∞"
        ub_str = isfinite(bounds[j][2]) ? @sprintf("%.4f", bounds[j][2]) : "+∞"
        println("  δD[$j] ∈ [$lb_str, $ub_str]")
    end
    if length(bounds) > 10
        println("  ... and $(length(bounds)-10) more")
    end
    
    # Run sensitivity analysis
    sensitivity = analyze_perturbation_sensitivity(basis; n_samples=500)
    
    println("\n=== END STABILITY ANALYSIS ===\n")
end