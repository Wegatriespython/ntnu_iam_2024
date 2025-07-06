# basis_extraction.jl
# Production-ready basis matrix extraction from HiGHS

using JuMP, HiGHS
using Printf, LinearAlgebra
import MathOptInterface as MOI

"""
    extract_basis_matrix(model::Model) -> (basic_vars, basis_inverse)

Extract the basis matrix inverse from a solved HiGHS model.

Returns:
- basic_vars: Vector{Int32} of basic variable indices (HiGHS convention)
- basis_inverse: Matrix{Float64} of basis inverse (num_rows × num_rows)
"""
function extract_basis_matrix(model::Model)
    # Validate model is solved
    if termination_status(model) != MOI.OPTIMAL
        error("Model must be solved optimally before extracting basis")
    end
    
    # Get HiGHS instance
    backend_obj = JuMP.backend(model)
    bridge_optimizer = backend_obj.optimizer
    highs_optimizer = bridge_optimizer.model
    highs_ptr = highs_optimizer.inner
    
    if highs_ptr == C_NULL
        error("Invalid HiGHS instance pointer")
    end
    
    # Get problem dimensions
    num_rows = MOI.get(model, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}}()) +
               MOI.get(model, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}}()) +
               MOI.get(model, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}}())
    
    num_cols = MOI.get(model, MOI.NumberOfVariables())
    
    println("Extracting basis from model with $num_rows rows, $num_cols columns")
    
    # Get basic variables
    basic_vars = Vector{Int32}(undef, num_rows)
    status = HiGHS.Highs_getBasicVariables(highs_ptr, basic_vars)
    
    if status != 0
        error("Failed to get basic variables: status $status")
    end
    
    # Validate basic variables
    col_basics = count(x -> x >= 0, basic_vars)
    row_basics = count(x -> x < 0, basic_vars)
    println("Basic variables: $col_basics columns, $row_basics rows")
    
    # Extract basis inverse matrix
    println("Extracting basis inverse matrix...")
    
    # Pre-allocate arrays
    row_vector = Vector{Float64}(undef, num_rows)
    row_indices = Vector{Int32}(undef, num_rows)
    row_num_nz = Ref{Int32}(0)
    
    # Initialize basis inverse matrix
    basis_inverse = zeros(Float64, num_rows, num_rows)
    
    # Extract each row of basis inverse
    for row in 0:(num_rows-1)
        if row % 100 == 0
            println("  Processing row $row...")
        end
        
        # Clear arrays
        fill!(row_vector, 0.0)
        fill!(row_indices, -1)
        row_num_nz[] = 0
        
        # Get this row of basis inverse
        status = HiGHS.Highs_getBasisInverseRow(highs_ptr, Int32(row), row_vector, row_num_nz, row_indices)
        
        if status != 0
            error("Failed to get basis inverse row $row: status $status")
        end
        
        # Store non-zero elements
        for i in 1:row_num_nz[]
            col_idx = row_indices[i] + 1  # Convert to 1-based indexing
            if col_idx > 0 && col_idx <= num_rows
                # Only store actually non-zero values
                if abs(row_vector[i]) > 1e-15
                    basis_inverse[row + 1, col_idx] = row_vector[i]
                end
            end
        end
    end
    
    # Print statistics
    total_nonzeros = count(x -> abs(x) > 1e-12, basis_inverse)
    sparsity = (1.0 - total_nonzeros / (num_rows * num_rows)) * 100
    matrix_rank = rank(basis_inverse)
    println("Basis inverse: $(num_rows)×$(num_rows) matrix, $total_nonzeros non-zeros ($(round(sparsity, digits=1))% sparse)")
    println("Matrix rank: $matrix_rank")
    
    return basic_vars, basis_inverse, matrix_rank
end

"""
    interpret_basic_variables(basic_vars::Vector{Int32}, num_cols::Int) -> (col_vars, row_vars)

Interpret HiGHS basic variable indices according to HiGHS convention.

Returns:
- col_vars: Vector of column indices (0-based)
- row_vars: Vector of row indices (0-based, converted from negative values)
"""
function interpret_basic_variables(basic_vars::Vector{Int32}, num_cols::Int)
    col_vars = Int32[]
    row_vars = Int32[]
    
    for idx in basic_vars
        if idx >= 0
            push!(col_vars, idx)
        else
            push!(row_vars, -idx - 1)  # Convert negative to positive row index
        end
    end
    
    return col_vars, row_vars
end

"""
    verify_basis_inverse(basis_inverse::Matrix{Float64}, A::Matrix{Float64}, basic_vars::Vector{Int32})

Verify that the extracted basis inverse is correct by checking B^(-1) * B = I.
"""
function verify_basis_inverse(basis_inverse::Matrix{Float64}, A::Matrix{Float64}, basic_vars::Vector{Int32})
    num_rows = size(basis_inverse, 1)
    
    # Extract basis matrix columns
    basis_matrix = zeros(Float64, num_rows, num_rows)
    
    for (i, var_idx) in enumerate(basic_vars)
        if var_idx >= 0
            # Column variable
            if var_idx + 1 <= size(A, 2)
                basis_matrix[:, i] = A[:, var_idx + 1]
            end
        else
            # Row variable (slack/surplus) - identity column
            row_idx = -var_idx - 1
            if row_idx + 1 <= num_rows
                basis_matrix[row_idx + 1, i] = 1.0
            end
        end
    end
    
    # Check B^(-1) * B ≈ I
    product = basis_inverse * basis_matrix
    identity_check = norm(product - I, Inf)
    
    println("Basis verification: ||B^(-1) * B - I||_∞ = $(identity_check)")
    
    return identity_check < 1e-10
end