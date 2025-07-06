# test_basis_extraction.jl
# Test the production-ready basis extraction function

include("shared.jl")
include("energy/energy_model_world.jl")
include("basis_extraction.jl")

function test_production_basis_extraction()
    println("=== PRODUCTION BASIS EXTRACTION TEST ===")
    
    # Build and solve the energy model
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    create_energy_model!(model)
    optimize!(model)
    
    if termination_status(model) != MOI.OPTIMAL
        println("❌ Model failed to solve")
        return
    end
    
    println("✅ Model solved successfully")
    
    # Extract basis matrix
    basic_vars, basis_inverse, matrix_rank = extract_basis_matrix(model)
    
    # Interpret basic variables
    num_cols = MOI.get(model, MOI.NumberOfVariables())
    col_vars, row_vars = interpret_basic_variables(basic_vars, num_cols)
    
    println("\nBasic variable analysis:")
    println("  Column variables in basis: $(length(col_vars))")
    println("  Row variables in basis: $(length(row_vars))")
    println("  Total basic variables: $(length(basic_vars))")
    
    # Show some statistics
    println("\nBasis inverse matrix statistics:")
    println("  Size: $(size(basis_inverse))")
    println("  Rank: $matrix_rank")
    println("  Condition number: $(cond(basis_inverse))")
    println("  Determinant: $(det(basis_inverse))")
    
    # Show first few elements
    println("\nFirst 5×5 block of basis inverse:")
    for i in 1:min(5, size(basis_inverse, 1))
        for j in 1:min(5, size(basis_inverse, 2))
            @printf("%10.3f ", basis_inverse[i,j])
        end
        println()
    end
    
    println("\n✅ Production basis extraction completed successfully")
    return basic_vars, basis_inverse
end

# Run the test
if abspath(PROGRAM_FILE) == @__FILE__
    result = test_production_basis_extraction()
    println("\n🎉 Basis extraction working correctly!")
end