# run_energy.jl - Standalone energy model using corrected GAMS-matching parameters

using JuMP, Ipopt

println("Starting standalone energy model (GAMS-equivalent)...")

# Load all definitions from energy_model_world.jl
include("shared.jl")
include("energy_model_world.jl")

# ------------------------------------------------------------------------------
# Create and solve model
# ------------------------------------------------------------------------------

model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 3)

# Variables (exact from GAMS lines 367-373)
@variable(model, ACT[technology, year_all] >= 0)  # technology activity
@variable(model, CAP_NEW[technology, year_all] >= 0)  # new technology capacity
@variable(model, EMISS[year_all])  # CO2 emissions
@variable(model, CUM_EMISS)  # cumulative CO2 emissions
@variable(model, TOTAL_COST)  # total discounted systems costs
@variable(model, COST_ANNUAL[year_all])  # costs per year

# Energy balance equations (exact from GAMS lines 396-398)
for e in energy, l in level, y in year_all
    if haskey(energy_level, (e, l)) && haskey(demand, (e, l))
        @constraint(model,
            sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology) -
            demand[(e, l)] * (gdp[y] / gdp[2020])^beta >= 0
        )
    elseif haskey(energy_level, (e, l))
        # No demand - net production must be zero
        @constraint(model,
            sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology) == 0
        )
    end
end

# Capacity balance (exact from GAMS lines 400-401)
for t in technology, y in year_all
    if haskey(hours, t) && haskey(lifetime, t)
        year_idx = findfirst(==(y), year_all)
        @constraint(model,
            ACT[t, y] <= sum(CAP_NEW[t, year_all[i]] * hours[t]
                           for i in 1:year_idx if (year_idx - i + 1) * period_length <= lifetime[t])
        )
    end
end

# Emissions (exact from GAMS lines 408-412)
@constraint(model, [y in year_all], 
    sum(ACT[t, y] * get(CO2_emission, t, 0) for t in technology) == EMISS[y])
@constraint(model, 
    sum(EMISS[y] * period_length for y in year_all) == CUM_EMISS)

# Diffusion constraints (exact from GAMS lines 414-415)
for t in technology, y in year_all[2:end]
    if haskey(diffusion_up, t)
        y_prev = year_all[findfirst(==(y), year_all) - 1]
        @constraint(model,
            CAP_NEW[t, y] <= CAP_NEW[t, y_prev] * (1 + diffusion_up[t])^period_length + get(startup, t, 0)
        )
    end
end

# Share constraints (exact from GAMS lines 417-421)
for s in share, y in year_all
    if haskey(share_up, s)
        lhs = [t for t in technology if haskey(tec_share, (s, t))]
        rhs = [t for t in technology if haskey(tec_share_rhs, (s, t))]
        @constraint(model, sum(ACT[t, y] for t in lhs) <= share_up[s] * sum(ACT[t, y] for t in rhs))
    end
    if haskey(share_lo, s)
        lhs = [t for t in technology if haskey(tec_share, (s, t))]
        rhs = [t for t in technology if haskey(tec_share_rhs, (s, t))]
        @constraint(model, sum(ACT[t, y] for t in lhs) >= share_lo[s] * sum(ACT[t, y] for t in rhs))
    end
end

# Cost equations (exact from GAMS lines 423-436)
year_idx = Dict(y => i for (i, y) in enumerate(year_all))
for y in year_all
    @constraint(model,
        sum(ACT[t, y] * get(vom, (t, y), 0) for t in technology) +
        sum(sum(CAP_NEW[t, y2] * get(cost_capacity, (t, y2), 0) 
               for y2 in year_all 
               if year_idx[y2] <= year_idx[y] && (year_idx[y] - year_idx[y2] + 1) * period_length <= get(lifetime, t, 0))
            for t in technology if haskey(lifetime, t))
        == COST_ANNUAL[y]
    )
end

@constraint(model,
    sum(COST_ANNUAL[y] * period_length * (1 - discount_rate)^(period_length * (year_idx[y] - 1)) for y in year_all)
    == TOTAL_COST
)

# Historical calibration (exact from GAMS lines 461-474)
fix(ACT["coal_ppl", 2020], 9.462; force=true)
fix(ACT["oil_ppl", 2020], 0.7; force=true)
fix(ACT["solar_PV_ppl", 2020], 0.839; force=true)
fix(ACT["gas_ppl", 2020], 6.36; force=true)
fix(ACT["nuclear_ppl", 2020], 2.68; force=true)
fix(ACT["hydro_ppl", 2020], 4.36; force=true)
fix(ACT["wind_ppl", 2020], 1.6; force=true)
fix(ACT["bio_ppl", 2020], 0.69; force=true)
set_lower_bound(ACT["other_ppl", 2020], 0.127)
fix(ACT["coal_nele", 2020], 10.7; force=true)
fix(ACT["oil_nele", 2020], 43.0; force=true)
fix(ACT["gas_nele", 2020], 18.7; force=true)
fix(ACT["bio_nele", 2020], 10.6; force=true)
set_lower_bound(ACT["other_nele", 2020], 0.28)

# Objective (minimize TOTAL_COST)
@objective(model, Min, TOTAL_COST)

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
    total_cost = value(TOTAL_COST)
    println("✅ Optimal solution found")
    println("Total system cost: $(round(total_cost, digits=2)) billion USD")
    
    println("\nKey results for 2020:")
    for t in ["coal_ppl", "gas_ppl", "oil_ppl", "bio_ppl", "hydro_ppl", "wind_ppl", "solar_PV_ppl", "nuclear_ppl"]
        if t in technology
            act = value(ACT[t, 2020])
            if act > 0.01
                println("  ACT[$t] = $(round(act, digits=2)) PWh")
            end
        end
    end
    
    println("\nNon-electric technologies 2020:")
    for t in ["coal_nele", "gas_nele", "oil_nele", "bio_nele", "solar_nele", "other_nele"]
        if t in technology
            act = value(ACT[t, 2020])
            if act > 0.01
                println("  ACT[$t] = $(round(act, digits=2)) PWh")
            end
        end
    end
    
    println("\nAnnual costs (billion USD):")
    for y in year_all
        annual_cost = value(COST_ANNUAL[y])
        println("  $y: $(round(annual_cost, digits=1)) billion USD")
    end
    
    println("\nEmissions (Mt CO2):")
    for y in year_all
        emiss = value(EMISS[y])
        println("  $y: $(round(emiss, digits=1)) Mt CO2")
    end
    
    # Compare with GAMS result
    gams_result = 97981.69
    difference = total_cost - gams_result
    println("\nComparison with GAMS:")
    println("  GAMS result: $(round(gams_result, digits=2)) billion USD")
    println("  Julia result: $(round(total_cost, digits=2)) billion USD")
    println("  Difference: $(round(difference, digits=2)) billion USD ($(round(difference/gams_result*100, digits=2))%)")
    
else
    println("❌ Model failed: $(termination_status(model))")
end

println("\n✓ Energy model completed")