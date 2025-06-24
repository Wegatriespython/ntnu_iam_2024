# energy_sub_problem.jl - Energy system optimization subproblem
# Can be solved standalone as LP least cost optimization or used in GBD decomposition

using JuMP, Ipopt, LinearAlgebra

# Create energy subproblem - flexible for standalone or GBD use
function create_energy_subproblem(S_bar::Union{Dict{Tuple{String,Int},Float64}, Nothing} = nothing; 
                                 config = nothing,
                                 standalone::Bool = false)
    
    # Use default config if not provided - fallback to basic solver settings
    if config === nothing
        solver_tolerance = 1e-5
        solver_max_iter = 2000
    else
        solver_tolerance = config.solver_tolerance
        solver_max_iter = config.solver_max_iter
    end
    
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", standalone ? 5 : 0)
    set_optimizer_attribute(model, "tol", solver_tolerance)
    set_optimizer_attribute(model, "max_iter", solver_max_iter)

    # Decision variables
    @variable(model, ACT[technology, year_all] >= 0)
    @variable(model, CAP_NEW[technology, year_all] >= 0)
    @variable(model, EMISS[year_all])
    @variable(model, CUM_EMISS)
    @variable(model, TOTAL_COST)                # **billion US$**
    @variable(model, COST_ANNUAL[year_all])
    @variable(model, demand_slack[sector, year_all] >= 0)

    # Bookkeeping for sector-split mapping
    n_map = Dict(s => sum(haskey(map_energy_sector, (e, l, s)) for e in energy, l in level)
                 for s in sector)

    # Energy balance constraints
    bal = Dict{Tuple{String,String,Int},ConstraintRef}()
    for e in energy, l in level, y in year_all
        # Check if this energy-level combination is valid
        if !haskey(energy_level, (e, l))
            continue
        end
        
        if standalone
            # Standalone mode: use standard demands or zero for intermediate commodities
            if haskey(demand, (e, l))
                # This has a standard demand (non-PHYSENE commodities)
                bal[(e, l, y)] = @constraint(model,
                    sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology)
                    >=
                    demand[(e, l)] * (gdp[y] / gdp[2020])^beta)
            elseif haskey(map_energy_sector, (e, l, "ELEC")) || haskey(map_energy_sector, (e, l, "NELE"))
                # For PHYSENE commodities in standalone mode, use demand_MESSAGE
                sector_mapped = haskey(map_energy_sector, (e, l, "ELEC")) ? "ELEC" : "NELE"
                if haskey(demand_MESSAGE, (sector_mapped, y))
                    bal[(e, l, y)] = @constraint(model,
                        sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology)
                        >=
                        demand_MESSAGE[(sector_mapped, y)])
                else
                    # No demand - net production must be zero
                    bal[(e, l, y)] = @constraint(model,
                        sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology)
                        ==
                        0)
                end
            else
                # No demand - net production must be zero (no free disposal)
                bal[(e, l, y)] = @constraint(model,
                    sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology)
                    ==
                    0)
            end
        else
            # GBD mode: use fixed service demands S_bar
            if S_bar === nothing
                error("S_bar must be provided for GBD mode (standalone=false)")
            end
            
            # Check if this maps to a sector (PHYSENE demand)
            if haskey(map_energy_sector, (e, l, "ELEC")) && haskey(S_bar, ("ELEC", y))
                # This is electricity useful energy - use ELEC demand
                bal[(e, l, y)] = @constraint(model,
                    sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology)
                    >=
                    S_bar[("ELEC", y)])
            elseif haskey(map_energy_sector, (e, l, "NELE")) && haskey(S_bar, ("NELE", y))
                # This is non-electric useful energy - use NELE demand  
                bal[(e, l, y)] = @constraint(model,
                    sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology)
                    >=
                    S_bar[("NELE", y)])
            elseif haskey(demand, (e, l))
                # This has a standard demand (non-PHYSENE commodities)
                bal[(e, l, y)] = @constraint(model,
                    sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology)
                    >=
                    demand[(e, l)] * (gdp[y] / gdp[2020])^beta)
            else
                # No demand - net production must be zero (no free disposal)
                bal[(e, l, y)] = @constraint(model,
                    sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology)
                    ==
                    0)
            end
        end
    end

    # Capacity & technical constraints
    for t in technology, y in year_all
        if haskey(hours, t) && haskey(lifetime, t)
            idx = findfirst(==(y), year_all)
            @constraint(model,
                ACT[t, y] <= sum(CAP_NEW[t, year_all[i]] * hours[t]
                               for i in 1:idx if (idx - i + 1) * period_length <= lifetime[t]))
        end
    end
    for t in technology, y in year_all[2:end]
        if !haskey(diffusion_up, t)
            continue
        end
        y_prev = year_all[findfirst(==(y), year_all)-1]
        @constraint(model, CAP_NEW[t, y] <= CAP_NEW[t, y_prev] * (1 + diffusion_up[t])^period_length + get(startup, t, 0))
    end

    # Share constraints
    for s in share, y in year_all
        lhs = [t for t in technology if haskey(tec_share, (s, t))]
        rhs = [t for t in technology if haskey(tec_share_rhs, (s, t))]
        if haskey(share_up, s)
            @constraint(model, sum(ACT[t, y] for t in lhs) <= share_up[s] * sum(ACT[t, y] for t in rhs))
        end
        if haskey(share_lo, s)
            @constraint(model, sum(ACT[t, y] for t in lhs) >= share_lo[s] * sum(ACT[t, y] for t in rhs))
        end
    end

    # Emissions
    @constraint(model, [y in year_all], sum(ACT[t, y] * get(CO2_emission, t, 0) for t in technology) == EMISS[y])
    @constraint(model, sum(EMISS[y] * period_length for y in year_all) == CUM_EMISS)

    # Cost accounting
    year_idx = Dict(y => i for (i, y) in enumerate(year_all))
    lifetime_r = Dict{Tuple{String,Int},Vector{Int}}()
    for t in technology, y in year_all
        if !haskey(lifetime, t)
            continue
        end
        idx, lt = year_idx[y], lifetime[t]
        lifetime_r[(t, y)] = [y2 for y2 in year_all if year_idx[y2] <= idx && (idx - year_idx[y2] + 1) * period_length <= lt]
    end
    disc = Dict(y => (1 - drate)^(period_length * (year_idx[y] - 1)) for y in year_all)
    @constraint(model, [y in year_all],
        sum(ACT[t, y] * get(vom, (t, y), 0) for t in technology) +
        sum(sum(CAP_NEW[t, y2] * get(cost_capacity, (t, y2), 0) for y2 in lifetime_r[(t, y)]) for t in keys(lifetime))
        ==
        COST_ANNUAL[y])
    @constraint(model, sum(COST_ANNUAL[y] * period_length * disc[y] for y in year_all) == TOTAL_COST)

    # Historical calibration
    for ((t, y), v) in energy_calibration
        fix(ACT[t, y], v; force=true)
    end
    for ((t, y), v) in energy_calibration_lo
        set_lower_bound(ACT[t, y], v)
    end

    # Objective: minimize cost with slack penalty
    slack_penalty = standalone ? 1e3 : 1e6  # Lower penalty for standalone to aid convergence
    @objective(model, Min, TOTAL_COST + slack_penalty * sum(demand_slack))
    
    return model, bal, demand_slack
end

# Standalone solver function
function solve_energy_standalone(;config = nothing, verbose::Bool = true)
    println("Solving standalone energy system optimization...")
    
    # Load required data if not already loaded
    if !@isdefined(year_all)
        include("model_config.jl")
        include("shared.jl")
        include("energy_model_world.jl")
        include("macro_data_load.jl")
    end
    
    model, bal, demand_slack = create_energy_subproblem(nothing; config=config, standalone=true)
    
    if verbose
        println("Model has $(num_variables(model)) variables")
        println("Model has $(num_constraints(model; count_variable_in_set_constraints=false)) constraints")
    end
    
    optimize!(model)
    status = termination_status(model)
    
    if verbose
        println("Termination status: $status")
    end
    
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        println("❌ Energy subproblem failed: $status")
        return nothing
    end
    
    total_cost = value(model[:TOTAL_COST])
    
    if verbose
        println("✅ Energy system solved successfully")
        println("Total system cost: $(round(total_cost, digits=2)) billion USD")
        
        # Check demand slack usage
        total_slack = sum(value(demand_slack[s, y]) for s in sector, y in year_all)
        if total_slack > 0.001
            println("⚠️  Demand slack used: $(round(total_slack, digits=3)) PWh")
        end
        
        # Show key results
        println("\nKey results for 2020:")
        for t in ["coal_ppl", "gas_ppl", "electricity_grid", "appliances"]
            if t in technology
                act = value(model[:ACT][t, 2020])
                if act > 0.01
                    println("  ACT[$t] = $(round(act, digits=2)) PWh")
                end
            end
        end
    end
    
    return model
end

# Driver for standalone execution
if abspath(PROGRAM_FILE) == @__FILE__
    solve_energy_standalone()
end