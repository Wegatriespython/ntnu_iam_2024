# NTNU Energy-Macro Integrated Assessment Model

## Project Overview

This project implements an integrated assessment model (IAM) for energy-economy interactions using Julia and JuMP. The model combines energy system optimization (MESSAGE framework) with macroeconomic modeling (MACRO framework) to analyze climate and energy policies.

## Model Components

### Core Modules

- **EnergyMacroModel.jl**: Main module containing all model components
- **Energy System Model**: Technology-rich energy system optimization with capacity constraints, emissions, and cost accounting
- **Macroeconomic Model**: CES production function with capital, labor, and energy services as inputs
- **Dynamic Programming MACRO**: Discrete-time dynamic programming formulation with vintage capital structure
- **Configuration System**: Centralized parameter management through YAML and Julia config files

### Key Files

- `julia_model/src/EnergyMacroModel.jl`: Module definition and exports
- `julia_model/shared.jl`: Common sets, parameters, and shared variable creation function
- `julia_model/energy/energy_model_world.jl`: Energy system model formulation (handles both standalone and integrated modes)
- `julia_model/macro/macro_core.jl`: MACRO model implementation
- `julia_model/macro/macro_data_load.jl`: Economic parameters and data loading
- `julia_model/macro/macro_presolve.jl`: Macro model presolving routines
- `julia_model/model_config.jl`: Centralized configuration management
- `julia_model/energy_parameters.yaml`: Energy system parameter definitions
- `julia_model/run_energy.jl`: Standalone energy model runner
- `julia_model/run_MACRO.jl`: Standalone MACRO model runner
- `julia_model/run_energy_macro.jl`: Integrated energy-macro model runner
- `julia_model/dp/dp_macro_model.jl`: Dynamic programming MACRO model with vintage capital structure
- `julia_model/dp/energy_cost_surrogate.jl`: Energy cost surrogate models (polynomial and AR processes)
- `julia_model/dp/dp_utilities.jl`: Utility functions for DP model operations
- `julia_model/run_dp_macro.jl`: Dynamic programming MACRO model runner

## Technical Implementation

### Model Structure

- **Time Horizon**: 2020-2080 in 10-year periods
- **Energy Sectors**: Electricity (ELEC) and Non-electricity (NELE)
- **Technologies**: 38 energy technologies including renewables, fossil fuels, and nuclear
- **Economic Framework**: Nested CES production function with autonomous energy efficiency improvement (AEEI)

### Scaling Conventions

- Energy system costs: Billion USD
- Macroeconomic variables: Trillion USD
- Energy services: Petawatt-hours (PWh)
- Automatic unit conversion in GBD interface (1000x scaling factor)

### Configuration Management

- **Centralized Configuration**: All model parameters managed through `model_config.jl` with `ModelConfig` struct
- **YAML Parameter Files**: Energy system parameters defined in `energy_parameters.yaml` for easy modification
- **Flexible Overrides**: Custom configurations support user-defined parameter modifications
- **Derived Parameters**: Automatic calculation of growth factors, AEEI factors, and economic coefficients

### Mathematical Formulation

- Read dp_mathematical_formulation.tex

#### Dynamic Programming Implementation

- **State Variables**: Capital stock (K), previous GDP (Y), production energy (PRODENE)
- **Control Variables**: Investment (I), new production (YN)
- **Bellman Equation**: Backward induction with value function interpolation
- **Envelope Conditions**: Analytical derivatives for first-order conditions
- **Energy Cost Function**: Surrogate models approximating energy system optimization
- **Grid Resolution**: 30×20×15 state space grid with 9,000 total states
- **Energy Endogeneity**: Fixed-point iteration solves for optimal energy allocation given production targets

## Current Status

### Development Infrastructure

- Git version control with systematic staging and commits
- Julia package structure with Project.toml dependencies
- Modular design preventing code duplication through include guards
- Comprehensive error handling and debugging capabilities

## Dependencies

### Required Packages

- JuMP.jl: Mathematical optimization modeling
- Ipopt.jl: Interior-point nonlinear solver
- LinearAlgebra.jl: Matrix operations
- CSV.jl, DataFrames.jl, DataFramesMeta.jl: Data handling and manipulation
- JSON.jl: Configuration and results export
- YAML.jl: Configuration file parsing
- GAMS.jl: Interface to GAMS optimization software
- ArgParse.jl: Command-line argument parsing
- Statistics.jl: Statistical functions
- Test.jl: Test infrastructure

### Development Tools

- SymPy (Python): Mathematical verification and analysis
- PackageCompiler.jl: Performance optimization through precompilation

## Data Sources

### Economic Parameters

- Base year GDP: 71 trillion USD (2020)
- Capital-GDP ratio: 2.8
- Depreciation rate: 5% annually
- Social discount rate: 5% annually
- Elasticity of substitution: 0.3

### Energy System Data

- Technology capacities, costs, and lifetimes
- Resource potentials and extraction costs
- Emissions factors and climate constraints
- Historical calibration data for base year

## File Organization

```
julia_model/
├── src/EnergyMacroModel.jl          # Main module
├── shared.jl                        # Common definitions and shared variables
├── model_config.jl                  # Centralized configuration management
├── energy_parameters.yaml           # Energy system parameter definitions
├── energy/                          # Energy model components
│   └── energy_model_world.jl        # Energy system optimization model
├── macro/                           # Macroeconomic model components
│   ├── macro_core.jl                # MACRO model implementation
│   ├── macro_data_load.jl           # Economic parameter loading
│   └── macro_presolve.jl            # Macro presolving routines
├── run_energy.jl                    # Standalone energy model runner
├── run_MACRO.jl                     # Standalone MACRO model runner
├── run_energy_macro.jl              # Integrated energy-macro runner
├── compile/                         # Performance optimization
│   ├── create_sysimage.jl           # System image creation
│   └── precompile_workload.jl       # Precompilation workload
├── dp/                              # Dynamic programming implementation
│   ├── dp_macro_model.jl            # DP MACRO model with vintage capital
│   ├── energy_cost_surrogate.jl     # Energy cost surrogate models
│   └── dp_utilities.jl              # DP utility functions
├── run_dp_macro.jl                  # DP MACRO model runner
└── Project.toml                     # Package configuration
```

## Implementation Details

### Dynamic Programming Model Structure

- **DPMacroModel**: Struct containing state grids, parameters, value functions, and envelope derivatives
- **BellmanProblem**: Container for model and energy cost function
- **State Space**: Three-dimensional grid over capital, GDP, and production energy
- **Value Function Storage**: Four-dimensional arrays V[K,Y,P,t] with envelope derivatives V_K, V_Y, V_P
- **Policy Functions**: Storage for optimal controls (consumption, investment, new production)

### Energy Cost Surrogate Models

- **Simple Linear**: Taylor expansion around MESSAGE base case using fixed shadow prices
- **Polynomial**: Higher-order polynomial fitting to energy system responses
- **AR Process**: Autoregressive time series models (AR1, AR2) with lagged cost terms
- **Integration**: Surrogate models called within DP Bellman operator for energy cost evaluation

### DP vs NLP Formulation Differences

- **NLP Model**: Endogenous energy system optimization within JuMP/GAMS framework
- **DP Model**: Exogenous energy cost calculation via surrogate functions
- **Energy Allocation**: DP uses fixed-point iteration to solve for optimal ELEC/NELE ratios using marginal prices from Taylor expansion
- **Price Dynamics**: DP extracts marginal costs ∂EC/∂ELEC and ∂EC/∂NELE from quadratic cost terms
- **Convergence**: DP achieves lower growth rates and higher energy cost shares compared to NLP

### Numerical Methods

- **Backward Induction**: Dynamic programming solution starting from terminal period
- **First-Order Conditions**: Analytical FOCs using envelope theorem for optimization
- **Interpolation**: Linear and cubic spline interpolation for value function approximation
- **Finite Differences**: Gradient approximation for energy cost derivatives
- **Newton Method**: FOC-based solver with grid search fallback

## Results Comparison

### NLP Results (run_MACRO.jl)
- **GDP Growth**: 71 → 373 trillion USD (2020-2080)
- **Utility**: 178.72
- **Energy Cost Share**: 7.1% → 1.3% of GDP
- **Growth Pattern**: Sustained growth throughout horizon

### DP Results (run_dp_macro.jl)
- **GDP Growth**: 71 → 135 trillion USD (2020-2080) 
- **Utility**: 163.18
- **Energy Cost Share**: 7.1% → 3.4% of GDP
- **Growth Pattern**: Stagnation after 2040, declining capital stock by 2080

## Command Line Memories

- `julia --project=@. filename.jl`: Run Julia script using project dependencies and environment
- `julia --project=@. run_dp_macro.jl`: Run dynamic programming MACRO model
- `julia --project=@. run_MACRO.jl`: Run NLP MACRO model

