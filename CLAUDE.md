# NTNU Energy-Macro Integrated Assessment Model

## Project Overview

This project implements an integrated assessment model (IAM) for energy-economy interactions using Julia and JuMP. The model combines energy system optimization (MESSAGE framework) with macroeconomic modeling (MACRO framework) to analyze climate and energy policies.

## Model Components

### Core Modules
- **EnergyMacroModel.jl**: Main module containing all model components
- **Energy System Model**: Technology-rich energy system optimization with capacity constraints, emissions, and cost accounting
- **Macroeconomic Model**: CES production function with capital, labor, and energy services as inputs
- **Generalized Benders Decomposition (GBD)**: Decomposition algorithm splitting energy and macro optimization
- **Configuration System**: Centralized parameter management through YAML and Julia config files
- **Sensitivity Analysis**: Systematic parameter perturbation and analysis framework

### Key Files
- `julia_model/src/EnergyMacroModel.jl`: Module definition and exports
- `julia_model/shared.jl`: Common sets, parameters, and shared variable creation function
- `julia_model/energy/energy_model_world.jl`: Energy system model formulation (handles both standalone and integrated modes) 
- `julia_model/energy_sub_problem.jl`: Standalone energy subproblem implementation
- `julia_model/macro/macro_core.jl`: MACRO model implementation
- `julia_model/macro/macro_data_load.jl`: Economic parameters and data loading
- `julia_model/macro/macro_presolve.jl`: Macro model presolving routines
- `julia_model/model_config.jl`: Centralized configuration management
- `julia_model/energy_parameters.yaml`: Energy system parameter definitions
- `julia_model/run_energy.jl`: Standalone energy model runner
- `julia_model/run_MACRO.jl`: Standalone MACRO model runner
- `julia_model/run_energy_macro.jl`: Integrated energy-macro model runner
- `julia_model/run_energy_macro_gbd.jl`: GBD algorithm implementation

### Test Infrastructure
- Comprehensive test suite with 890+ tests covering parameter validation, component functionality, and integration
- Mathematical verification using SymPy for dimensional consistency and formulation validation
- Systematic debugging infrastructure for identifying infeasibility issues

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

### Time Parameter Definitions
**CRITICAL**: The model uses two distinct time parameters that must not be confused:

- **`duration_period`** (10 years): MACRO model time discretization
  - Used for: Macroeconomic dynamics, depreciation, growth calculations, utility discounting
  - Domain: Economic model (capital accumulation, labor growth, AEEI factors)
  - Files: `macro/macro_data_load.jl`, `macro/macro_core.jl`, economic parts of GBD

- **`period_length`** (10 years): MESSAGE model time discretization  
  - Used for: Technology lifetime constraints, capacity planning, energy system accounting
  - Domain: Energy model (technology constraints, emissions, energy costs)
  - Files: `energy/energy_model_world.jl`, `src/EnergyMacroModel.jl`, energy parts of GBD

**Note**: While both currently equal 10 years, they represent conceptually different time scales and should be used consistently within their respective domains to maintain model integrity.

### Mathematical Formulation
- **Energy Subproblem**: Minimizes total system cost subject to service demand constraints
- **Master Problem**: Maximizes utility subject to macroeconomic constraints and Benders cuts
- **Dual Value Extraction**: Shadow prices from energy balance constraints converted to marginal service costs
- **Convergence**: Relative gap tolerance between upper bound (energy cost) and lower bound (master objective)

## Current Status

### Working Components
- Parameter loading and validation (197/197 tests passing)
- Module structure with proper exports and imports
- Master problem construction and basic solving (203/227 tests passing)
- Benders cut generation and mathematical framework

### Identified Issues
- Energy subproblem infeasibility (24/30 tests failing)
- GBD algorithm convergence problems
- Numerical stability concerns in specific parameter ranges

### Recent Fixes
- Fixed duplicate COST_ANNUAL variable definition by creating shared variables function in shared.jl
- Energy-macro integrated model now runs successfully with optimal utility value of 179.50
- Results match GAMS model output for all macro variables (GDP, C, I) and utility

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

## Model Validation

### Mathematical Verification
- Dimensional consistency across billion/trillion USD scaling
- Production function calibration to base year data
- Energy-economy linkage through AEEI factors
- Dual variable extraction and Benders cut formulation

### Test Coverage
- Parameter validation: All economic and energy parameters
- Component testing: Energy subproblem, master problem, Benders cuts
- Integration testing: Full GBD algorithm execution
- Numerical stability: Solver tolerances, scaling factors, initialization
- Performance testing: Construction and solve times

## File Organization

```
julia_model/
├── src/EnergyMacroModel.jl          # Main module
├── shared.jl                        # Common definitions and shared variables
├── model_config.jl                  # Centralized configuration management
├── energy_parameters.yaml           # Energy system parameter definitions
├── energy/                          # Energy model components
│   └── energy_model_world.jl        # Energy system optimization model
├── energy_sub_problem.jl            # Standalone energy subproblem
├── macro/                           # Macroeconomic model components
│   ├── macro_core.jl                # MACRO model implementation
│   ├── macro_data_load.jl           # Economic parameter loading
│   └── macro_presolve.jl            # Macro presolving routines
├── run_energy.jl                    # Standalone energy model runner
├── run_MACRO.jl                     # Standalone MACRO model runner
├── run_energy_macro.jl              # Integrated energy-macro runner
├── run_energy_macro_gbd.jl          # GBD algorithm implementation
├── sensitivity/                     # Sensitivity analysis framework
│   ├── run_sensitivity.jl           # Sensitivity analysis runner
│   ├── analyze_sensitivity.jl       # Analysis and visualization
│   └── data/                        # Generated sensitivity data
├── compile/                         # Performance optimization
│   ├── create_sysimage.jl           # System image creation
│   └── precompile_workload.jl       # Precompilation workload
├── test/                            # Test suite
│   ├── runtests.jl                  # Test runner
│   ├── test_parameters.jl           # Parameter validation
│   ├── test_energy_subproblem.jl    # Energy model tests
│   ├── test_master_problem.jl       # MACRO model tests
│   ├── test_benders_cuts.jl         # GBD implementation tests
│   ├── test_numerical_stability.jl  # Stability analysis
│   └── test_integration.jl          # Full algorithm tests
├── gbd_docs/                        # Mathematical documentation
│   └── gbd_mathematical_verification.py # SymPy analysis
└── Project.toml                     # Package configuration
```

## New Features and Capabilities

### Sensitivity Analysis Framework
- **Systematic Parameter Perturbation**: Automated generation of perturbed parameter sets
- **Monte Carlo Analysis**: Multiple runs with statistical analysis of results
- **Technology Mix Analysis**: Detailed breakdown of energy technology contributions
- **Elasticity Analysis**: Comprehensive sensitivity metrics and visualizations

### Performance Optimization
- **Precompilation Support**: System image creation for faster startup times
- **Modular Architecture**: Clean separation of energy, macro, and sensitivity components
- **Configuration-Driven**: Easy parameter modification without code changes

### Enhanced Data Management
- **CSV Export**: Comprehensive results export for external analysis
- **YAML Configuration**: Human-readable parameter files with documentation
- **Structured Outputs**: Organized data files with timestamped results

## Command Line Memories

- `julia --project=@. filename.jl`: Run Julia script using the project's dependencies and environment
- `julia --project=@. sensitivity/run_sensitivity.jl`: Run sensitivity analysis with default parameters
- `julia --project=@. sensitivity/analyze_sensitivity.jl`: Analyze and visualize sensitivity results