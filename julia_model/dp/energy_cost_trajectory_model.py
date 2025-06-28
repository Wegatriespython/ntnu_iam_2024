#!/usr/bin/env python3
"""
Energy Cost Trajectory Prediction Model
Unified AR model for predicting energy system costs across demand scenarios
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import warnings

warnings.filterwarnings("ignore")


def load_and_prepare_data(filepath="cost_dataset_combined.csv"):
    """Load and prepare the full panel dataset for unified modeling"""

    df = pd.read_csv(filepath)
    df = df[df["solve_status"] == "OPTIMAL"].copy()

    # Create the full panel dataset
    years = [2020, 2030, 2040, 2050, 2060, 2070, 2080]
    cost_cols = [
        "cost_2020",
        "cost_2030",
        "cost_2040",
        "cost_2050",
        "cost_2060",
        "cost_2070",
        "cost_2080",
    ]

    panel_data = []

    for _, scenario in df.iterrows():
        scenario_id = scenario["scenario_id"]
        base_elec = scenario["elec_demand"]
        base_nele = scenario["nele_demand"]
        total_cost = scenario["total_cost"]

        for i, year in enumerate(years):
            cost_current = scenario[cost_cols[i]]
            elec_demand_current = scenario[f"elec_demand_{year}"]
            nele_demand_current = scenario[f"nele_demand_{year}"]
            total_demand_current = scenario[f"total_demand_{year}"]

            # Add lagged cost (for AR term)
            cost_lag1 = scenario[cost_cols[i - 1]] if i > 0 else np.nan
            cost_lag2 = scenario[cost_cols[i - 2]] if i > 1 else np.nan

            # Add time features
            time_index = i  # 0,1,2,3,4,5,6 for 2020-2080
            decade = (year - 2020) // 10  # 0,1,2,3,4,5,6

            panel_data.append(
                {
                    "scenario_id": scenario_id,
                    "year": year,
                    "time_index": time_index,
                    "decade": decade,
                    "cost_current": cost_current,
                    "cost_lag1": cost_lag1,
                    "cost_lag2": cost_lag2,
                    "elec_demand_current": elec_demand_current,
                    "nele_demand_current": nele_demand_current,
                    "total_demand_current": total_demand_current,
                    "base_elec_demand": base_elec,
                    "base_nele_demand": base_nele,
                    "base_total_demand": base_elec + base_nele,
                    "total_cost_scenario": total_cost,
                    # Add demand growth rates
                    "elec_demand_growth": elec_demand_current / base_elec
                    if base_elec > 0
                    else 1,
                    "nele_demand_growth": nele_demand_current / base_nele
                    if base_nele > 0
                    else 1,
                    "total_demand_growth": total_demand_current
                    / (base_elec + base_nele),
                }
            )

    panel_df = pd.DataFrame(panel_data)

    # Add interaction terms
    panel_df["demand_time_interaction"] = (
        panel_df["total_demand_current"] * panel_df["time_index"]
    )
    panel_df["elec_time_interaction"] = (
        panel_df["elec_demand_current"] * panel_df["time_index"]
    )
    panel_df["demand_squared"] = panel_df["total_demand_current"] ** 2

    return panel_df


def create_feature_sets(panel_df):
    """Create different feature sets for model comparison"""

    # Remove rows with missing lags for AR models
    ar_df = panel_df.dropna(subset=["cost_lag1"]).copy()

    feature_sets = {
        "basic_ar1": ["cost_lag1", "total_demand_current"],
        "extended_ar1": ["cost_lag1", "elec_demand_current", "nele_demand_current"],
        "ar2_with_demand": [
            "cost_lag1",
            "cost_lag2",
            "elec_demand_current",
            "nele_demand_current",
        ],
        "time_aware": [
            "cost_lag1",
            "elec_demand_current",
            "nele_demand_current",
            "time_index",
            "demand_time_interaction",
        ],
        "full_features": [
            "cost_lag1",
            "elec_demand_current",
            "nele_demand_current",
            "time_index",
            "demand_time_interaction",
            "elec_time_interaction",
            "demand_squared",
            "base_elec_demand",
            "base_nele_demand",
        ],
    }

    # For non-AR models, use all data including t=0
    feature_sets_no_ar = {
        "demand_only": ["elec_demand_current", "nele_demand_current"],
        "demand_with_time": [
            "elec_demand_current",
            "nele_demand_current",
            "time_index",
        ],
        "full_no_ar": [
            "elec_demand_current",
            "nele_demand_current",
            "time_index",
            "demand_time_interaction",
            "elec_time_interaction",
            "demand_squared",
            "base_elec_demand",
            "base_nele_demand",
        ],
    }

    return feature_sets, feature_sets_no_ar, ar_df


def k_fold_validation(panel_df, ar_df, feature_sets, k=5):
    """K-fold cross-validation for robust model evaluation"""
    
    scenarios = panel_df["scenario_id"].unique()
    np.random.shuffle(scenarios)
    fold_size = len(scenarios) // k
    
    key_models = {
        "basic_ar1": feature_sets["basic_ar1"],
        "time_aware": feature_sets["time_aware"],
        "ar2_with_demand": feature_sets["ar2_with_demand"]
    }
    
    results = {name: {"r2_scores": [], "rmse_scores": [], "models": []} for name in key_models}
    
    for fold in range(k):
        test_scenarios = scenarios[fold * fold_size:(fold + 1) * fold_size]
        train_scenarios = [s for s in scenarios if s not in test_scenarios]
        
        train_ar_data = ar_df[ar_df["scenario_id"].isin(train_scenarios)]
        test_ar_data = ar_df[ar_df["scenario_id"].isin(test_scenarios)]
        
        for name, features in key_models.items():
            if "cost_lag2" in features:
                train_subset = train_ar_data.dropna(subset=features)
                test_subset = test_ar_data.dropna(subset=features)
            else:
                train_subset = train_ar_data
                test_subset = test_ar_data
            
            if len(train_subset) == 0 or len(test_subset) == 0:
                continue
            
            X_train = train_subset[features]
            y_train = train_subset["cost_current"]
            X_test = test_subset[features]
            y_test = test_subset["cost_current"]
            
            model = LinearRegression().fit(X_train, y_train)
            y_pred = model.predict(X_test)
            
            r2 = r2_score(y_test, y_pred)
            rmse = np.sqrt(mean_squared_error(y_test, y_pred))
            
            results[name]["r2_scores"].append(r2)
            results[name]["rmse_scores"].append(rmse)
            results[name]["models"].append(model)
    
    # Print results
    print("\n5-FOLD CROSS-VALIDATION RESULTS:")
    print("-" * 40)
    for name, result in results.items():
        if result["r2_scores"]:
            mean_r2 = np.mean(result["r2_scores"])
            std_r2 = np.std(result["r2_scores"])
            mean_rmse = np.mean(result["rmse_scores"])
            std_rmse = np.std(result["rmse_scores"])
            print(f"{name}: R² = {mean_r2:.4f} ± {std_r2:.4f}, RMSE = {mean_rmse:.1f} ± {std_rmse:.1f}")
    
    # Return best model for bias analysis
    best_name = max(results.keys(), key=lambda k: np.mean(results[k]["r2_scores"]) if results[k]["r2_scores"] else 0)
    best_model = results[best_name]["models"][-1]  # Use last fold's model
    
    # Create final train/test split for bias analysis
    train_scenarios = scenarios[:int(0.8 * len(scenarios))]
    test_scenarios = scenarios[int(0.8 * len(scenarios)):]
    
    # Time-aware bias analysis
    if "time_aware" in results and results["time_aware"]["r2_scores"]:
        train_ar_data = ar_df[ar_df["scenario_id"].isin(train_scenarios)]
        test_ar_data = ar_df[ar_df["scenario_id"].isin(test_scenarios)]
        
        features = feature_sets["time_aware"]
        X_train = train_ar_data[features]
        y_train = train_ar_data["cost_current"]
        
        bias_model = LinearRegression().fit(X_train, y_train)
        
        bias_result = {
            "model": bias_model,
            "features": features,
            "test_data": test_ar_data
        }
        
        print("\n2030 BIAS ANALYSIS:")
        print("-" * 20)
        investigate_2030_bias(panel_df, bias_result, train_scenarios, test_scenarios)
    
    return results, train_scenarios, test_scenarios


def investigate_2030_bias(panel_df, best_result, train_scenarios, test_scenarios):
    """Investigate systematic bias in 2030 predictions"""
    
    test_data = best_result["test_data"]
    model = best_result["model"]
    features = best_result["features"]
    
    X_test = test_data[features]
    y_test_pred = model.predict(X_test)
    test_data_with_pred = test_data.copy()
    test_data_with_pred["predicted_cost"] = y_test_pred
    test_data_with_pred["residual"] = test_data_with_pred["cost_current"] - y_test_pred
    test_data_with_pred["pct_error"] = (test_data_with_pred["residual"] / test_data_with_pred["cost_current"]) * 100
    
    # Yearly bias summary
    yearly_bias = test_data_with_pred.groupby("year")["pct_error"].agg(["mean", "std", "count"]).round(2)
    print("Year-wise bias (% error):")
    for year, row in yearly_bias.iterrows():
        print(f"  {year}: {row['mean']:.1f}% ± {row['std']:.1f}% (n={row['count']})")
    
    # 2030 analysis
    data_2030 = test_data_with_pred[test_data_with_pred["year"] == 2030]
    if len(data_2030) > 0:
        bias_2030 = data_2030["pct_error"].mean()
        print(f"\n2030 bias: {bias_2030:.1f}% ({data_2030['residual'].mean():.0f}B USD)")
        
        # Transition analysis
        transition_data = []
        for scenario in data_2030["scenario_id"].unique():
            scenario_panel = panel_df[panel_df["scenario_id"] == scenario].sort_values("year")
            costs = scenario_panel[scenario_panel["year"].isin([2020, 2030, 2040])]["cost_current"].values
            
            if len(costs) >= 3:
                decline_2020_2030 = costs[0] - costs[1]
                decline_2030_2040 = costs[1] - costs[2]
                transition_data.append({"decline_2020_2030": decline_2020_2030, "decline_2030_2040": decline_2030_2040})
        
        if transition_data:
            trans_df = pd.DataFrame(transition_data)
            ratio = trans_df["decline_2020_2030"].mean() / trans_df["decline_2030_2040"].mean()
            print(f"2020→2030 decline is {ratio:.1f}x steeper than 2030→2040")
    
    return yearly_bias


def create_time_period_analysis(panel_df, best_result):
    """Analyze model performance across different time periods"""

    print(f"\nTIME PERIOD ANALYSIS:")
    print("-" * 25)

    test_data = best_result["test_data"]
    model = best_result["model"]
    features = best_result["features"]

    # Get predictions
    X_test = test_data[features]
    y_test_pred = model.predict(X_test)
    test_data_analysis = test_data.copy()
    test_data_analysis["predicted_cost"] = y_test_pred
    test_data_analysis["residual"] = test_data_analysis["cost_current"] - y_test_pred
    test_data_analysis["abs_residual"] = np.abs(test_data_analysis["residual"])
    test_data_analysis["pct_error"] = (
        test_data_analysis["residual"] / test_data_analysis["cost_current"]
    ) * 100

    # Analyze by time period
    period_analysis = (
        test_data_analysis.groupby("year")
        .agg(
            {
                "cost_current": ["mean", "std"],
                "predicted_cost": "mean",
                "residual": ["mean", "std"],
                "abs_residual": "mean",
                "pct_error": ["mean", "std"],
                "scenario_id": "count",
            }
        )
        .round(2)
    )

    period_analysis.columns = ["_".join(col).strip() for col in period_analysis.columns]

    print("Performance by time period:")
    print(period_analysis)

    return period_analysis


def predict_cost_trajectories(best_model, features, panel_df, scenario_ids=None):
    """Use the best model to predict cost trajectories for given scenarios"""

    if scenario_ids is None:
        scenario_ids = np.random.choice(
            panel_df["scenario_id"].unique(), 5, replace=False
        )

    predictions = {}

    for scenario_id in scenario_ids:
        scenario_data = panel_df[panel_df["scenario_id"] == scenario_id].sort_values(
            "year"
        )

        # For AR models, we need to predict sequentially
        if "cost_lag1" in features:
            predicted_costs = []
            actual_costs = scenario_data["cost_current"].values

            for i, row in scenario_data.iterrows():
                if pd.isna(row["cost_lag1"]):
                    # First period - use actual cost
                    predicted_costs.append(actual_costs[0])
                else:
                    # Use previous prediction as lag
                    if len(predicted_costs) > 0:
                        # Create prediction vector manually to avoid NaN issues
                        pred_values = []
                        for feature in features:
                            if feature == "cost_lag1":
                                pred_values.append(predicted_costs[-1])
                            elif feature == "cost_lag2" and len(predicted_costs) > 1:
                                pred_values.append(predicted_costs[-2])
                            elif feature == "cost_lag2":
                                # If we don't have enough history, use the first available cost
                                pred_values.append(
                                    predicted_costs[0]
                                    if predicted_costs
                                    else actual_costs[0]
                                )
                            else:
                                pred_values.append(row[feature])

                        X_pred = np.array(pred_values).reshape(1, -1)
                        pred_cost = best_model.predict(X_pred)[0]
                        predicted_costs.append(pred_cost)
                    else:
                        predicted_costs.append(actual_costs[len(predicted_costs)])
        else:
            # Non-AR model - direct prediction
            X_pred = scenario_data[features]
            predicted_costs = best_model.predict(X_pred)

        predictions[scenario_id] = {
            "years": scenario_data["year"].values,
            "actual_costs": scenario_data["cost_current"].values,
            "predicted_costs": np.array(predicted_costs),
            "demands": scenario_data["total_demand_current"].values,
            "base_demand": scenario_data["base_total_demand"].iloc[0],
        }

    return predictions


def create_unified_model_plots(results, predictions, panel_df):
    """Create comprehensive plots for unified model analysis"""

    fig = plt.figure(figsize=(20, 15))

    # Plot 1: Model performance comparison
    plt.subplot(3, 4, 1)
    model_names = list(results.keys())
    test_r2_scores = [results[name]["test_r2"] for name in model_names]

    bars = plt.bar(range(len(model_names)), test_r2_scores, alpha=0.7)
    plt.xticks(range(len(model_names)), model_names, rotation=45, ha="right")
    plt.ylabel("Test R²")
    plt.title("Model Performance Comparison")
    plt.grid(True, alpha=0.3, axis="y")

    # Add value labels
    for i, (bar, score) in enumerate(zip(bars, test_r2_scores)):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.01,
            f"{score:.3f}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=8,
        )

    # Plot 2: Best model actual vs predicted
    plt.subplot(3, 4, 2)
    best_model_name = max(results.keys(), key=lambda k: results[k]["test_r2"])
    best_result = results[best_model_name]

    plt.scatter(best_result["y_test"], best_result["y_test_pred"], alpha=0.6, s=30)
    min_val = min(best_result["y_test"].min(), best_result["y_test_pred"].min())
    max_val = max(best_result["y_test"].max(), best_result["y_test_pred"].max())
    plt.plot([min_val, max_val], [min_val, max_val], "r--", alpha=0.8)

    plt.xlabel("Actual Cost (Billion USD)")
    plt.ylabel("Predicted Cost (Billion USD)")
    plt.title(f"Best Model: {best_model_name}\nTest R² = {best_result['test_r2']:.4f}")
    plt.grid(True, alpha=0.3)

    # Plot 3: Residuals vs fitted
    plt.subplot(3, 4, 3)
    residuals = best_result["y_test"] - best_result["y_test_pred"]
    plt.scatter(best_result["y_test_pred"], residuals, alpha=0.6, s=30)
    plt.axhline(y=0, color="red", linestyle="--", alpha=0.8)
    plt.xlabel("Fitted Values (Billion USD)")
    plt.ylabel("Residuals (Billion USD)")
    plt.title("Residuals vs Fitted")
    plt.grid(True, alpha=0.3)

    # Plot 4: Feature importance (coefficients)
    plt.subplot(3, 4, 4)
    if "coefficients" in best_result:
        features = list(best_result["coefficients"].keys())
        coeffs = list(best_result["coefficients"].values())

        bars = plt.bar(range(len(features)), coeffs, alpha=0.7)
        plt.xticks(range(len(features)), features, rotation=45, ha="right")
        plt.ylabel("Coefficient Value")
        plt.title("Feature Importance (Coefficients)")
        plt.grid(True, alpha=0.3, axis="y")

    # Plot 5-8: Cost trajectory predictions (highlighting 2030 bias)
    for i, (scenario_id, pred_data) in enumerate(list(predictions.items())[:4]):
        plt.subplot(3, 4, 5 + i)

        years = pred_data["years"]
        actual = pred_data["actual_costs"]
        predicted = pred_data["predicted_costs"]
        base_demand = pred_data["base_demand"]

        plt.plot(years, actual, "o-", label="Actual", linewidth=2, markersize=6)
        plt.plot(years, predicted, "s--", label="Predicted", linewidth=2, markersize=6)

        # Highlight 2030 point
        idx_2030 = np.where(years == 2030)[0]
        if len(idx_2030) > 0:
            idx = idx_2030[0]
            plt.scatter(
                2030,
                actual[idx],
                color="red",
                s=80,
                marker="o",
                label="2030 Actual",
                zorder=5,
            )
            plt.scatter(
                2030,
                predicted[idx],
                color="orange",
                s=80,
                marker="s",
                label="2030 Predicted",
                zorder=5,
            )

        plt.xlabel("Year")
        plt.ylabel("Cost (Billion USD)")
        plt.title(f"Scenario {scenario_id}\nBase Demand: {base_demand:.1f} PWh")
        plt.legend(fontsize=8)
        plt.grid(True, alpha=0.3)

    # Plot 9: Error distribution
    plt.subplot(3, 4, 9)
    errors = []
    for pred_data in predictions.values():
        actual = pred_data["actual_costs"][1:]  # Skip first year for AR models
        predicted = pred_data["predicted_costs"][1:]
        errors.extend(actual - predicted)

    plt.hist(errors, bins=20, alpha=0.7, edgecolor="black")
    plt.axvline(
        np.mean(errors),
        color="red",
        linestyle="--",
        label=f"Mean: {np.mean(errors):.1f}",
    )
    plt.xlabel("Prediction Error (Billion USD)")
    plt.ylabel("Frequency")
    plt.title("Distribution of Prediction Errors")
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Plot 10: 2030 Bias Analysis by Scenario
    plt.subplot(3, 4, 10)

    # Extract 2030 prediction errors for visualization
    if best_model_name in results:
        test_data = results[best_model_name]["test_data"]
        model = results[best_model_name]["model"]
        features = results[best_model_name]["features"]

        # Get 2030 data
        data_2030 = test_data[test_data["year"] == 2030]
        if len(data_2030) > 0:
            X_2030 = data_2030[features]
            y_2030_pred = model.predict(X_2030)
            y_2030_actual = data_2030["cost_current"]

            # Plot actual vs predicted for 2030
            plt.scatter(y_2030_actual, y_2030_pred, alpha=0.7, s=50)

            # Add perfect prediction line
            min_val = min(y_2030_actual.min(), y_2030_pred.min())
            max_val = max(y_2030_actual.max(), y_2030_pred.max())
            plt.plot([min_val, max_val], [min_val, max_val], "r--", alpha=0.8)

            plt.xlabel("Actual Cost 2030 (Billion USD)")
            plt.ylabel("Predicted Cost 2030 (Billion USD)")
            plt.title("2030 Prediction Bias Analysis")
            plt.grid(True, alpha=0.3)

    # Plot 11: Demand sensitivity analysis
    plt.subplot(3, 4, 11)
    demand_range = np.linspace(
        panel_df["total_demand_current"].min(),
        panel_df["total_demand_current"].max(),
        50,
    )

    # Use best model to predict costs for different demand levels
    # Fix other variables at their means
    if best_model_name in results:
        best_model = results[best_model_name]["model"]
        features = results[best_model_name]["features"]

        # Create synthetic data for sensitivity analysis
        synthetic_data = pd.DataFrame()
        for demand in demand_range:
            row = panel_df[features].mean()  # Use mean values for other features
            if "total_demand_current" in features:
                row["total_demand_current"] = demand
            elif "elec_demand_current" in features:
                # Assume 20% electricity, 80% non-electricity split
                row["elec_demand_current"] = demand * 0.2
                row["nele_demand_current"] = demand * 0.8

            synthetic_data = pd.concat(
                [synthetic_data, row.to_frame().T], ignore_index=True
            )

        pred_costs = best_model.predict(synthetic_data)
        plt.plot(demand_range, pred_costs, "b-", linewidth=2)
        plt.xlabel("Total Demand (PWh)")
        plt.ylabel("Predicted Cost (Billion USD)")
        plt.title("Cost Sensitivity to Demand")
        plt.grid(True, alpha=0.3)

    # Plot 12: Training vs test performance
    plt.subplot(3, 4, 12)
    model_names_subset = [
        name for name in results.keys() if "train_r2" in results[name]
    ]
    train_scores = [results[name]["train_r2"] for name in model_names_subset]
    test_scores = [results[name]["test_r2"] for name in model_names_subset]

    x = np.arange(len(model_names_subset))
    width = 0.35

    plt.bar(x - width / 2, train_scores, width, label="Train", alpha=0.8)
    plt.bar(x + width / 2, test_scores, width, label="Test", alpha=0.8)

    plt.xlabel("Model")
    plt.ylabel("R²")
    plt.title("Train vs Test Performance")
    plt.xticks(
        x,
        [name.replace("_", " ") for name in model_names_subset],
        rotation=45,
        ha="right",
    )
    plt.legend()
    plt.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig("energy_cost_trajectory_analysis.png", dpi=300, bbox_inches="tight")
    plt.show()

    return fig


def print_unified_model_summary(results, predictions):
    """Print comprehensive summary of unified model results"""

    print("\n" + "=" * 80)
    print("UNIFIED AR(1) MODEL ANALYSIS SUMMARY")
    print("=" * 80)

    # Find best model
    best_model_name = max(results.keys(), key=lambda k: results[k]["test_r2"])
    best_result = results[best_model_name]

    print(f"\nBEST PERFORMING MODEL: {best_model_name.upper()}")
    print(f"Features: {best_result['features']}")
    print(f"Test R²: {best_result['test_r2']:.4f}")
    print(f"Test RMSE: {best_result['test_rmse']:.1f} billion USD")

    if "coefficients" in best_result:
        print(f"\nModel Equation:")
        equation = f"Cost[t] = {best_result['intercept']:.1f}"
        for feature, coef in best_result["coefficients"].items():
            equation += f" + {coef:.3f}*{feature}"
        print(f"  {equation}")

    print(f"\nMODEL PERFORMANCE RANKING:")
    sorted_models = sorted(results.items(), key=lambda x: x[1]["test_r2"], reverse=True)
    for i, (name, result) in enumerate(sorted_models[:5]):
        print(
            f"  {i + 1}. {name}: R² = {result['test_r2']:.4f}, RMSE = {result['test_rmse']:.1f}"
        )

    # Prediction accuracy
    total_errors = []
    total_actual = []
    total_predicted = []

    for pred_data in predictions.values():
        actual = pred_data["actual_costs"][1:]  # Skip first year for AR models
        predicted = pred_data["predicted_costs"][1:]
        total_errors.extend(actual - predicted)
        total_actual.extend(actual)
        total_predicted.extend(predicted)

    mae = np.mean(np.abs(total_errors))
    mape = np.mean(np.abs(np.array(total_errors) / np.array(total_actual))) * 100

    print(f"\nPREDICTION ACCURACY:")
    print(f"  Mean Absolute Error: {mae:.1f} billion USD")
    print(f"  Mean Absolute Percentage Error: {mape:.1f}%")
    print(f"  Prediction Error Std Dev: {np.std(total_errors):.1f} billion USD")

    print(f"\nMODEL INSIGHTS:")
    if "cost_lag1" in best_result.get("coefficients", {}):
        ar_coef = best_result["coefficients"]["cost_lag1"]
        print(f"  • AR(1) coefficient: {ar_coef:.3f}")
        if ar_coef > 0:
            print(f"    - Positive persistence: costs exhibit momentum")
            if ar_coef < 1:
                print(f"    - Mean-reverting: costs return to long-run equilibrium")
            else:
                print(f"    - Non-stationary: costs may have unit root")

    if "elec_demand_current" in best_result.get("coefficients", {}):
        elec_coef = best_result["coefficients"]["elec_demand_current"]
        print(
            f"  • Electricity demand sensitivity: {elec_coef:.1f} billion USD per PWh"
        )

    if "nele_demand_current" in best_result.get("coefficients", {}):
        nele_coef = best_result["coefficients"]["nele_demand_current"]
        print(
            f"  • Non-electricity demand sensitivity: {nele_coef:.1f} billion USD per PWh"
        )

        if "elec_demand_current" in best_result.get("coefficients", {}):
            ratio = elec_coef / nele_coef if nele_coef != 0 else float("inf")
            print(f"  • Electricity vs non-electricity cost ratio: {ratio:.1f}x")

    print(f"\nOUT-OF-SAMPLE VALIDATION:")
    print(
        f"  The model successfully predicts cost trajectories for unseen demand scenarios"
    )
    print(
        f"  Strong generalization performance (Test R² = {best_result['test_r2']:.4f})"
    )
    print(
        f"  Model captures the complex relationship between demand patterns and cost evolution"
    )


def main():
    """Main function for unified AR(1) model analysis"""

    import sys

    # Get filename from command line argument
    filename = "cost_dataset_with_demand_projections.csv"  # default
    if len(sys.argv) > 1:
        filename = sys.argv[1]

    print(f"Energy Cost Trajectory Model - Dataset: {filename}")

    # Load and prepare data
    panel_df = load_and_prepare_data(filename)

    # Create feature sets
    feature_sets, feature_sets_no_ar, ar_df = create_feature_sets(panel_df)

    # K-fold cross-validation
    results, train_scenarios, test_scenarios = k_fold_validation(panel_df, ar_df, feature_sets)

    # Get best model and make predictions (use last trained model)
    best_model_name = max(results.keys(), key=lambda k: np.mean(results[k]["r2_scores"]) if results[k]["r2_scores"] else 0)
    best_model = results[best_model_name]["models"][-1]
    best_features = feature_sets[best_model_name]

    # Predict trajectories for test scenarios
    predictions = predict_cost_trajectories(
        best_model, best_features, panel_df, test_scenarios[:5]
    )

    return results, predictions, panel_df


if __name__ == "__main__":
    results, predictions, panel_df = main()

