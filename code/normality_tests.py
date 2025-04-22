import pandas as pd

def run(data: pd.DataFrame, column_name: str):   
    # Import necessary libraries
    from scipy.stats import shapiro, normaltest

    # Remove potential NaN values to avoid errors in calculations
    data = data.dropna()

    # Normality Test for Column (Shapiro-Wilk & D'Agostino-Pearson)
    shapiro_test = shapiro(data[column_name])
    dagostino_test = normaltest(data[column_name])

    # Display results
    print("Normality Test Results")
    print(f"\n {column_name}:")
    print(f"  - Shapiro-Wilk test: W = {shapiro_test.statistic:.4f}, p-value = {shapiro_test.pvalue:.4g}")
    print(f"  - D'Agostino-Pearson test: K^2 = {dagostino_test.statistic:.4f}, p-value = {dagostino_test.pvalue:.4g}")

    # Interpretation
    alpha = 0.05
    if shapiro_test.pvalue < alpha or dagostino_test.pvalue < alpha:
        print(f"\n The {column_name} data is NOT normally distributed (p < 0.05).")
    else:
        print(f"\n The {column_name} data appears to be normally distributed (p > 0.05).")
    print('\n')

def run_for_all(data):   
    # Import necessary libraries
    from scipy.stats import shapiro, normaltest

    # Remove potential NaN values to avoid errors in calculations
    data = data.dropna()

    # Normality Test for Observed Mutations (Shapiro-Wilk & D'Agostino-Pearson)
    shapiro_test = shapiro(data["M. tuber"])
    dagostino_test = normaltest(data["M. tuber"])

    # Normality Test for Mutation Ways
    shapiro_test_ways = shapiro(data["Mutation_Pathways"])
    dagostino_test_ways = normaltest(data["Mutation_Pathways"])

    # Normality Test for Observed/Frequency
    shapiro_test_freq = shapiro(data["LeftSide"])
    dagostino_test_freq = normaltest(data["LeftSide"])

    # Display results
    print("Normality Test Results")
    print("\n Observed Mutations (M. tuber):")
    print(f"  - Shapiro-Wilk test: W = {shapiro_test.statistic:.4f}, p-value = {shapiro_test.pvalue:.4g}")
    print(f"  - D'Agostino-Pearson test: K^2 = {dagostino_test.statistic:.4f}, p-value = {dagostino_test.pvalue:.4g}")

    print("\n Possible Mutation Ways (Mutation_Pathways):")
    print(f"  - Shapiro-Wilk test: W = {shapiro_test_ways.statistic:.4f}, p-value = {shapiro_test_ways.pvalue:.4g}")
    print(f"  - D'Agostino-Pearson test: K^2 = {dagostino_test_ways.statistic:.4f}, p-value = {dagostino_test_ways.pvalue:.4g}")

    print("\n Observed mutations / Frequency:")
    print(f"  - Shapiro-Wilk test: W = {shapiro_test_freq.statistic:.4f}, p-value = {shapiro_test_freq.pvalue:.4g}")
    print(f"  - D'Agostino-Pearson test: K^2 = {dagostino_test_freq.statistic:.4f}, p-value = {dagostino_test_freq.pvalue:.4g}")

    # Interpretation
    alpha = 0.05
    if shapiro_test.pvalue < alpha or dagostino_test.pvalue < alpha:
        print("\n The observed mutations data is NOT normally distributed (p < 0.05).")
    else:
        print("\n The observed mutations data appears to be normally distributed (p > 0.05).")

    if shapiro_test_ways.pvalue < alpha or dagostino_test_ways.pvalue < alpha:
        print("\n The possible mutation ways data is NOT normally distributed (p < 0.05).")
    else:
        print("\n The possible mutation ways data appears to be normally distributed (p > 0.05).")

    if shapiro_test_freq.pvalue < alpha or dagostino_test_freq.pvalue < alpha:
        print("\n The observed mutations/frequency data is NOT normally distributed (p < 0.05).")
    else:
        print("\n The observed mutations/frequency data appears to be normally distributed (p > 0.05).")
