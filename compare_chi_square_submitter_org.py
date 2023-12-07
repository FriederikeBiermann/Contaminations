import pandas as pd
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt

# Load your data
filename = "Data/gx_details_refseq.20230416_contaminations_from_prokaryotes.csv"
df = pd.read_csv(filename)

# Get unique lists of families, phyla, and organizations
families = df["Contamination_Genus"].unique()
phyla = df["Host_Phylum"].unique()
organizations = df["submitter_org"].unique()

pseudomonas_data = df[df["Contamination_Genus"] == "Pseudomonas"]

# Create contingency tables
cont_table_org = pd.crosstab(
    pseudomonas_data["submitter_org"], pseudomonas_data["Contamination_Genus"]
)
cont_table_phylum = pd.crosstab(
    pseudomonas_data["Host_Phylum"], pseudomonas_data["Contamination_Genus"]
)
print(cont_table_org)

# Perform chi-squared tests
chi2_org, p_org, _, _ = (
    chi2_contingency(cont_table_org) if cont_table_org.shape[0] > 1 else (0, 1, 0, 0)
)
chi2_phylum, p_phylum, _, _ = (
    chi2_contingency(cont_table_phylum)
    if cont_table_phylum.shape[0] > 1
    else (0, 1, 0, 0)
)

print(f"Chi-Squared for Pseudomonas - Organization: {chi2_org}, p-value: {p_org}")
print(f"Chi-Squared for Pseudomonas - Host Phylum: {chi2_phylum}, p-value: {p_phylum}")

# Prepare a dictionary to hold the results
chi_squared_results = {}

# Iterate over each family to perform chi-squared tests
for family in families:
    # Filter data for the current family
    family_data = df[df["Contamination_Genus"] == family]

    # Create contingency tables
    cont_table_org = pd.crosstab(
        family_data["submitter_org"], family_data["Contamination_Genus"]
    )
    cont_table_phylum = pd.crosstab(
        family_data["Host_Phylum"], family_data["Contamination_Genus"]
    )

    # Perform chi-squared tests
    chi2_org, p_org, _, _ = (
        chi2_contingency(cont_table_org)
        if cont_table_org.shape[0] > 1
        else (0, 1, 0, 0)
    )
    chi2_phylum, p_phylum, _, _ = (
        chi2_contingency(cont_table_phylum)
        if cont_table_phylum.shape[0] > 1
        else (0, 1, 0, 0)
    )

    # Store results
    chi_squared_results[family] = {
        "org_chi2": chi2_org,
        "org_p": p_org,
        "phylum_chi2": chi2_phylum,
        "phylum_p": p_phylum,
    }

# Calculate the difference in Chi-Squared values and sort
difference_in_chi_squared = {
    family: abs(info["org_chi2"] - info["phylum_chi2"])
    for family, info in chi_squared_results.items()
}
sorted_families = sorted(
    difference_in_chi_squared, key=difference_in_chi_squared.get, reverse=True
)

# Select top N families with the biggest differences
top_n = 5  # Adjust this number as needed
top_families = sorted_families[:top_n]
top_differences = [difference_in_chi_squared[family] for family in top_families]

print(top_families, top_differences)
### Step 2: Plot These Families as a Bar Graph
plt.figure(figsize=(10, 6))
plt.bar(top_families, top_differences, color="skyblue")
plt.xlabel("Contamination Genus")
plt.ylabel("Difference in Chi-Squared Statistics")
plt.title(f"Top {top_n} Contamination Genera with Largest Differences in Associations")
plt.xticks(rotation=45)
plt.show()

# Extract Chi-Squared values associated with Host_Phylum and sort
phylum_chi_squared = {
    family: info["phylum_chi2"] for family, info in chi_squared_results.items()
}
sorted_by_phylum_association = sorted(phylum_chi_squared, key=phylum_chi_squared.get)

# Select families with the lowest association (smallest Chi-Squared values)
bottom_n = 5  # Adjust this number as needed
least_dependent_families = sorted_by_phylum_association[:bottom_n]
lowest_chi_squared = [phylum_chi_squared[family] for family in least_dependent_families]
print(least_dependent_families, lowest_chi_squared)
### Step 2: Plot These Families as a Bar Graph
plt.figure(figsize=(10, 6))
plt.bar(least_dependent_families, lowest_chi_squared, color="lightgreen")
plt.xlabel("Contamination Family")
plt.ylabel("Chi-Squared with Host Genus")
plt.title(f"Top {bottom_n} Contamination Genera with Least Dependence on Host Phylum")
plt.xticks(rotation=45)
plt.show()
