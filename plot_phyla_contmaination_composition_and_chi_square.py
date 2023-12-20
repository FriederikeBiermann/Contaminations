import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
import matplotlib.patches as mpatches

# Load your data
filename = "Data/gx_details_refseq.20230416_contaminations_from_prokaryotes.csv"
df = pd.read_csv(filename)


# Define color mappings
kingdom_colors = {
    "Archaea": "#3b3a42",  # Darker and more saturated
    "Bacteria": "#455a91",  # Darker and more saturated
    "Eukaryota": "#5e5e5e",  # Darker and more saturated
    "Fungi": "#625c78",  # Darker and more saturated
    "Metazoa": "#7d3e51",  # Darker and more saturated
    "Viridiplantae": "#557455",  # Darker and more saturated
}


contamination_phylum_colors = {
    "Bacillota": "#625c78",  # Crimson
    "Pseudomonadota": "#557455",  # Dark Blue
    "Actinomycetota": "#7d3e51",  # Green
    "Bacteroidota": "#455a91",  # Amber
}

df["Host_Kingdom_Group"] = df["Host_Kingdom"].fillna(df["Host_Superkingdom"])
# Map Host_Phylum to Host_Kingdom_Group and Contamination_Genus to Contamination_Phylum
df["Host_Kingdom_Color"] = df["Host_Kingdom_Group"].map(kingdom_colors)
df["Contamination_Phylum_Color"] = df["Contamination_Phylum"].map(
    contamination_phylum_colors
)

host_phylum_to_kingdom = df.drop_duplicates(subset="Host_Phylum").set_index(
    "Host_Phylum"
)["Host_Kingdom_Group"]
host_phylum_colors = host_phylum_to_kingdom.map(kingdom_colors)

Contamination_Genus_to_phylum = df.drop_duplicates(
    subset="Contamination_Genus"
).set_index("Contamination_Genus")["Contamination_Phylum"]
Contamination_Genus_colors = Contamination_Genus_to_phylum.map(
    contamination_phylum_colors
)

# Select the relevant columns
data = df[["Host_Phylum", "Contamination_Genus"]]

# Calculate the overall percentage of each Contamination_Genus
total_counts = data["Contamination_Genus"].value_counts()
total_percentage = total_counts / total_counts.sum()

# Identify the top 10 families based on overall percentage
top_families = total_percentage.nlargest(30).index
print([Contamination_Genus_to_phylum[genus] for genus in top_families])

# Filter the DataFrame to include only these top 10 families
filtered_data = data[data["Contamination_Genus"].isin(top_families)]

# Create a contingency table with the filtered data
contingency_table = pd.crosstab(
    filtered_data["Host_Phylum"], filtered_data["Contamination_Genus"]
)

# Perform the Chi-Square Test
chi2, p, dof, expected = chi2_contingency(contingency_table)
print(f"Chi-Square Statistic: {chi2}, P-value: {p}")

# Convert to relative values (percentages) within each Host_Phylum
contingency_table_percentage = contingency_table.div(
    contingency_table.sum(axis=1), axis=0
)

# Plotting the heatmap
plt.figure(figsize=(14, 10))

# Define custom diverging colormap
cmap = sns.diverging_palette(240, 10, n=7, center="light", as_cmap=True)

# Create the heatmap
ax = sns.heatmap(
    contingency_table_percentage, cmap=cmap, linewidths=0.5, linecolor="grey"
)
# Modify y-tick labels
yticklabels = [
    label.get_text().replace("Candidatus ", "") for label in plt.gca().get_yticklabels()
]
plt.gca().set_yticklabels(yticklabels)

# Modify x-tick labels
xticklabels = [
    label.get_text().replace("Candidatus ", "") for label in plt.gca().get_xticklabels()
]
plt.gca().set_xticklabels(xticklabels)

# Color the y-tick labels by Host_Kingdom_Group
for label in plt.gca().get_yticklabels():
    label.set_color(host_phylum_colors.get(label.get_text(), "#455a91"))

# Color the x-tick labels by Contamination_Phylum
for label in plt.gca().get_xticklabels():
    label.set_color(Contamination_Genus_colors.get(label.get_text(), "black"))


# Set plot title and labels
plt.title(f"Heatmap of host phyla and top 30 contamination genera")
plt.xlabel("Contamination Family")
plt.ylabel("Host Phylum")

# Create legends
kingdom_patches = [mpatches.Patch(color=color, label=kingdom) for kingdom, color in kingdom_colors.items()]
phylum_patches = [mpatches.Patch(color=color, label=phylum) for phylum, color in contamination_phylum_colors.items()]

# Add legends to the plot
plt.legend(handles=kingdom_patches, title="Host Kingdoms", loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
plt.gca().add_artist(plt.legend(handles=phylum_patches, title="Contamination Phyla", loc='lower left', bbox_to_anchor=(1.05, 0), borderaxespad=0.))


plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig("heatmap_phyla_genera.svg")
