import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Define the colors and labels for the legend
kingdom_legend_colors = [color for color in kingdom_colors.values()]
kingdom_legend_labels = [kingdom for kingdom in kingdom_colors.keys()]

# Remove "Candidatus" from phylum labels
phylum_legend_colors = [color for color in phylum_colors.values()]
phylum_legend_labels = [phylum.replace("Candidatus ", "") for phylum in phylum_colors.keys()]

# Create a figure and axis for the legends
fig, ax = plt.subplots(1, 2, figsize=(12, 4))

# Create a legend for kingdoms
kingdom_legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(kingdom_legend_colors, kingdom_legend_labels)]
ax[0].legend(handles=kingdom_legend_patches, loc='center', fontsize=10)
ax[0].set_title("Kingdom Colors")

# Create a legend for phyla
phylum_legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(phylum_legend_colors, phylum_legend_labels)]
ax[1].legend(handles=phylum_legend_patches, loc='center', fontsize=10)
ax[1].set_title("Phylum Colors")

# Remove x and y ticks for a clean legend plot
ax[0].axis('off')
ax[1].axis('off')

# Show the legends
plt.tight_layout()
plt.show()
plt.save_fig("legends_sanky.svg") 

