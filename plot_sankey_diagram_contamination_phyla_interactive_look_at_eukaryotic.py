import pandas as pd
import plotly.graph_objects as go
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np

svg_width = 500
svg_height = 700

# Assuming df is your DataFrame
filename = "Data/gx_details_refseq.20230416_contaminations_from_prokaryotes.csv"


def interpolate_color(value, min_val, max_val, color1, color2, color3):
    # Normalize the value
    norm_value = (value - min_val) / (max_val - min_val)

    # Create a custom colormap
    cmap = LinearSegmentedColormap.from_list("custom_cmap", [color1, color2, color3])

    # Use the normalized value to get the interpolated color
    return mcolors.to_hex(cmap(norm_value))


def hex_to_rgba(hex_color, alpha):
    rgb = mcolors.to_rgb(hex_color)  # Convert hex to RGB
    return (
        f"rgba({int(rgb[0] * 255)}, {int(rgb[1] * 255)}, {int(rgb[2] * 255)}, {alpha})"
    )


df_raw = pd.read_csv(filename)
df = df_raw[df_raw["Contamination_Superkingdom"].isin(["Bacteria", "Archaea"])]
# TO only show eukayrotic hosts
df = df_raw[~df_raw["Host_Superkingdom"].isin(["Bacteria", "Archaea"])]
df["contam_length"] = df["end_pos"] - df["start_pos"]

# Use Host_Kingdom if available, otherwise use Host_Superkingdom
df["Host_Kingdom_Group"] = df["Host_Kingdom"].fillna(df["Host_Superkingdom"])

# Calculate the counts for each combination
counts = (
    df.groupby(
        [
            "Host_Kingdom_Group",
            "Host_Phylum",
            "Contamination_Phylum",
            "Contamination_Family",
        ]
    )
    .size()
    .reset_index(name="counts")
)

# Sort the data first by Host_Kingdom_Group then by Host_Phylum
counts.sort_values(by=["Host_Kingdom_Group", "Host_Phylum"], inplace=True)
print(counts["Host_Kingdom_Group"].unique())
print(counts["Contamination_Phylum"].unique())
# Define color mapping for each Host_Kingdom_Group
kingdom_colors = {
    "Archaea": "#44424D",
    "Bacteria": "#516DA4",
    "Eukaryota": "#727272",
    "Fungi": "#73668F",
    "Metazoa": "#924F60",
    "Viridiplantae": "#658F65",
}

# phylum_colors = {
#     "Bacillota": "#835e5e",
#     "Pseudomonadota": "#83645e",
#     "Thermodesulfobacteriota": "#836a5e",
#     "Acidobacteriota": "#83705e",
#     "Actinomycetota": "#83765e",
#     "Bacteroidota": "#837c5e",
#     "Candidatus Cloacimonadota": "#83825e",
#     "Deinococcota": "#7e835e",
#     "Mycoplasmatota": "#78835e",
#     "Planctomycetota": "#72835e",
#     "Synergistota": "#6c835e",
#     "Verrucomicrobiota": "#66835e",
#     "Nitrospinota": "#60835e",
#     "Campylobacterota": "#5e8362",
#     "Candidatus Saccharibacteria": "#5e8368",
#     "Cyanobacteriota": "#5e836e",
#     "Fusobacteriota": "#5e8374",
#     "Spirochaetota": "#5e837a",
#     "Balneolota": "#5e8380",
#     "Chlamydiota": "#5e8083",
#     "Ignavibacteriota": "#5e7a83",
#     "Rhodothermota": "#5e7483",
#     "Myxococcota": "#5e6e83",
#     "Armatimonadota": "#5e6883",
#     "Chloroflexota": "#5e6283",
#     "Aquificota": "#605e83",
#     "Coprothermobacterota": "#665e83",
#     "Deferribacterota": "#6c5e83",
#     "Fibrobacterota": "#725e83",
#     "Lentisphaerota": "#785e83",
#     "Nitrospirota": "#7e5e83",
#     "Thermotogota": "#835e82",
#     "Candidatus Moduliflexota": "#835e7c",
#     "Chrysiogenota": "#835e76",
#     "Chlorobiota": "#835e70",
#     "Kiritimatiellota": "#835e6a",
#     "Bdellovibrionota": "#835e64",
phylum_colors = {
    "Bacillota": "#8B5A5A",  # Darker Shade of Red
    "Pseudomonadota": "#8B6B5A",  # Earthy Orange
    "Thermodesulfobacteriota": "#8B765A",  # Muted Gold
    "Acidobacteriota": "#8B825A",  # Olive Green
    "Actinomycetota": "#8B8B5A",  # Dark Khaki
    "Bacteroidota": "#7A8B5A",  # Moss Green
    "Candidatus Cloacimonadota": "#698B5A",  # Greenish Grey
    "Deinococcota": "#5A8B5A",  # Medium Sea Green
    "Mycoplasmatota": "#5A8B6B",  # Seafoam Green
    "Planctomycetota": "#5A8B7A",  # Muted Teal
    "Synergistota": "#5A8B8B",  # Dull Cyan
    "Verrucomicrobiota": "#5A7A8B",  # Slate Blue
    "Nitrospinota": "#5A698B",  # Steel Blue
    "Campylobacterota": "#5A5A8B",  # Royal Blue
    "Candidatus Saccharibacteria": "#6B5A8B",  # Purple
    "Cyanobacteriota": "#7A5A8B",  # Lavender
    "Fusobacteriota": "#8B5A8B",  # Plum
    "Spirochaetota": "#8B5A7A",  # Pale Violet Red
    "Balneolota": "#8B5A69",  # Indian Red
    "Chlamydiota": "#8B5A5A",  # Brick Red
    "Ignavibacteriota": "#8B5A5A",  # Dark Red
    "Rhodothermota": "#8B6969",  # Rosy Brown
    "Myxococcota": "#8B7575",  # Dark Salmon
    "Armatimonadota": "#8B8282",  # Light Slate Gray
    "Chloroflexota": "#8B8B8B",  # Gray
    "Aquificota": "#82828B",  # Slate Gray
    "Coprothermobacterota": "#75758B",  # Light Steel Blue
    "Deferribacterota": "#69698B",  # Medium Slate Blue
    "Fibrobacterota": "#5A5A8B",  # Cornflower Blue
    "Lentisphaerota": "#5A698B",  # Dodger Blue
    "Nitrospirota": "#5A758B",  # Light Sky Blue
    "Thermotogota": "#5A828B",  # Deep Sky Blue
    "Candidatus Moduliflexota": "#5A8B8B",  # Cadet Blue
    "Chrysiogenota": "#5A8B75",  # Medium Aquamarine
    "Chlorobiota": "#5A8B69",  # Medium Sea Green
    "Kiritimatiellota": "#5A8B5A",  # Forest Green
    "Bdellovibrionota": "#698B5A",  # Olive Drab
}


# Create mappings for Host_Phylum and Contamination_Family to numerical ids
host_phylum_mapping = {
    phylum: i for i, phylum in enumerate(counts["Host_Phylum"].unique())
}
contam_class_mapping = {
    cls: i + len(host_phylum_mapping)
    for i, cls in enumerate(counts["Contamination_Family"].unique())
}

# Define the source, target, and value for the Sankey diagram
source = [host_phylum_mapping[phylum] for phylum in counts["Host_Phylum"]]
target = [contam_class_mapping[cls] for cls in counts["Contamination_Family"]]
value = counts["counts"].tolist()

# Calculate y-positions for each node based on their sorted order within kingdoms
y_positions = []
current_y = 0
complete_entries = counts["counts"].sum()
spacer = 0.012
height = 1 - (len(counts["Host_Phylum"].unique()) * spacer)
for kingdom in counts["Host_Kingdom_Group"].unique():
    kingdom_phyla = counts[counts["Host_Kingdom_Group"] == kingdom][
        "Host_Phylum"
    ].unique()
    for phylum in kingdom_phyla:
        phylum_index = host_phylum_mapping[phylum]
        phylum_entries = counts[
            (counts["Host_Kingdom_Group"] == kingdom)
            & (counts["Host_Phylum"] == phylum)
        ]["counts"].sum()
        phylum_space = phylum_entries / complete_entries * height + spacer
        y_positions.append(current_y + phylum_space * 0.5)
        current_y += phylum_space

# Also do this for the contamination axis
current_y = 0
complete_entries = counts["counts"].sum()
spacer = 0.008
height = 1 - (len(counts["Contamination_Family"].unique()) * spacer)
original_order = {key: "" for key in counts["Contamination_Family"].unique()}
for phylum in counts["Contamination_Phylum"].unique():
    phylum_class = counts[counts["Contamination_Phylum"] == phylum][
        "Contamination_Family"
    ].unique()
    print(phylum_class)
    for cls in phylum_class:
        cls_index = contam_class_mapping[cls]
        cls_entries = counts[
            (counts["Contamination_Phylum"] == phylum)
            & (counts["Contamination_Family"] == cls)
        ]["counts"].sum()
        cls_space = cls_entries / complete_entries * height + spacer
        original_order[cls] = current_y + cls_space * 0.5

        current_y += cls_space

y_positions.extend([value for key, value in original_order.items()])

# Apply color mappings
node_colors = [
    kingdom_colors.get(
        df[df["Host_Phylum"] == phylum]["Host_Kingdom_Group"].iloc[0], "grey"
    )
    for phylum in counts["Host_Phylum"].unique()
] + [
    phylum_colors.get(
        df[df["Contamination_Family"] == cls]["Contamination_Phylum"].iloc[0], "grey"
    )
    for cls in counts["Contamination_Family"].unique()
]

host_phylum_labels = [
    phylum.replace("Candidatus ", "") for phylum in host_phylum_mapping.keys()
]
contam_class_labels = [
    cls.replace("Candidatus ", "") for cls in contam_class_mapping.keys()
]

# Ctreate Data for the links
# Calculate the average contam_length and normalize for color mapping
counts["median_contam_length"] = (
    df.groupby(
        [
            "Host_Kingdom_Group",
            "Host_Phylum",
            "Contamination_Phylum",
            "Contamination_Family",
        ]
    )["contam_length"]
    .median()
    .reset_index(name="median_contam_length")["median_contam_length"]
)

min_median = np.percentile(counts["median_contam_length"], 1)
max_median = np.percentile(counts["median_contam_length"], 60)

counts["link_color"] = counts["median_contam_length"].apply(
    lambda x: interpolate_color(
        x, min_median, max_median, "#924F60", "#73668F", "#516DA4"
    )
)

counts["opacity"] = (counts["counts"] / counts["counts"].max()) * 0.6 + 0.4

link_colors_rgba = [
    hex_to_rgba(color, opacity)
    for color, opacity in zip(counts["link_color"], counts["opacity"])
]


# Create the Sankey diagram
fig = go.Figure(
    data=[
        go.Sankey(
            node=dict(
                pad=20,
                thickness=25,
                line=dict(color=node_colors, width=0.0),
                label=host_phylum_labels + contam_class_labels,
                color=node_colors,
                x=[0.001] * len(host_phylum_mapping)
                + [0.999] * len(contam_class_mapping),
                y=y_positions,
            ),
            link=dict(
                source=source, target=target, value=value, color=link_colors_rgba
            ),
        )
    ]
)

# Update layout for the plot
fig.update_layout(
    title_text="Host Phylum to Contamination Class Flow",
    font_size=4,
    font_family="Arial",
    font_color="black",
    annotations=[
        go.layout.Annotation(
            x=0.5,
            y=-0.1,
            showarrow=False,
            xref="paper",
            yref="paper",
        )
    ],
)


# Save the figure as an SVG file and HTML file
file_name_svg = "Host_Contamination_Flow_Phylum_to_Class_Color_just_eukaryotic.svg"
fig.write_image(file_name_svg, width=svg_width, height=svg_height)
print(f"Plot saved as {file_name_svg}")

file_name_html = "Host_Contamination_Flow_Phylum_to_Class_just_eukaryotic.html"
fig.write_html(file_name_html)
print(f"Plot saved as {file_name_html}")


# Do the legends

kingdom_legend_colors = list(kingdom_colors.values())
kingdom_legend_labels = list(kingdom_colors.keys())
phylum_legend_colors = list(phylum_colors.values())
phylum_legend_labels = [
    phylum.replace("Candidatus ", "") for phylum in phylum_colors.keys()
]

# Define colors for the gradient legend for median contam_length
colors = ["#924F60", "#73668F", "#516DA4"]
cmap = mcolors.LinearSegmentedColormap.from_list("contam_length_cmap", colors)
min_median = counts["median_contam_length"].min()
max_median = counts["median_contam_length"].max()

# Create Matplotlib figure for legends
fig, ax = plt.subplots(1, 3, figsize=(18, 4))

# Kingdom legend
kingdom_legend_patches = [
    mpatches.Patch(color=color, label=label)
    for color, label in zip(kingdom_legend_colors, kingdom_legend_labels)
]
ax[0].legend(handles=kingdom_legend_patches, loc="center", fontsize=10)
ax[0].set_title("Kingdom Colors")
ax[0].axis("off")

# Phylum legend
phylum_legend_patches = [
    mpatches.Patch(color=color, label=label)
    for color, label in zip(phylum_legend_colors, phylum_legend_labels)
]
ax[1].legend(handles=phylum_legend_patches, loc="center", fontsize=10)
ax[1].set_title("Phylum Colors")
ax[1].axis("off")

# Gradient legend for median contam_length
sm = ScalarMappable(cmap=cmap, norm=Normalize(vmin=min_median, vmax=max_median))
sm.set_array([])
cb = plt.colorbar(sm, ax=ax[2], orientation="horizontal", fraction=1)
cb.set_label("Median Contamination Length")
ax[2].set_title("Median Contamination Gradient")

plt.tight_layout()
# Save the figure with legends
fig.savefig("legends_sanky_just_eukaryotic.svg")
