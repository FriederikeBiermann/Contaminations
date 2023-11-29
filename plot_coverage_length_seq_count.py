import plotly.graph_objs as go
import pandas as pd
import numpy as np

# Load your data
data = pd.read_csv("Data/gx_details_genbank.20230416_long_contaminations_from_prokaryotes_coverage_less_100.csv")

# Define ranges for minimum sequence length (log scale) and maximum coverage (10 steps)
min_lengths = np.logspace(np.log10(data['length_contamination'].min()), np.log10(data['length_contamination'].max()), 50)
max_coverages = np.linspace(data['coverage'].min(), data['coverage'].max(), 10)

# Create a Plotly figure
fig = go.Figure()

for max_cov in max_coverages:
    counts = [data[(data['length_contamination'] >= min_len) & (data['coverage'] <= max_cov)].shape[0] for min_len in min_lengths]
    
    # Add a trace for each coverage
    fig.add_trace(go.Scatter(
        x=min_lengths, 
        y=counts, 
        mode='lines+markers',
        name=f'Coverage ≤ {max_cov:.2f}',
        text=[f'Min Length: {min_len:.2f}, Coverage ≤ {max_cov:.2f}, Rows: {count}' for min_len, count in zip(min_lengths, counts)],
        hoverinfo='text'
    ))

# Update layout for log scale and titles
fig.update_layout(
    title='Number of Rows vs. Min Sequence Length for Various Max Coverages',
    xaxis_title='Min Sequence Length',
    yaxis_title='Number of Rows',
    xaxis_type='log',
    yaxis_type='log'
)

# Show the figure

fig.write_html("gx_details_genbank.20230416_long_contaminations_from_prokaryotes_coverage_less_100_length_coverage_sequences.html")
fig.show()
