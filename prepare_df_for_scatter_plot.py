import pandas as pd
data = (
    "Data/gx_details_genbank.20230416_contaminations_from_prokaryotes_with_antismash.csv"
)
df = pd.read_csv(data)
df = df.dropna(subset=["Host_Phylum", "Contamination_Family"])

# filter only ones with BGCs
df = df[df["Product_class"].notna()]
df.to_csv("Data/gx_details_genbank.20230416_contaminations_from_prokaryotes_with_antismash_80_coverage_only_hits.csv")
