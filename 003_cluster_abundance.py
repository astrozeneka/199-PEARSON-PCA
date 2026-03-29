
import pandas as pd

if __name__ == '__main__':
    relative_abundance_stacked = None

    df = pd.read_csv("data/cell_coordinates_with_sample.csv")
    samples = df["SampleId"].unique()
    for sample in samples:
        print(f"Processing {sample}...")
        df_sample = df[df["SampleId"] == sample]
        cell_type_counts = df_sample["cluster"].value_counts().sort_index()
        cell_type_percentage = cell_type_counts / cell_type_counts.sum()
        cell_type_percentage.name = sample
        if relative_abundance_stacked is None:
            relative_abundance_stacked = cell_type_percentage
        else:
            relative_abundance_stacked = pd.concat([relative_abundance_stacked, cell_type_percentage], axis=1)

    # Save the relative abundance to csv
    relative_abundance_stacked.to_csv("splitted/overall_abundance/overall_abundance.csv")
    print("done")
