from glob import glob
from os.path import basename
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

CLINICAL_DATA_PATH = '../125-WHOLE-SLIDE-ANALYSIS/density_with_cd/Clinical_data.csv'
SPLITTED_META_DIR = "../193-SEURAT-HPC/splitted/meta/*.csv"

if __name__ == '__main__':

    files = glob(SPLITTED_META_DIR)
    meta_frames = [
        pd.read_csv(f, usecols=[0], index_col=0).assign(SampleId=basename(f).replace(".csv", ""))
        for f in files
    ]
    cell_to_sample = pd.concat(meta_frames)["SampleId"].to_dict()

    df = pd.read_csv("data/cell_coordinates.csv")
    df["SampleId"] = df["cell"].map(cell_to_sample)
    for f in files:
        slug = basename(f).replace(".csv", "")
        sub_df = df[df["SampleId"] == slug]
        # Do a scatter plot of x and
        plt.figure(figsize=(6, 6))
        sns.scatterplot(data=sub_df, x="x", y="y", s=1)
        plt.title(f"Sample: {slug}")
        plt.axis("equal")
        plt.savefig(f"splitted/raw_scatter_fixed/{slug}.png", dpi=300)
        plt.close()
    print("Done plotting scatter plots for each sample.")
