from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

if __name__ == '__main__':
    # df = pd.read_csv("data/cell_coordinates.csv")
    # df = pd.read_csv("D:\\Ryan\\197-CELL-ANNOTATOR\\public\\cell-annotator-data-v2\\cell_coordinates_with_sample.csv")
    df = pd.read_csv("data/cell_coordinates_with_sample.csv")

    # Fetch all cell_list of the sample S19_31776B1
    cell_meta_df = pd.read_csv("splitted/roi_annotated/LARC_A_S19_31776B1_annotated.csv", index_col=0)
    cell_list = cell_meta_df.index.tolist()

    x_col = "x"
    y_col = "y"

    # filter the 317 sample
    #df = df[df["cell"].isin(cell_list)]
    df = df[df["SampleId"] == "S19_31776B1"]

    fig, ax = plt.subplots()
    sns.scatterplot(data=df, x=x_col, y=y_col, s=1, ax=ax)
    plt.title("All samples")
    plt.axis("equal")
    plt.savefig("splitted/raw_scatter_v2/S317_fixed_from_ccws2.png", dpi=300)
    plt.close()
    print("Done plotting scatter plot for all samples.")