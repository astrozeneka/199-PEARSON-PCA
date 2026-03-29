from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd


if __name__ == '__main__':
    #df = pd.read_csv("splitted/roi_annotated/LARC_A_S19_31776B1_annotated.csv")
    df = pd.read_csv("data/cell_coordinates_S19_31776B1.csv")

    x_col = "x_slide_mm" if "x_slide_mm" in df.columns else "x"
    y_col = "y_slide_mm" if "y_slide_mm" in df.columns else "y"

    # ONly slect with selection_mask is 1
    # df = df[df["selection_mask"] == 1]
    # do a scatterplot using seaborn and
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.scatterplot(data=df, x=x_col, y=y_col, s=1, ax=ax)
    plt.title("Sample: LARC_A_S19_31776B1")
    plt.axis("equal")
    plt.savefig("splitted/raw_scatter_v2/LARC_A_S19_31776B1_v2.png", dpi=300)
    plt.close()
    print("Done plotting scatter plot for LARC_A_S19_31776B1.")