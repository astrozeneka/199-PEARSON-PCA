from glob import glob
from os.path import basename
import pandas as pd

CLINICAL_DATA_PATH = '../125-WHOLE-SLIDE-ANALYSIS/density_with_cd/Clinical_data.csv'
SPLITTED_META_DIR = "../193-SEURAT-HPC/splitted/meta/*.csv"

if __name__ == '__main__':

    clinical_df = pd.read_csv(CLINICAL_DATA_PATH)
    clinical_df = clinical_df.dropna(subset=['SampleId'])
    clinical_df.set_index('SampleId', inplace=True)
    patient_response = clinical_df['pCR'].to_dict()

    files = glob(SPLITTED_META_DIR)
    meta_frames = [
        pd.read_csv(f, usecols=[0], index_col=0).assign(SampleId=basename(f).replace(".csv", ""))
        for f in files
    ]
    cell_to_sample = pd.concat(meta_frames)["SampleId"].to_dict()

    # read the main df
    df = pd.read_csv("data/cell_coordinates.csv")
    df["SampleId"] = df["cell"].map(cell_to_sample)
    df.to_csv("data/cell_coordinates_with_sample.csv", index=False)
    print("Saved: data/cell_coordinates_with_sample.csv")