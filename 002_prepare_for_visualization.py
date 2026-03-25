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

    sample_by_cell = {}
    for file in glob(SPLITTED_META_DIR):
        patient_id = basename(file).replace(".csv", "")
        #if patient_id not in patient_response:
        #    print(f"WARNING: no clinical data for {patient_id} — skipping.")
        #    continue
        # response = patient_response[patient_id]
        cell_ids = pd.read_csv(file, usecols=[0], index_col=0).index
        for cell_id in cell_ids:
            sample_by_cell[cell_id] = patient_id

    # read the main df
    df = pd.read_csv("data/cell_coordinates.csv")
    df["SampleId"] = df["cell_id"].map(sample_by_cell)
    df.to_csv("data/cell_coordinates_with_sample.csv", index=False)
    print("Saved: data/cell_coordinates_with_sample.csv")