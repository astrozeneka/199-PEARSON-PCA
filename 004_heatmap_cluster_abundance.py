import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import os
import argparse

from plot_utils import plot_density_stats

CLINICAL_DATA_PATH = '../125-WHOLE-SLIDE-ANALYSIS/density_with_cd/Clinical_data.csv'
OUTPUT_DIR = 'splitted/overall_abundance'

parser = argparse.ArgumentParser(description='Plot cell type abundance heatmap.')
parser.add_argument('--normalization', type=str, default='zscore',
                    help="Normalization method: 'zscore', 'minmax', 'robust', 'l1', 'l2'")
parser.add_argument('--log_scale', action='store_true',
                    help='Apply log transformation to normalized values')
args = parser.parse_args()



# --- Clinical metadata ---

def categorize_age(age):
    if age < 60:
        return '<60'
    elif age <= 69:
        return '60-69'
    return '70+'

if __name__ == '__main__':

    clinical_df = pd.read_csv(CLINICAL_DATA_PATH)
    clinical_df = clinical_df.dropna(subset=['SampleId'])
    clinical_df.set_index('SampleId', inplace=True)

    patients_sex_info = clinical_df['Sex'].map({'Male': 'M', 'Female': 'F'}).to_dict()
    patients_age_info = clinical_df['Age_Dx'].apply(categorize_age).to_dict()
    patients_response_info = clinical_df['pCR'].to_dict()
    patients_label_info = clinical_df['Patients'].str.replace('_', '', regex=False).to_dict()
    ctnm_label_info = clinical_df['cTNM'].apply(lambda x: str(int(x))).to_dict()
    location_label_info = clinical_df['cT_Locate'].to_dict()


    # --- Per-patient cell type relative abundances (pre-computed) ---

    abundance_df = pd.read_csv(f'{OUTPUT_DIR}/overall_abundance.csv', index_col='cluster')
    abundance_df = abundance_df.drop(
        columns=[c for c in abundance_df.columns if c == 'UNIDENTIFIED' or c.strip() == '']
    )

    # Warn about patients without clinical data
    known = set(patients_response_info.keys())
    unmatched = [pid for pid in abundance_df.columns if pid not in known]
    if unmatched:
        print(f"WARNING: no clinical data for {unmatched} — excluded from heatmap.")

    responders = [pid for pid in abundance_df.columns
                  if patients_response_info.get(pid) == 'Responder']
    non_responders = [pid for pid in abundance_df.columns
                      if patients_response_info.get(pid) == 'Non-responder']

    print(f"Responders    ({len(responders)}): {responders}")
    print(f"Non-responders ({len(non_responders)}): {non_responders}")


    # --- Build stats DataFrame expected by plot_density_stats ---
    #
    # Rows correspond to cell types (Combination).
    # 'Responder dict' / 'Non-responder dict' are string repr of {patient_id: abundance}.
    # P-value is from a two-sided Mann-Whitney U test.

    records = []
    for cell_type in abundance_df.index:
        r_vals = np.nan_to_num(abundance_df.loc[cell_type, responders].values.astype(float), nan=0.0)
        nr_vals = np.nan_to_num(abundance_df.loc[cell_type, non_responders].values.astype(float), nan=0.0)
        r_dict = {pid: float(v) for pid, v in zip(responders, r_vals)}
        nr_dict = {pid: float(v) for pid, v in zip(non_responders, nr_vals)}
        try:
            _, p_val = mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
        except ValueError:
            p_val = float('nan')
        records.append({
            'Combination': cell_type,
            'Responder dict': str(r_dict),
            'Non-responder dict': str(nr_dict),
            'P-value': p_val,
        })

    stats_df = pd.DataFrame(records)


    # --- Plot ---

    fig, ax = plt.subplots(figsize=(14, 10))

    plot_density_stats(
        stats_df, ax,
        title='Cell type relative abundance',
        patient_sex=patients_sex_info,
        patient_age=patients_age_info,
        patient_response=patients_response_info,
        patient_ctnm=ctnm_label_info,
        patient_location=location_label_info,
        patient_labels=patients_label_info,
        display_legends=True,
        normalization=args.normalization,
        log_scale=args.log_scale,
        manual_group_patterns=None,  # cluster all rows symmetrically; no predefined ordering
        map_patient_id=False,
    )

    plt.tight_layout()
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    scale_suffix = 'log' if args.log_scale else 'raw'
    out_path = f"{OUTPUT_DIR}/abundance_heatmap_{args.normalization}_{scale_suffix}.png"
    plt.savefig(out_path, dpi=330, transparent=False, bbox_inches='tight')
    print(f"Plot saved to {out_path}")
