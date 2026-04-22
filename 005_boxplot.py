import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

CLINICAL_DATA_PATH = '../125-WHOLE-SLIDE-ANALYSIS/density_with_cd/Clinical_data.csv'
ABUNDANCE_PATH = 'splitted/overall_abundance/overall_abundance.csv'
CLUSTER = 8
OUTPUT_PATH = 'splitted/overall_abundance/boxplot_cluster8.png'


if __name__ == '__main__':

    clinical_df = pd.read_csv(CLINICAL_DATA_PATH)
    clinical_df = clinical_df.dropna(subset=['SampleId'])
    clinical_df.set_index('SampleId', inplace=True)
    patients_response_info = clinical_df['pCR'].to_dict()

    abundance_df = pd.read_csv(ABUNDANCE_PATH, index_col='cluster')
    abundance_df = abundance_df.drop(
        columns=[c for c in abundance_df.columns if c == 'UNIDENTIFIED' or c.strip() == '']
    )

    responders = [pid for pid in abundance_df.columns
                  if patients_response_info.get(pid) == 'Responder']
    non_responders = [pid for pid in abundance_df.columns
                      if patients_response_info.get(pid) == 'Non-responder']

    r_vals = abundance_df.loc[CLUSTER, responders].astype(float).values
    nr_vals = abundance_df.loc[CLUSTER, non_responders].astype(float).values

    rows = (
        [{'group': 'Responder', 'abundance': v} for v in r_vals] +
        [{'group': 'Non-responder', 'abundance': v} for v in nr_vals]
    )
    long_df = pd.DataFrame(rows)

    group_order = ['Responder', 'Non-responder']
    palette = {'Responder': '#009e73', 'Non-responder': '#ed1558'}

    fig, ax = plt.subplots(figsize=(4, 5))

    sns.boxplot(data=long_df, x='group', y='abundance', order=group_order,
                palette=palette, showfliers=False, width=0.5, ax=ax)
    sns.stripplot(data=long_df, x='group', y='abundance', order=group_order,
                  color='black', alpha=0.7, size=5, jitter=True, ax=ax)

    # Mann-Whitney p-value
    _, pval = mannwhitneyu(r_vals, nr_vals, alternative='two-sided')

    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    bracket_y = y_max + y_range * 0.05
    ax.plot([0, 0, 1, 1],
            [bracket_y, bracket_y + y_range * 0.02, bracket_y + y_range * 0.02, bracket_y],
            'k-', linewidth=0.8)
    label = f"p={pval:.2e}" if pval < 0.001 else f"p={pval:.3f}"
    ax.text(0.5, bracket_y + y_range * 0.04, label, ha='center', fontsize=8,
            color='red' if pval < 0.05 else 'black')
    ax.set_ylim(y_min, bracket_y + y_range * 0.12)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('')
    ax.set_ylabel('Relative abundance')
    ax.set_title(f'Cluster {CLUSTER} abundance')

    plt.tight_layout()
    plt.savefig(OUTPUT_PATH, dpi=330, transparent=False)
    print(f"Saved plot to {OUTPUT_PATH}")
