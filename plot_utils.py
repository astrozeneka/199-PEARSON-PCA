import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, leaves_list


def zscore_with_nan(row):
    """Apply zscore while preserving NaN values"""
    mask = ~np.isnan(row)
    if mask.sum() < 2:  # Need at least 2 values for zscore
        return row * 0  # Return zeros if insufficient data
    result = row.copy()
    valid_values = row[mask]
    if valid_values.std() == 0:  # No variation
        result[mask] = 0
    else:
        result[mask] = (valid_values - valid_values.mean()) / valid_values.std()
    return result


def plot_density_stats(df, ax, title="",
                       patient_sex=None, patient_age=None, patient_response=None,
                       patient_ctnm=None, patient_location=None,
                       patient_labels=None, display_legends=True,
                       xtick_fontsize=11, normalization='zscore', log_scale=False,
                       manual_group_patterns=None, map_patient_id=True):
    """
    Plot density/abundance statistics as a heatmap with optional annotations.

    Parameters:
    -----------
    df : DataFrame
        Input data with columns: 'Combination', 'Responder dict', 'Non-responder dict', 'P-value'.
        'Responder dict' and 'Non-responder dict' are string representations of
        {patient_id: value} dicts.
    ax : matplotlib axis
        Axis to plot on.
    title : str
        Plot title.
    patient_sex : dict, optional
        {patient_id: 'M' or 'F'}
    patient_age : dict, optional
        {patient_id: '<60', '60-69', or '70+'}
    patient_response : dict, optional
        {patient_id: 'Responder' or 'Non-responder'}
    patient_ctnm : dict, optional
        {patient_id: '2' or '3'}
    patient_location : dict, optional
        {patient_id: 'Lower rectum', 'Mid rectum', or 'Upper rectum'}
    patient_labels : dict, optional
        {patient_id: display_label} — replaces x-tick labels.
    manual_group_patterns : list of str, optional
        Ordered list of row names used to define the primary group for column clustering:
        responder columns are clustered using these rows, non-responder columns are clustered
        using the remaining rows. When None, all rows are used for clustering both groups
        and no predefined row ordering is applied.
    map_patient_id : bool, optional (default True)
        When True, x-tick labels are replaced using patient_labels.
        When False, the raw sample IDs extracted from column names are shown instead.
        Color coding and group brackets are applied either way.
    """
    group_prefixes = {'group1': 'R_', 'group2': 'NR_'}
    label_colors = {'group1': 'green', 'group2': 'red'}
    legend_labels = {'group1': 'Responder', 'group2': 'Non-Responder'}

    # Annotation colors
    sex_colors = {'M': 'steelblue', 'F': 'lightcoral'}
    age_colors = {'<60': '#FFB000', '60-69': '#FE6100', '70+': '#DC267F'}
    ctnm_colors = {'2': '#FFB000', '3': '#f23349'}
    location_colors = {'Lower rectum': '#9c20e8', 'Mid rectum': '#FE6100', 'Upper rectum': '#FFB000'}
    response_colors = {'Responder': '#009E73', 'Non-responder': '#ed1558'}

    # Prepare the region_df
    region_df = []
    for _, row in df.iterrows():
        responder_dict = eval(row['Responder dict'])
        non_responder_dict = eval(row['Non-responder dict'])
        region_df.append({
            'Combination': row['Combination'],
            **{f"{group_prefixes['group1']}{k}": v for k, v in responder_dict.items()},
            **{f"{group_prefixes['group2']}{k}": v for k, v in non_responder_dict.items()},
            'P-value': row['P-value']
        })
    region_df = pd.DataFrame(region_df)

    # Extract p-values and combination names
    p_values = region_df['P-value'].values

    # Prepare heatmap data
    heatmap_data = region_df.set_index('Combination').drop(columns=['P-value'])

    # Normalization
    if normalization == 'zscore':
        normalized = heatmap_data.apply(zscore_with_nan, axis=1)
    elif normalization == 'minmax':
        lower = heatmap_data.min(axis=1)
        upper = heatmap_data.max(axis=1)
        normalized = (heatmap_data.sub(lower, axis=0)).div(upper - lower, axis=0)
    elif normalization == 'robust':
        q1 = heatmap_data.quantile(0.25, axis=1)
        q3 = heatmap_data.quantile(0.75, axis=1)
        iqr = q3 - q1
        normalized = heatmap_data.sub(q1, axis=0).div(iqr, axis=0)
    elif normalization == 'l1':
        row_sums = heatmap_data.abs().sum(axis=1)
        normalized = heatmap_data.div(row_sums, axis=0)
    elif normalization == 'l2':
        row_squares = (heatmap_data ** 2).sum(axis=1)
        normalized = heatmap_data.div(np.sqrt(row_squares), axis=0)
    else:
        raise ValueError(f"Unsupported normalization method '{normalization}', "
                         f"supported: 'zscore', 'minmax', 'robust', 'l1', 'l2'")

    # Apply log transformation if requested
    if log_scale:
        min_val = normalized.min().min()
        if min_val <= 0:
            normalized = normalized - min_val + 1e-10
        normalized = np.log1p(normalized)

    # Calculate symmetric colormap bounds
    max_abs_value = max(abs(normalized.min().min()), abs(normalized.max().max()))
    vmin, vmax = -max_abs_value, max_abs_value

    # Split columns into Responder (R_) and Non-Responder (NR_) groups
    r_cols = [col for col in normalized.columns if col.startswith(group_prefixes['group1'])]
    nr_cols = [col for col in normalized.columns if col.startswith(group_prefixes['group2'])]

    # Hierarchical clustering — dual-axis: R cols clustered independently from NR cols.
    # When manual_group_patterns is provided, R cols are clustered using only those rows
    # and NR cols are clustered using the remaining rows (original block-specific logic).
    # When None, all rows are used for clustering both groups.
    if manual_group_patterns is not None:
        manual_group_mask = normalized.index.isin(manual_group_patterns)
        r_cluster_data = normalized.loc[manual_group_mask, r_cols]
        nr_cluster_data = normalized.loc[~manual_group_mask, nr_cols]
    else:
        r_cluster_data = normalized[r_cols]
        nr_cluster_data = normalized[nr_cols]

    if len(r_cluster_data) > 0 and len(r_cols) > 1:
        r_col_linkage = linkage(r_cluster_data.T.fillna(0), method='ward')
        r_col_order = [r_cols[i] for i in leaves_list(r_col_linkage)]
    else:
        r_col_order = r_cols

    if len(nr_cluster_data) > 0 and len(nr_cols) > 1:
        nr_col_linkage = linkage(nr_cluster_data.T.fillna(0), method='ward')
        nr_col_order = [nr_cols[i] for i in leaves_list(nr_col_linkage)]
    else:
        nr_col_order = nr_cols

    # Reorder rows and columns
    ordered_cols = r_col_order + nr_col_order
    if manual_group_patterns is not None:
        ordered_rows = [r for r in manual_group_patterns if r in normalized.index]
        remaining = [r for r in normalized.index if r not in ordered_rows]
        normalized = normalized.loc[ordered_rows + remaining, ordered_cols]
    else:
        normalized = normalized[ordered_cols]

    # Reorder p-values to match row order
    p_values_dict = dict(zip(region_df['Combination'], p_values))
    p_values_clustered = [p_values_dict[combo] for combo in normalized.index]

    # Custom colormap
    custom_cmap = LinearSegmentedColormap.from_list('custom_diverging',
                                                    ['#0071ff', '#ffffff', '#ff004f'])

    sns.heatmap(normalized, cmap=custom_cmap, center=0, annot=False,
                vmin=vmin, vmax=vmax,
                cbar_kws={'label': 'Z-score', 'shrink': 0.2, 'pad': 0.15, 'aspect': 10,
                          'anchor': (0.0, 0.08), 'location': 'left'},
                ax=ax)
    ax.set_title(title, pad=1)
    ax.set_ylabel("")

    # Replace underscores/dots with spaces in ytick labels
    ax.set_yticklabels([label.get_text().replace('_', ' ').replace('.', ' ')
                        for label in ax.get_yticklabels()])

    # Update xtick labels if patient_labels provided
    if patient_labels is not None:
        column_names = normalized.columns
        new_labels = []
        for col_name in column_names:
            patient_id = col_name.split('_', 1)[1] if '_' in col_name else col_name
            label = patient_labels.get(patient_id, patient_id) if map_patient_id else patient_id
            new_labels.append(label)
        ax.set_xticklabels(new_labels, rotation=45, ha='right', fontsize=xtick_fontsize)

        # Color x-axis labels by group
        labels = ax.get_xticklabels()
        for i, label in enumerate(labels):
            col_name = column_names[i]
            if col_name.startswith(group_prefixes['group1']):
                label.set_color(response_colors['Responder'])
            elif col_name.startswith(group_prefixes['group2']):
                label.set_color(response_colors['Non-responder'])

        # Add brackets for response groups
        responder_indices = [i for i, col in enumerate(column_names)
                             if col.startswith(group_prefixes['group1'])]
        nonresponder_indices = [i for i, col in enumerate(column_names)
                                if col.startswith(group_prefixes['group2'])]

        bracket_offset_y = -0.13
        bracket_offset = -0.6
        bracket_gap = 0.2
        bracket_height = 0.03
        text_offset = 0.02

        if responder_indices:
            x_start = responder_indices[0] - 0.5 + bracket_offset
            x_end = responder_indices[-1] + 0.5 + bracket_offset - bracket_gap / 2
            ax.plot([x_start, x_start, x_end, x_end],
                    [bracket_offset_y, bracket_offset_y - bracket_height,
                     bracket_offset_y - bracket_height, bracket_offset_y],
                    color='black', lw=1.5, clip_on=False,
                    transform=ax.get_xaxis_transform())
            ax.text((x_start + x_end) / 2, bracket_offset_y - bracket_height - text_offset,
                    'Responder', ha='center', va='top',
                    color=response_colors['Responder'], fontsize=10, fontweight='bold',
                    transform=ax.get_xaxis_transform())

        if nonresponder_indices:
            x_start = nonresponder_indices[0] - 0.5 + bracket_offset + bracket_gap / 2
            x_end = nonresponder_indices[-1] + 0.5 + bracket_offset
            ax.plot([x_start, x_start, x_end, x_end],
                    [bracket_offset_y, bracket_offset_y - bracket_height,
                     bracket_offset_y - bracket_height, bracket_offset_y],
                    color='black', lw=1.5, clip_on=False,
                    transform=ax.get_xaxis_transform())
            ax.text((x_start + x_end) / 2, bracket_offset_y - bracket_height - text_offset,
                    'Non-responder', ha='center', va='top',
                    color=response_colors['Non-responder'], fontsize=10, fontweight='bold',
                    transform=ax.get_xaxis_transform())
    else:
        ax.tick_params(axis='x', labelsize=xtick_fontsize)

    # Add p-values on the right side
    x_pos_pval = normalized.shape[1] + 1.5
    ax.text(x_pos_pval, -0.5, r'$p$-value', fontsize=9,
            va='center', ha='center', clip_on=False)
    for i, p_val in enumerate(p_values_clustered):
        if pd.notna(p_val):
            color = 'red' if p_val < 0.05 else 'black'
            ax.text(x_pos_pval, i + 0.5, f"{p_val:.3f}",
                    color=color, fontsize=8, va='center', ha='center')

    # Add annotation bars above heatmap
    if any(x is not None for x in [patient_sex, patient_age, patient_ctnm, patient_location]):
        from matplotlib.patches import Rectangle, Patch
        from matplotlib.lines import Line2D

        column_names = normalized.columns
        n_rows = len(normalized)

        sex_bar_y = -2.0
        ctnm_bar_y = -1.4
        location_bar_y = -0.8
        bar_height = 0.5

        for col_idx, col_name in enumerate(column_names):
            patient_id = col_name.split('_', 1)[1] if '_' in col_name else col_name

            if patient_sex is not None and patient_id in patient_sex:
                sex = patient_sex[patient_id]
                rect = Rectangle((col_idx, sex_bar_y), 1, bar_height,
                                  facecolor=sex_colors.get(sex, 'gray'),
                                  edgecolor='white', linewidth=0.5, clip_on=False)
                ax.add_patch(rect)

            if patient_ctnm is not None and patient_id in patient_ctnm:
                ctnm = patient_ctnm[patient_id]
                rect = Rectangle((col_idx, ctnm_bar_y), 1, bar_height,
                                  facecolor=ctnm_colors.get(ctnm, 'gray'),
                                  edgecolor='white', linewidth=0.5, clip_on=False)
                ax.add_patch(rect)

            if patient_location is not None and patient_id in patient_location:
                location = patient_location[patient_id]
                rect = Rectangle((col_idx, location_bar_y), 1, bar_height,
                                  facecolor=location_colors.get(location, 'gray'),
                                  edgecolor='white', linewidth=0.5, clip_on=False)
                ax.add_patch(rect)

        # Row labels for annotation bars
        if patient_sex is not None:
            ax.text(-0.5, sex_bar_y + bar_height / 2, 'Sex',
                    ha='right', va='center', fontsize=9, fontweight='bold', clip_on=False)
        if patient_ctnm is not None:
            ax.text(-0.5, ctnm_bar_y + bar_height / 2, 'Stage',
                    ha='right', va='center', fontsize=9, fontweight='bold', clip_on=False)
        if patient_location is not None:
            ax.text(-0.5, location_bar_y + bar_height / 2, 'Location',
                    ha='right', va='center', fontsize=9, fontweight='bold', clip_on=False)

        ax.set_ylim(n_rows + 0.5, -2.9)

        # Legend
        from matplotlib.lines import Line2D as _Line2D

        def _spacer():
            return _Line2D([0], [0], color='none', label=' ')

        def _header(label):
            return _Line2D([0], [0], color='none', label=label)

        legend_elements = []
        if patient_sex is not None:
            legend_elements.extend([
                _header('Sex:'),
                Patch(facecolor=sex_colors['M'], label='  Male'),
                Patch(facecolor=sex_colors['F'], label='  Female'),
            ])
        if patient_ctnm is not None:
            if legend_elements:
                legend_elements.append(_spacer())
            legend_elements.extend([
                _header('Stage:'),
                Patch(facecolor=ctnm_colors['2'], label='  2'),
                Patch(facecolor=ctnm_colors['3'], label='  3'),
            ])
        if patient_location is not None:
            if legend_elements:
                legend_elements.append(_spacer())
            legend_elements.extend([
                _header('Location:'),
                Patch(facecolor=location_colors['Lower rectum'], label='  Lower rectum'),
                Patch(facecolor=location_colors['Mid rectum'], label='  Mid rectum'),
                Patch(facecolor=location_colors['Upper rectum'], label='  Upper rectum'),
            ])

        if legend_elements and display_legends:
            legend = ax.legend(handles=legend_elements, loc='upper left',
                               bbox_to_anchor=(1.15, 1), frameon=False, fontsize=8,
                               handlelength=1)
            for text in legend.get_texts():
                if text.get_text() in ['Sex:', 'Stage:', 'Location:']:
                    text.set_fontweight('bold')
