import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Journal styling setup
sns.set_theme(style="ticks")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']

# 1. Load background full proteomics data and extract metrics
full_df = pd.read_excel('1.5倍-Diff_ResultSummary.xlsx')  # Use read_excel for .xlsx files
full_df.columns = full_df.columns.str.strip()  # Remove leading/trailing whitespace from column names

ratio_col = 'NPC43-Re/NPC43-Ctl'
q_col = 'Qvalue'

full_df[ratio_col] = pd.to_numeric(full_df[ratio_col], errors='coerce')
full_df[q_col] = pd.to_numeric(full_df[q_col], errors='coerce')
full_df = full_df.dropna(subset=[ratio_col, q_col])
full_df = full_df[(full_df[q_col] > 0) & (full_df[ratio_col] > 0)]

full_df['log2FC'] = np.log2(full_df[ratio_col])
full_df['negLog10Q'] = -np.log10(full_df[q_col])

# 2. Curate Core 11 Biphasic Signature
core_genes = ['FN1', 'S100P', 'EPHA2', 'ITGA3', 'THBS1', 'S100A10', 'STAT3', 'STMN1', 'STAT6', 'STAT1', 'PFN2']
barplot_data = full_df[full_df['Gene'].str.upper().isin([g.upper() for g in core_genes])].copy()
barplot_data = barplot_data.drop_duplicates(subset=['Gene']).sort_values(by='log2FC', ascending=True)

# -------------------------------------------------------------------------
# Plot 1: Bidirectional Diverging Barplot
# -------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 5), dpi=300)
colors = ['#2166AC' if x < 0 else '#B2182B' for x in barplot_data['log2FC']]

bars = ax.barh(barplot_data['Gene'], barplot_data['log2FC'], color=colors, edgecolor='black', linewidth=0.5, height=0.55)
ax.axvline(x=0, color='black', linewidth=1)
ax.axvline(x=np.log2(1.5), color='gray', linewidth=0.6, linestyle='--')
ax.axvline(x=-np.log2(1.5), color='gray', linewidth=0.6, linestyle='--')

ax.set_xlabel('log2(Fold Change) [NPC43-Re / NPC43-Ctl]', fontsize=10, fontweight='bold', labelpad=10)
ax.set_title('Biphasic Microenvironment Remodeling Hubs', fontsize=11, fontweight='bold', pad=15)
sns.despine(top=True, right=True)
ax.grid(axis='x', linestyle=':', alpha=0.5)

for bar in bars:
    width = bar.get_width()
    if width >= 0:
        ax.text(width + 0.08, bar.get_y() + bar.get_height()/2, f'+{width:.2f}', 
                va='center', ha='left', fontsize=8, fontweight='bold', color='#B2182B')
    else:
        ax.text(width - 0.08, bar.get_y() + bar.get_height()/2, f'{width:.2f}', 
                va='center', ha='right', fontsize=8, fontweight='bold', color='#2166AC')

ax.set_xlim(barplot_data['log2FC'].min() - 0.6, barplot_data['log2FC'].max() + 0.6)
plt.tight_layout()
plt.savefig('Bidirectional_Diverging_Barplot.pdf', bbox_inches='tight')
plt.savefig('Bidirectional_Diverging_Barplot.png', bbox_inches='tight', dpi=300)  # 可选：同时保留栅格版
plt.close()

# -------------------------------------------------------------------------
# Plot 2: Targeted Volcano Plot
# -------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 7), dpi=300)
ax.scatter(full_df['log2FC'], full_df['negLog10Q'], color='#EAEAEA', alpha=0.6, s=10, label='All Detected Proteins', zorder=1)

up_targets = barplot_data[barplot_data['log2FC'] > 0]
down_targets = barplot_data[barplot_data['log2FC'] < 0]

ax.scatter(up_targets['log2FC'], up_targets['negLog10Q'], color='#B2182B', s=60, label='High_CP Track (Curated Up)', edgecolors='black', linewidths=0.6, zorder=5)
ax.scatter(down_targets['log2FC'], down_targets['negLog10Q'], color='#2166AC', s=60, label='Low_CP Track (Curated Down)', edgecolors='black', linewidths=0.6, zorder=5)

ax.axhline(y=-np.log10(0.05), color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
ax.axvline(x=np.log2(1.5), color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
ax.axvline(x=-np.log2(1.5), color='gray', linestyle='--', linewidth=0.8, alpha=0.6)

# Non-overlapping static alignment offsets mapping
custom_offsets = {
    'S100P': (0.3, -1.5), 'FN1': (0.3, 0.8), 'EPHA2': (-0.7, 1.2), 'ITGA3': (0.3, -0.5),
    'THBS1': (-0.6, -1.5), 'S100A10': (0.2, 0.8), 'STMN1': (-0.7, 1.2), 'STAT1': (0.2, 1.2),
    'STAT3': (-0.6, -1.5), 'STAT6': (0.2, -1.2), 'PFN2': (0.2, 0.5)
}

for i, row in barplot_data.iterrows():
    gene = row['Gene']
    x, y = row['log2FC'], row['negLog10Q']
    ox, oy = custom_offsets.get(gene, (0.2, 0.2))
    ax.annotate(gene, (x, y), xytext=(x + ox, y + oy), fontsize=8.5, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75),
                arrowprops=dict(arrowstyle="-", color='black', lw=0.5, alpha=0.5), zorder=6)

ax.set_xlabel('log2(Fold Change)', fontsize=10, fontweight='bold')
ax.set_ylabel('-log10(Q-value)', fontsize=10, fontweight='bold')
ax.set_title('Targeted Volcano Plot of Lytic Remodeling Cascade', fontsize=11, fontweight='bold', pad=15)
ax.legend(loc='upper right', frameon=True, fontsize=8.5, facecolor='white', edgecolor='none')
sns.despine()
ax.grid(linestyle=':', alpha=0.3)

plt.tight_layout()
plt.savefig('Targeted_Volcano_Plot.pdf', bbox_inches='tight')
plt.savefig('Targeted_Volcano_Plot.png', bbox_inches='tight', dpi=300)
plt.close()