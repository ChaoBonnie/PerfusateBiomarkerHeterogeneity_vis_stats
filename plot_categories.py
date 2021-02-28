import pandas as pd
import seaborn as sns; sns.set(style='white', context='paper')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from statannot import add_stat_annotation


### Graph Swarm Plots of Biomarker Expressions in Each Perfusate Sample ###


ys = ['IL-1β (pg/mL)', 'IL-6 (pg/mL)', 'IL-8 (pg/mL)', 'IL-10 (pg/mL)']

location_map = {'RU': 'Right Upper', 'RL': 'Right Lower', 'LU': 'Left Upper', 'LL': 'Left Lower'}

# Annotate significance levels from an imported table #

def annotate_anova(ax, data, y, anova_path, anova_sheet):
    df = pd.read_excel(anova_path, sheet_name=anova_sheet, index_col=0)
    df = df[y.split(' ')[0]]
    pvalues = []
    box_pairs = []
    for x in df.index:
        p = df[x]
        if p < 0.05:
            pvalues.append(p)
            box_pairs.append(((x, 'RU'), (x, 'LL')))
    add_stat_annotation(ax, data=data, x='EVLP ID', y=y, hue='Location',
                        box_pairs=box_pairs, pvalues=pvalues, perform_stat_test=False,
                        loc='outside', verbose=0)

# Plot horizontal lines on existing graphs #

def horizontal_lines(data, ax, x, y, **kwargs):
    data = data.copy()
    xticks = ax.get_xticklabels()
    xticks = {tick.get_text(): i for i, tick in enumerate(xticks)}
    data['xmin'] = data[x].apply(lambda xval: xticks[str(xval)] - 0.4)
    data['xmax'] = data[x].apply(lambda xval: xticks[str(xval)] + 0.4)
    ax.hlines(y=data[y], xmin=data['xmin'], xmax=data['xmax'], **kwargs)
    return ax

# Make swarm plots with horizontal lines showing mean vluaes and significance annotations #

def make_graph(data_path, data_sheet, anova_path=None, anova_sheet=None, show_fig=False):
    raw_df = pd.read_excel(data_path, sheet_name=data_sheet)
    raw_df = raw_df.drop(columns=['Kit ID'])
    df = pd.DataFrame(raw_df)

    fig = plt.figure(figsize=(15, 3))
    gs = GridSpec(1, 4)
    axs = []
    mypalette = ["royalblue", "orange", "#95a5a6", "#e74c3c", "#34495e", "#3498db", "#2ecc71"]
    sns.set_palette(mypalette)

    for i in range(len(ys)):
        axs.append(fig.add_subplot(gs[0, i]))

    for ax, y in zip(axs, ys):
        swarm = sns.swarmplot(x='EVLP ID', y=y, hue='Location', data=df, ax=ax, size=3.1)
        swarm.axhline(df[y].mean(), color='green', alpha=0.5, linewidth=2, linestyle='dashed')
        df_subj_mean = df.groupby('EVLP ID', as_index=False)[y].mean()
        horizontal_lines(data=df_subj_mean, ax=swarm, x='EVLP ID', y=y, color='black', linewidth=2)
        ax.xaxis.labelpad = 15
        ax.yaxis.labelpad = 5
        ax.set_xlabel(ax.get_xlabel(), fontsize=12)
        ax.set_ylabel(ax.get_ylabel(), fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=10.5)
        ax.legend().remove()
        if anova_path is not None:
            annotate_anova(ax, df, y, anova_path, anova_sheet)

    handles, labels = axs[0].get_legend_handles_labels()
    labels = [location_map[l] for l in labels]
    fig.legend(handles, labels, loc=(0.805, 0.68))

    fig.tight_layout()

    if show_fig:
        fig.show()
    else:
        fig.savefig(f'Cat Plots in a row without ANOVA.png', dpi=200)

make_graph(r'C:\Users\chaob\Documents\Isolated Lung Perfusate Sampling Summary.xlsx', '4 Locations')
make_graph(r'C:\Users\chaob\Documents\Isolated Lung Perfusate Sampling Summary.xlsx', '4 Locations',
           'ANOVA Result Table.xlsx', '4 Locations')
make_graph(r'C:\Users\chaob\Documents\Isolated Lung Perfusate Sampling Summary.xlsx', '2 Locations',
           'ANOVA Result Table.xlsx', '2 Locations')


### Coefficient of Variance Bar Plots from Cytokine and Biochemistry Data ###


df = pd.read_excel(r'C:\Users\chaob\Documents\Perfusate Heterogeneity %CV Plotting.xlsx', 'CV Plotting')
print(df)

y_cytokine = ['IL-6', 'IL-8', 'sTNF-R1', 'sTREM-1',
      'ET-1', 'GM-CSF', 'IL-10', 'IL-1β']

fig_cytokine = plt.figure(figsize=(32, 16))
gs = GridSpec(2, 4)
axs = []
for i in range(len(y_cytokine)):
    row = i // gs.ncols
    col = i % gs.ncols
    axs.append(fig_cytokine.add_subplot(gs[row, col]))

for ax, y in zip(axs, y_cytokine):
    graph_cytokine = sns.barplot(x='EVLP ID', y=y, data=df, ax=ax, palette='Set2')
    graph_cytokine.axhline(df[y].mean())
    # sns.lineplot(x='EVLP ID', y=df[y].mean(), data=df, ax=ax)

fig_cytokine.tight_layout()
fig_cytokine.savefig('Cytokine_%CV.png', dpi=200)


y_biochem = ['pH', 'dPO2', 'HCO3', 'BE', 'Na',
      'K', 'Ca', 'Cl', 'Glu', 'Lac']

fig_biochem = plt.figure(figsize=(32, 20))
gs = GridSpec(2, 5)
axs = []
for i in range(len(y_biochem)):
    row = i // gs.ncols
    col = i % gs.ncols
    axs.append(fig_biochem.add_subplot(gs[row, col]))

for ax, y in zip(axs, y_biochem):
    graph_biochem = sns.barplot(x='EVLP ID', y=y, data=df, ax=ax, palette='Set2')
    graph_biochem.axhline(df[y].mean())
    # sns.lineplot(x='EVLP ID', y=df[y].mean(), data=df, ax=ax)

fig_biochem.tight_layout()
fig_biochem.savefig('Biochem_%CV.png', dpi=200)


### Correlations of dPO2 with Other Important EVLP Parameters ###


df = pd.read_excel(r'C:\Users\chaob\Documents\Perfusate Heterogeneity %CV Plotting.xlsx', 'dPO2 Correlations')
# markers = ['o', '^', 's', 'P', '*', 'D', 'v', 'X', '>']

fig = plt.figure(figsize=(12, 6))
sns.swarmplot(x='EVLP ID', y='Correlation with dPO2', data=df, hue='Biomarker', size=10)
plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
fig.tight_layout()
fig.savefig('dPO2_Correlation_Scatterplot.png', dpi=200)


### Intra-class Correlations ###


import pingouin as pg

sheet_names = ['#389', '#393', '#398', '#406', '#422', '#426', '#472', '#474']

for s in sheet_names:
    df = pd.read_excel(r'C:\Users\chaob\Documents\Perfusate Protein Data by Donor.xlsx', s)
    icc = pg.intraclass_corr(data=df, targets='Target Cytokine', raters='Location',
                             ratings='Expression').round(3)
    print(s, '\n', icc)