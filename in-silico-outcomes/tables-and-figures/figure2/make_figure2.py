import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

tableA = pd.read_csv('raw-data-figure-2/mace.csv')
tableB = pd.read_csv('raw-data-figure-2/cv-death.csv')
tableC = pd.read_csv('raw-data-figure-2/mi.csv')
tableD = pd.read_csv('raw-data-figure-2/isc-stroke.csv')
tableE = pd.read_csv('raw-data-figure-2/male.csv')
table = pd.concat([tableA, tableB, tableC, tableD, tableE])
table['Cumulative incidence %'] = 100 * table['KM_estimate']

fig = plt.figure(figsize=(9, 8), layout='tight')
ax = fig.add_subplot()
ax = sns.relplot(table, x='timeline', y='Cumulative incidence %', hue='Arm', hue_order=['Placebo', 'Inclisiran'],
                 kind='line', lw=1.2, col='Event', col_order=['3P-MACE', 'CvDeath', 'MI', 'IscStroke', 'MALE'],
                 aspect=1.3, col_wrap=2, height=3.8, legend=True,
                 facet_kws={'legend_out': False, 'sharey': False, 'sharex': False})
sns.move_legend(ax, 'lower left', bbox_to_anchor=(60, 0), frameon=True)
plt.savefig('figure2.pdf')
