import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

table = pd.read_csv('sirus-trial1-ldlc-dynamics-median-IQR.csv')

fig = plt.figure(figsize=(8, 6), layout='tight')
ax = fig.add_subplot()
ax = sns.lineplot(table, ax=ax, x='Time', y='Value', hue='Arm', style='Arm', hue_order=['Placebo', 'Inclisiran'],
                  markers=True, markersize=15, lw=4, errorbar=None, estimator=None, n_boot=0)
ax.set_ylim(0, None)
sns.move_legend(ax, 'lower left', bbox_to_anchor=(0, 0), frameon=True)
plt.savefig('figure1.pdf')
