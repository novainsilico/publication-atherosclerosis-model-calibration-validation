"""
Processing raw calibration-simulation output data
(big CSV with all lipoproteins, CRP and PCKS9 measured at all timepoints)
(The script can be adapted to only plot LDL-C median dynamics using "scalars-trimmed-for-fourier-ldlc")
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def set_calib_globals():
    global base_folder, simu_path, ordered_arms, arm_names

    base_folder = 'vpop-lipo-data/'
    simu_path = 'scalars-FOURIER-calib-all-batches-merged.csv'

    # always put the placebo arm first
    ordered_arms = ['Placebo', 'Evolocumab']
    arm_names = {'FOURIER-placebo': 'Placebo', 'FOURIER-evolocumab': 'Evolocumab', 'ORION10-inclisiran-300mg': 'Inclisiran'}


def gen_lipo_median_dynamics_data(has_ref, variable_name, descriptor_name, csv_path):
    simu_table = pd.read_csv(os.path.join(base_folder, simu_path))
    if has_ref:
        ref_path = 'literature-lipo-data/FOURIER-entities-dynamics.csv'
        ref_table = pd.read_csv(os.path.join(base_folder, ref_path))
        ref_table = ref_table[ref_table['Condition'].isna()]

    simu_times = [0, 4, 12, 24, 48, 72, 96, 120, 144]
    descriptors = [f'{descriptor_name}W{timepoint}' for timepoint in simu_times]

    long_table = pd.DataFrame()
    for arm, arm_name in arm_names.items():
        if arm_name == 'Inclisiran':
            continue
        table = {'Value': [], 'Time': []}
        for timepoint, descriptor in zip(simu_times, descriptors):
            table['Value'].append(np.median(simu_table[f'{arm}_{descriptor}']))
            table['Time'].append(timepoint)
        table = pd.DataFrame(table)
        table['Variable'] = variable_name
        table['Arm'] = arm_name
        table['Origin'] = 'Simulated data (Vpop)'
        long_table = pd.concat([long_table, table], ignore_index=True)

        if has_ref:
            sub_table = ref_table[ref_table['Variable'].apply(lambda x: x.startswith(arm))]
            table = pd.DataFrame({'Value': sub_table['Value'], 'Time': sub_table['TimeWeeks']})
            table['Variable'] = variable_name
            table['Arm'] = arm_name
            table['Origin'] = 'Observed data (FOURIER)'
            long_table = pd.concat([long_table, table], ignore_index=True)
    long_table.to_csv(csv_path, index=False)


def plot_median_dynamics():
    cats = ['LDL-C (mg/dL)', 'Lp(a) (nmol/L)', 'PCSK9 (ng/mL)']
    table1 = pd.read_csv('vpop-lipo-data/ldlc-median-dynamics.csv')
    table1['Cat'] = cats[0]
    table2 = pd.read_csv('vpop-lipo-data/lpa-median-dynamics.csv')
    table2['Cat'] = cats[1]
    table3 = pd.read_csv('vpop-lipo-data/pcsk9-median-dynamics.csv')
    table3['Cat'] = cats[2]
    long_table = pd.concat([table1, table2, table3], ignore_index=True)
    long_table['\nOrigin'] = long_table['Origin']
    g = sns.relplot(long_table, x='Time', y='Value', hue='Arm', hue_order=ordered_arms,
                    kind='line', style='\nOrigin',
                    markers=True, markersize=8, lw=2, errorbar=None, estimator=None, n_boot=0,
                    col='Cat', col_order=cats,
                    aspect=1.6, col_wrap=2, height=3.6, legend=True,
                    facet_kws={'legend_out': False, 'sharey': False, 'sharex': True})

    g.set_titles('A').tight_layout(w_pad=3).set_ylabels('').set_xlabels('Time (days)')
    for i, ax in enumerate(g.axes):
        ax.set_title(f'({chr(65+i)}) Median {cats[i].split()[0]} over time')
        ax.set_ylabel(cats[i])
    g.axes[0].set_ylim(-1, 130)
    g.axes[1].set_ylim(-1, 60)
    g.axes[2].set_ylim(-5, 500)
    sns.move_legend(g, 'upper left', bbox_to_anchor=(.66, .40), frameon=True)
    plt.savefig('vpop-lipo-data/FOURIER-lipo-median-dynamics.png')


set_calib_globals()
gen_lipo_median_dynamics_data(True, 'LDL-C (mg/dL)', 'valueLDLcUncorrectedMass', 'vpop-lipo-data/FOURIER-ldlc-median-dynamics.csv')
gen_lipo_median_dynamics_data(False, 'Lp(a) (nmol/L)', 'valueLpaSubstance', 'vpop-lipo-data/FOURIER-lpa-median-dynamics.csv')
gen_lipo_median_dynamics_data(False, 'PCSK9 (ng/mL)', 'valueFreePCSK9', 'vpop-lipo-data/FOURIER-pcsk9-median-dynamics.csv')
plot_median_dynamics()
