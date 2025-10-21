import numpy as np
import pandas as pd
from lifelines import CoxPHFitter


def get_bootstrap_interactions(n, nb_boot_repetitions, sample_size, variable):
    m = CoxPHFitter()
    table1 = table.loc[table[variable] == 1, ['Arm', 'timeOfMACE', 'isMACE']]
    m.fit(table1, duration_col='timeOfMACE', event_col='isMACE')
    p1 = np.array(m.params_)[0]
    table2 = table.loc[table[variable] == 0, ['Arm', 'timeOfMACE', 'isMACE']]
    m.fit(table2, duration_col='timeOfMACE', event_col='isMACE')
    p2 = np.array(m.params_)[0]
    point_estimate = p2 - p1

    all_indices = np.array(range(n))
    boot_estimates = []
    for _ in range(nb_boot_repetitions):
        indices_arm0 = np.random.choice(all_indices, size=sample_size, replace=True)
        boot_table_arm0 = table_arm0.iloc[indices_arm0, :]

        indices_arm1 = np.random.choice(all_indices, size=sample_size, replace=True)
        boot_table_arm1 = table_arm1.iloc[indices_arm1, :]

        boot_table = pd.concat([boot_table_arm0, boot_table_arm1])

        m.fit(boot_table.loc[boot_table[variable] == 1, ['Arm', 'timeOfMACE', 'isMACE']],
            duration_col='timeOfMACE', event_col='isMACE')
        p1 = np.array(m.params_)[0]

        m.fit(boot_table.loc[boot_table[variable] == 0, ['Arm', 'timeOfMACE', 'isMACE']],
            duration_col='timeOfMACE', event_col='isMACE')
        p2 = np.array(m.params_)[0]

        boot_estimates.append(p2 - p1)

    results = float(point_estimate), float(np.percentile(boot_estimates, 2.5)), float(np.percentile(boot_estimates, 97.5))
    return variable + ',' + ','.join(map(lambda x: str(round(x, 3)), results))


np.random.seed(42)

table = pd.read_csv('simu-data.csv')

nb_repetitions = 100
size_FOURIER = 13780

table['Arm'] = np.array(table['Arm'] == 'RWLLT-ICL', dtype=np.int64)
table['isAgeAbove65'] = np.array(table['ageInit'] > 65, dtype=np.int64)
table['isAgeAbove70'] = np.array(table['ageInit'] > 70, dtype=np.int64)
table['isEgfrBelow60'] = np.array(table['eGFR'] < 60, dtype=np.int64)
table['isMoreThanOneBed'] = np.array(table['nbPlaques'] > 1, dtype=np.int64)
table['isMoreThanOneRF'] = np.array(table['nbRiskFactors'] > 1, dtype=np.int64)

table_arm0 = table[table['Arm'] == 0]
table_arm1 = table[table['Arm'] == 1]

print(get_bootstrap_interactions(204691, nb_repetitions, size_FOURIER, 'isAgeAbove65'))
print(get_bootstrap_interactions(204691, nb_repetitions, size_FOURIER, 'isAgeAbove70'))
print(get_bootstrap_interactions(204691, nb_repetitions, size_FOURIER, 'isMaleSex'))
print(get_bootstrap_interactions(204691, nb_repetitions, size_FOURIER, 'isDiabetes'))
print(get_bootstrap_interactions(204691, nb_repetitions, size_FOURIER, 'isSmoking'))
print(get_bootstrap_interactions(204691, nb_repetitions, size_FOURIER, 'isUncontrolledBloodPressure'))
print(get_bootstrap_interactions(204691, nb_repetitions, size_FOURIER, 'isEgfrBelow60'))
print(get_bootstrap_interactions(204691, nb_repetitions, size_FOURIER, 'isPriorChd'))
print(get_bootstrap_interactions(204691, nb_repetitions, size_FOURIER, 'isPriorCbvd'))
print(get_bootstrap_interactions(204691, nb_repetitions, size_FOURIER, 'isPriorPad'))
print(get_bootstrap_interactions(204691, nb_repetitions, size_FOURIER, 'isMoreThanOneBed'))
print(get_bootstrap_interactions(204691, nb_repetitions, size_FOURIER, 'isMoreThanOneRF'))


table = pd.read_csv('simu-data-bis.csv')

table['Arm'] = np.array(table['Arm'] == 'RWLLT-ICL', dtype=np.int64)
table['isLDLcAboveMedian'] = np.array(table['baselineLDLcQuartile'] >= 2, dtype=np.int64)

table_arm0 = table[table['Arm'] == 0]
table_arm1 = table[table['Arm'] == 1]

print(get_bootstrap_interactions(204691, nb_repetitions, size_FOURIER, 'isLDLcAboveMedian'))
