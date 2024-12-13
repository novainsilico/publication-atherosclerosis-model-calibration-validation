import os
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter


base_folder = 'vpop-event-data/'
input_csv = 'scalars-trimmed-for-events.csv'

arms = {'FOURIER-placebo': 0,  'FOURIER-evolocumab': 1}

ldlc = 'FOURIER-placebo_baselineLDLcUncorrectedMass'
lpa = 'FOURIER-placebo_baselineLpaSubstance'
his = 'FOURIER-placebo_isHighIntStatinMono'
hiseze = 'FOURIER-placebo_isHighIntStatinEze'

table = pd.read_csv(os.path.join(base_folder, input_csv))


def print_hr(big_table, t, e):
    small_table = pd.DataFrame()
    for arm, arm_value in arms.items():
        arm_table = pd.DataFrame({'time': big_table[f'{arm}_{t}'], 'status': big_table[f'{arm}_{e}']})
        arm_table['Arm'] = arm_value
        small_table = pd.concat([small_table, arm_table], ignore_index=True)

    mdl = CoxPHFitter()
    mdl.fit(small_table, duration_col='time', event_col='status')
    print(e, '{:.2f}'.format(np.array(mdl.hazard_ratios_)[0]))


def overall():
    print_hr(table, 'timeOfMACE', 'isMACE')
    print_hr(table, 'timeOfCvDeath', 'isCvDeath')
    print_hr(table, 'obsTimeOfMI', 'isFatalCHD')
    print_hr(table, 'obsTimeOfIscStroke', 'isFatalIscStroke')
    print_hr(table, 'timeOfCvDeath', 'isOtherCvDeath')
    print_hr(table, 'obsTimeOfMI', 'isMI')
    print_hr(table, 'obsTimeOfIscStroke', 'isIscStroke')
    print_hr(table, 'obsTimeOfMajorAdverseLimbEvent2_5Y', 'isMajorAdverseLimbEvent2_5Y')


def key_subgroups():
    print_hr(table[table['ageInit'] < 65], 'timeOfMACE', 'isMACE')
    print_hr(table[table['ageInit'] >= 65], 'timeOfMACE', 'isMACE')
    print_hr(table[table['isMaleSex'] == 0], 'timeOfMACE', 'isMACE')
    print_hr(table[table['isMaleSex'] == 1], 'timeOfMACE', 'isMACE')
    print_hr(table[(table['isPriorChd'] == 1) & (table['bedInvolvement'] == 1)], 'timeOfMACE', 'isMACE')
    print_hr(table[(table['isPriorCbvd'] == 1) & (table['bedInvolvement'] == 1)], 'timeOfMACE', 'isMACE')
    print_hr(table[(table['isPriorPad'] == 1) & (table['bedInvolvement'] == 1)], 'timeOfMACE', 'isMACE')
    print_hr(table[table['bedInvolvement'] > 1.5], 'timeOfMACE', 'isMACE')
    print_hr(table[table[ldlc] < 80], 'timeOfMACE', 'isMACE')
    print_hr(table[(table[ldlc] >= 80) & (table[ldlc] < 92)], 'timeOfMACE', 'isMACE')
    print_hr(table[(table[ldlc] >= 92) & (table[ldlc] < 109)], 'timeOfMACE', 'isMACE')
    print_hr(table[table[ldlc] >= 109], 'timeOfMACE', 'isMACE')
    print_hr(table[(table[his] == 1) | (table[hiseze] == 1)], 'timeOfMACE', 'isMACE')
    print_hr(table[(table[his] == 0) & (table[hiseze] == 0)], 'timeOfMACE', 'isMACE')
    print_hr(table[table[hiseze] == 1], 'timeOfMACE', 'isMACE')
    print_hr(table[table[hiseze] == 0], 'timeOfMACE', 'isMACE')

    print_hr(table[table['ageInit'] <= 56], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[(table['ageInit'] > 56) & (table['ageInit'] <= 63)], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[(table['ageInit'] > 63) & (table['ageInit'] <= 69)], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[table['ageInit'] > 69], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[table['isMaleSex'] == 1], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[table['isMaleSex'] == 0], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[table['isDiabetes'] == 1], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[table['isDiabetes'] == 0], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[table['eGFR'] < 60], 'timeOfCvDeath2_5Y', 'isCvDeath2_5Y')
    print_hr(table[(table['eGFR'] >= 60) & (table['eGFR'] <= 90)], 'timeOfCvDeath2_5Y', 'isCvDeath2_5Y')
    print_hr(table[table['eGFR'] >= 90], 'timeOfCvDeath2_5Y', 'isCvDeath2_5Y')

    print_hr(table[table['ageInit'] <= 56], 'obsTimeOfMI3Y', 'isFatalCHD3Y')
    print_hr(table[(table['ageInit'] > 56) & (table['ageInit'] <= 63)], 'obsTimeOfMI3Y', 'isFatalCHD3Y')
    print_hr(table[(table['ageInit'] > 63) & (table['ageInit'] <= 69)], 'obsTimeOfMI3Y', 'isFatalCHD3Y')
    print_hr(table[table['ageInit'] > 69], 'obsTimeOfMI3Y', 'isFatalCHD3Y')
    print_hr(table[table['isMaleSex'] == 1], 'obsTimeOfMI3Y', 'isFatalCHD3Y')
    print_hr(table[table['isMaleSex'] == 0], 'obsTimeOfMI3Y', 'isFatalCHD3Y')
    print_hr(table[table['isDiabetes'] == 1], 'obsTimeOfMI3Y', 'isFatalCHD3Y')
    print_hr(table[table['isDiabetes'] == 0], 'obsTimeOfMI3Y', 'isFatalCHD3Y')
    print_hr(table[table['eGFR'] < 60], 'obsTimeOfMI2_5Y', 'isFatalCHD2_5Y')
    print_hr(table[(table['eGFR'] >= 60) & (table['eGFR'] <= 90)], 'obsTimeOfMI2_5Y', 'isFatalCHD2_5Y')
    print_hr(table[table['eGFR'] >= 90], 'obsTimeOfMI2_5Y', 'isFatalCHD2_5Y')

    print_hr(table[table['ageInit'] <= 56], 'obsTimeOfIscStroke3Y', 'isFatalIscStroke3Y')
    print_hr(table[(table['ageInit'] > 56) & (table['ageInit'] <= 63)], 'obsTimeOfIscStroke3Y', 'isFatalIscStroke3Y')
    print_hr(table[(table['ageInit'] > 63) & (table['ageInit'] <= 69)], 'obsTimeOfIscStroke3Y', 'isFatalIscStroke3Y')
    print_hr(table[table['ageInit'] > 69], 'obsTimeOfIscStroke3Y', 'isFatalIscStroke3Y')
    print_hr(table[table['isMaleSex'] == 1], 'obsTimeOfIscStroke3Y', 'isFatalIscStroke3Y')
    print_hr(table[table['isMaleSex'] == 0], 'obsTimeOfIscStroke3Y', 'isFatalIscStroke3Y')
    print_hr(table[table['isDiabetes'] == 1], 'obsTimeOfIscStroke3Y', 'isFatalIscStroke3Y')
    print_hr(table[table['isDiabetes'] == 0], 'obsTimeOfIscStroke3Y', 'isFatalIscStroke3Y')
    print_hr(table[table['eGFR'] < 60], 'obsTimeOfIscStroke2_5Y', 'isFatalIscStroke2_5Y')
    print_hr(table[(table['eGFR'] >= 60) & (table['eGFR'] <= 90)], 'obsTimeOfIscStroke2_5Y', 'isFatalIscStroke2_5Y')
    print_hr(table[table['eGFR'] >= 90], 'obsTimeOfIscStroke2_5Y', 'isFatalIscStroke2_5Y')


def by_age():
    print_hr(table[table['ageInit'] <= 56], 'timeOfMACE3Y', 'isMACE3Y')
    print_hr(table[(table['ageInit'] > 56) & (table['ageInit'] <= 63)], 'timeOfMACE3Y', 'isMACE3Y')
    print_hr(table[(table['ageInit'] > 63) & (table['ageInit'] <= 69)], 'timeOfMACE3Y', 'isMACE3Y')
    print_hr(table[table['ageInit'] > 69], 'timeOfMACE3Y', 'isMACE3Y')
    print_hr(table[table['ageInit'] <= 56], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[(table['ageInit'] > 56) & (table['ageInit'] <= 63)], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[(table['ageInit'] > 63) & (table['ageInit'] <= 69)], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[table['ageInit'] > 69], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[table['ageInit'] <= 56], 'obsTimeOfMI3Y', 'isMI3Y')
    print_hr(table[(table['ageInit'] > 56) & (table['ageInit'] <= 63)], 'obsTimeOfMI3Y', 'isMI3Y')
    print_hr(table[(table['ageInit'] > 63) & (table['ageInit'] <= 69)], 'obsTimeOfMI3Y', 'isMI3Y')
    print_hr(table[table['ageInit'] > 69], 'obsTimeOfMI3Y', 'isMI3Y')
    print_hr(table[table['ageInit'] <= 56], 'obsTimeOfIscStroke3Y', 'isIscStroke3Y')
    print_hr(table[(table['ageInit'] > 56) & (table['ageInit'] <= 63)], 'obsTimeOfIscStroke3Y', 'isIscStroke3Y')
    print_hr(table[(table['ageInit'] > 63) & (table['ageInit'] <= 69)], 'obsTimeOfIscStroke3Y', 'isIscStroke3Y')
    print_hr(table[table['ageInit'] > 69], 'obsTimeOfIscStroke3Y', 'isIscStroke3Y')


def by_sex():
    print_hr(table[table['isMaleSex'] == 1], 'timeOfMACE3Y', 'isMACE3Y')
    print_hr(table[table['isMaleSex'] == 0], 'timeOfMACE3Y', 'isMACE3Y')
    print_hr(table[table['isMaleSex'] == 1], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[table['isMaleSex'] == 0], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[table['isMaleSex'] == 1], 'obsTimeOfMI3Y', 'isMI3Y')
    print_hr(table[table['isMaleSex'] == 0], 'obsTimeOfMI3Y', 'isMI3Y')
    print_hr(table[table['isMaleSex'] == 1], 'obsTimeOfIscStroke3Y', 'isIscStroke3Y')
    print_hr(table[table['isMaleSex'] == 0], 'obsTimeOfIscStroke3Y', 'isIscStroke3Y')


def by_egfr():
    print_hr(table[table['eGFR'] < 60], 'timeOfMACE2_5Y', 'isMACE2_5Y')
    print_hr(table[(table['eGFR'] >= 60) & (table['eGFR'] < 90)], 'timeOfMACE2_5Y', 'isMACE2_5Y')
    print_hr(table[table['eGFR'] >= 90], 'timeOfMACE2_5Y', 'isMACE2_5Y')
    print_hr(table[table['eGFR'] < 60], 'timeOfCvDeath2_5Y', 'isCvDeath2_5Y')
    print_hr(table[(table['eGFR'] >= 60) & (table['eGFR'] < 90)], 'timeOfCvDeath2_5Y', 'isCvDeath2_5Y')
    print_hr(table[table['eGFR'] >= 90], 'timeOfCvDeath2_5Y', 'isCvDeath2_5Y')
    print_hr(table[table['eGFR'] < 60], 'obsTimeOfMI2_5Y', 'isMI2_5Y')
    print_hr(table[(table['eGFR'] >= 60) & (table['eGFR'] < 90)], 'obsTimeOfMI2_5Y', 'isMI2_5Y')
    print_hr(table[table['eGFR'] >= 90], 'obsTimeOfMI2_5Y', 'isMI2_5Y')
    print_hr(table[table['eGFR'] < 60], 'obsTimeOfIscStroke2_5Y', 'isIscStroke2_5Y')
    print_hr(table[(table['eGFR'] >= 60) & (table['eGFR'] < 90)], 'obsTimeOfIscStroke2_5Y', 'isIscStroke2_5Y')
    print_hr(table[table['eGFR'] >= 90], 'obsTimeOfIscStroke2_5Y', 'isIscStroke2_5Y')


def by_diabetes():
    print_hr(table[table['isDiabetes'] == 1], 'timeOfMACE3Y', 'isMACE3Y')
    print_hr(table[table['isDiabetes'] == 0], 'timeOfMACE3Y', 'isMACE3Y')
    print_hr(table[table['isDiabetes'] == 1], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[table['isDiabetes'] == 0], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[table['isDiabetes'] == 1], 'obsTimeOfMI3Y', 'isMI3Y')
    print_hr(table[table['isDiabetes'] == 0], 'obsTimeOfMI3Y', 'isMI3Y')
    print_hr(table[table['isDiabetes'] == 1], 'obsTimeOfIscStroke3Y', 'isIscStroke3Y')
    print_hr(table[table['isDiabetes'] == 0], 'obsTimeOfIscStroke3Y', 'isIscStroke3Y')


def by_recent_remote_mi():
    print_hr(table[(table['isEventInPastYearAtStart'] == 1) & (table['isPriorChd'] == 1)], 'timeOfMACE3Y', 'isMACE3Y')
    print_hr(table[(table['isEventInPastYearAtStart'] == 0) & (table['isPriorChd'] == 1)], 'timeOfMACE3Y', 'isMACE3Y')
    print_hr(table[(table['isEventInPastYearAtStart'] == 1) & (table['isPriorChd'] == 1)], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[(table['isEventInPastYearAtStart'] == 0) & (table['isPriorChd'] == 1)], 'timeOfCvDeath3Y', 'isCvDeath3Y')
    print_hr(table[(table['isEventInPastYearAtStart'] == 1) & (table['isPriorChd'] == 1)], 'obsTimeOfMI3Y', 'isMI3Y')
    print_hr(table[(table['isEventInPastYearAtStart'] == 0) & (table['isPriorChd'] == 1)], 'obsTimeOfMI3Y', 'isMI3Y')
    print_hr(table[(table['isEventInPastYearAtStart'] == 1) & (table['isPriorChd'] == 1)], 'obsTimeOfIscStroke3Y', 'isIscStroke3Y')
    print_hr(table[(table['isEventInPastYearAtStart'] == 0) & (table['isPriorChd'] == 1)], 'obsTimeOfIscStroke3Y', 'isIscStroke3Y')


def by_prior_stroke():
    print_hr(table[table['isPriorCbvd'] == 1], 'timeOfMACE', 'isMACE')
    print_hr(table[table['isPriorCbvd'] == 0], 'timeOfMACE', 'isMACE')
    print_hr(table[table['isPriorCbvd'] == 1], 'timeOfCvDeath', 'isCvDeath')
    print_hr(table[table['isPriorCbvd'] == 0], 'timeOfCvDeath', 'isCvDeath')
    print_hr(table[table['isPriorCbvd'] == 1], 'obsTimeOfMI', 'isMI')
    print_hr(table[table['isPriorCbvd'] == 0], 'obsTimeOfMI', 'isMI')
    print_hr(table[table['isPriorCbvd'] == 1], 'obsTimeOfIscStroke', 'isIscStroke')
    print_hr(table[table['isPriorCbvd'] == 0], 'obsTimeOfIscStroke', 'isIscStroke')


def by_prior_pad():
    print_hr(table[table['isPriorPad'] == 1], 'timeOfMACE2_5Y', 'isMACE2_5Y')
    print_hr(table[(table['isPriorPad'] == 1) & (table['bedInvolvement'] == 1)], 'timeOfMACE2_5Y', 'isMACE2_5Y')
    print_hr(table[table['isPriorPad'] == 1], 'timeOfCvDeath2_5Y', 'isCvDeath2_5Y')
    print_hr(table[(table['isPriorPad'] == 1) & (table['bedInvolvement'] == 1)], 'timeOfCvDeath2_5Y', 'isCvDeath2_5Y')
    print_hr(table[table['isPriorPad'] == 1], 'obsTimeOfMI2_5Y', 'isMI2_5Y')
    print_hr(table[(table['isPriorPad'] == 1) & (table['bedInvolvement'] == 1)], 'obsTimeOfMI2_5Y', 'isMI2_5Y')
    print_hr(table[table['isPriorPad'] == 1], 'obsTimeOfIscStroke2_5Y', 'isIscStroke2_5Y')
    print_hr(table[(table['isPriorPad'] == 1) & (table['bedInvolvement'] == 1)], 'obsTimeOfIscStroke2_5Y', 'isIscStroke2_5Y')
    print_hr(table[table['isPriorPad'] == 1], 'obsTimeOfMajorAdverseLimbEvent2_5Y', 'isMajorAdverseLimbEvent2_5Y')
    print_hr(table[(table['isPriorPad'] == 1) & (table['bedInvolvement'] == 1)], 'obsTimeOfMajorAdverseLimbEvent2_5Y', 'isMajorAdverseLimbEvent2_5Y')


def by_statin_potency():
    print_hr(table[(table[his] == 1) | (table[hiseze] == 1)], 'timeOfMACE', 'isMACE')
    print_hr(table[(table[his] == 0) & (table[hiseze] == 0)], 'timeOfMACE', 'isMACE')
    print_hr(table[(table[his] == 1) | (table[hiseze] == 1)], 'timeOfCvDeath', 'isCvDeath')
    print_hr(table[(table[his] == 0) & (table[hiseze] == 0)], 'timeOfCvDeath', 'isCvDeath')
    print_hr(table[(table[his] == 1) | (table[hiseze] == 1)], 'obsTimeOfMI', 'isMI')
    print_hr(table[(table[his] == 0) & (table[hiseze] == 0)], 'obsTimeOfMI', 'isMI')
    print_hr(table[(table[his] == 1) | (table[hiseze] == 1)], 'obsTimeOfIscStroke', 'isIscStroke')
    print_hr(table[(table[his] == 0) & (table[hiseze] == 0)], 'obsTimeOfIscStroke', 'isIscStroke')


def by_ldlc_70():
    print_hr(table[table[ldlc] < 70], 'timeOfMACE', 'isMACE')
    print_hr(table[table[ldlc] >= 70], 'timeOfMACE', 'isMACE')
    print_hr(table[table[ldlc] < 70], 'timeOfCvDeath', 'isCvDeath')
    print_hr(table[table[ldlc] >= 70], 'timeOfCvDeath', 'isCvDeath')
    print_hr(table[table[ldlc] < 70], 'obsTimeOfMI', 'isMI')
    print_hr(table[table[ldlc] >= 70], 'obsTimeOfMI', 'isMI')
    print_hr(table[table[ldlc] < 70], 'obsTimeOfIscStroke', 'isIscStroke')
    print_hr(table[table[ldlc] >= 70], 'obsTimeOfIscStroke', 'isIscStroke')


def by_lpa():
    print_hr(table[table[lpa] <= 37], 'timeOfMACE3Y', 'isMACE3Y')
    print_hr(table[table[lpa] > 37], 'timeOfMACE3Y', 'isMACE3Y')
    print_hr(table[table[lpa] <= 37], 'obsTimeOfMI3Y', 'isMI3Y')
    print_hr(table[table[lpa] > 37], 'obsTimeOfMI3Y', 'isMI3Y')
    print_hr(table[table[lpa] <= 37], 'obsTimeOfIscStroke3Y', 'isIscStroke3Y')
    print_hr(table[table[lpa] > 37], 'obsTimeOfIscStroke3Y', 'isIscStroke3Y')

    print_hr(table[table[lpa] <= 120], 'timeOfMACE3Y', 'isMACE3Y')
    print_hr(table[table[lpa] > 120], 'timeOfMACE3Y', 'isMACE3Y')
    print_hr(table[table[lpa] <= 120], 'obsTimeOfMI3Y', 'isMI3Y')
    print_hr(table[table[lpa] > 120], 'obsTimeOfMI3Y', 'isMI3Y')
    print_hr(table[table[lpa] <= 120], 'obsTimeOfIscStroke3Y', 'isIscStroke3Y')
    print_hr(table[table[lpa] > 120], 'obsTimeOfIscStroke3Y', 'isIscStroke3Y')


overall()
key_subgroups()
by_age()
by_sex()
by_egfr()
by_diabetes()
by_recent_remote_mi()
by_prior_stroke()
by_prior_pad()
by_statin_potency()
by_ldlc_70()
by_lpa()

