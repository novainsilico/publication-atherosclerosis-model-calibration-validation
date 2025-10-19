"""
Processing raw calibration-simulation output data
and produces CSVs for making Kaplan-Meier curves for each endpoint x subgroup
"""
import os
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter

kmf = KaplanMeierFitter()
study_name = 'FOURIER'

def gen_data_FOURIER():
    global base_folder, full_scalars, output_path, arms, arm_names
    # always put the placebo arm first
    arms = ['FOURIER-placebo', 'FOURIER-evolocumab']
    arm_names = ['Placebo', 'Evolocumab']

    base_folder = 'vpop-event-data'
    input_path = os.path.join(base_folder, 'scalars-trimmed-for-events.csv')
    output_path = os.path.join(base_folder)

    full_scalars = pd.read_csv(input_path)

    gen_data_entire_population()
    mace_by_age()
    mace_by_sex()
    mace_by_egfr()
    mace_by_diabetes()
    mace_by_recent_remote_mi()
    mace_stroke_by_prior_stroke()
    mace_male_by_prior_pad()
    mace_by_prior_pad_co_history()
    mace_by_ldlc_70()
    mace_by_max_statin()
    mace_by_lpa_quartiles()
    mace_by_lpa_37_120()


def gen_km_table(simu_table, events, ref_paths, output_name):
    long_table = pd.DataFrame()
    for event, (event_scalar, time_scalar) in events.items():
        event_table = pd.DataFrame()

        for i, arm in enumerate(arms):
            event_happened = simu_table[f'{arm}_{event_scalar}'] == 1
            event_time = simu_table[f'{arm}_{time_scalar}']

            kmf.fit(event_time, event_observed=event_happened)
            table = kmf.cumulative_density_.reset_index()
            table['Arm'] = arm_names[i]
            table['Origin'] = 'Simulated data (Vpop)'
            event_table = pd.concat([table, event_table], ignore_index=True)

        if event in ref_paths:
            ref_table = pd.read_csv(os.path.join(base_folder, ref_paths[event]))
            ref_table.rename({'x': 'timeline', 'y': 'KM_estimate'}, axis=1, inplace=True)
            ref_table.replace({arm: name for arm, name in zip(arms, arm_names)}, inplace=True)
        else:
            ref_table = pd.DataFrame({'timeline': np.nan, 'KM_estimate': np.nan, 'Arm': arm_names})
        ref_table['Origin'] = f'Observed data ({study_name})'
        event_table = pd.concat([ref_table, event_table], ignore_index=True)
        event_table['Event'] = event
        print(event_table.tail(1))
        long_table = pd.concat([event_table, long_table], ignore_index=True)
    long_table.to_csv(os.path.join(output_path, output_name), index=False)
    print('Done writing', output_name)


def gen_data_entire_population():
    events = {'MACE': ['isMACE3Y', 'timeOfMACE3Y'],
              'CV death': ['isCvDeath3Y', 'timeOfCvDeath3Y'],
              'Fatal or nonfatal MI': ['isMI3Y', 'obsTimeOfMI3Y'],
              'Fatal or nonfatal ischemic stroke': ['isIscStroke3Y', 'obsTimeOfIscStroke3Y'],
              'MALE': ['isMajorAdverseLimbEvent2_5Y', 'obsTimeOfMajorAdverseLimbEvent2_5Y']}
    ref_paths = {'MACE': 'literature-KM-data/kmScaledData-FOURIER-Sabatine-2017-MACE.csv',
                 'Fatal or nonfatal ischemic stroke': 'literature-KM-data/kmScaledData-FOURIER-Giugliano-2020-IscStroke.csv',
                 'MALE': 'literature-KM-data/kmScaledData-FOURIER-Bonaca-2018-MALE.csv'}

    gen_km_table(full_scalars, events, ref_paths, 'entire-population-all-endpoints.csv')


def mace_by_age():
    events = {'MACE': ['isMACE3Y', 'timeOfMACE3Y']}
    ref_paths = {}
    sub_table = full_scalars[full_scalars["ageInit"] < 56]
    gen_km_table(sub_table, events, ref_paths, 'MACE-age-below-56.csv')
    sub_table = full_scalars[(full_scalars["ageInit"] < 63) & (full_scalars["ageInit"] >= 56)]
    gen_km_table(sub_table, events, ref_paths, 'MACE-age-56-to-63.csv')
    sub_table = full_scalars[(full_scalars["ageInit"] < 69) & (full_scalars["ageInit"] >= 63)]
    gen_km_table(sub_table, events, ref_paths, 'MACE-age-63-to-69.csv')
    sub_table = full_scalars[full_scalars["ageInit"] >= 69]
    gen_km_table(sub_table, events, ref_paths, 'MACE-age-above-69.csv')


def mace_by_sex():
    events = {'MACE': ['isMACE3Y', 'timeOfMACE3Y']}
    ref_paths = {}
    sub_table = full_scalars[full_scalars["isMaleSex"] == 1]
    gen_km_table(sub_table, events, ref_paths, 'MACE-male.csv')
    sub_table = full_scalars[full_scalars["isMaleSex"] == 0]
    gen_km_table(sub_table, events, ref_paths, 'MACE-female.csv')


def mace_by_egfr():
    events = {'MACE': ['isMACE2_5Y', 'timeOfMACE2_5Y']}
    ref_paths = {}
    sub_table = full_scalars[full_scalars["eGFR"] < 60]
    gen_km_table(sub_table, events, ref_paths, 'MACE-eGFR-below-60.csv')
    sub_table = full_scalars[(full_scalars["eGFR"] < 90) & (full_scalars["eGFR"] >= 60)]
    gen_km_table(sub_table, events, ref_paths, 'MACE-eGFR-60-to-90.csv')
    sub_table = full_scalars[full_scalars["eGFR"] >= 90]
    gen_km_table(sub_table, events, ref_paths, 'MACE-eGFR-above-90.csv')


def mace_by_diabetes():
    events = {'MACE': ['isMACE3Y', 'timeOfMACE3Y']}
    ref_paths = {'MACE': 'literature-KM-data/kmScaledData-FOURIER-Sabatine-diab-MACE-Diabetes.csv'}
    sub_table = full_scalars[full_scalars["isDiabetes"] == 1]
    gen_km_table(sub_table, events, ref_paths, 'MACE-diabetes.csv')
    ref_paths = {'MACE': 'literature-KM-data/kmScaledData-FOURIER-Sabatine-diab-MACE-NoDiabetes.csv'}
    sub_table = full_scalars[full_scalars["isDiabetes"] == 0]
    gen_km_table(sub_table, events, ref_paths, 'MACE-no-diabetes.csv')


def mace_by_recent_remote_mi():
    events = {'MACE': ['isMACE3Y', 'timeOfMACE3Y']}
    ref_paths = {'MACE': 'literature-KM-data/kmScaledData-FOURIER-Gencer-2020-MACE-RecentMI.csv'}
    sub_table = full_scalars[(full_scalars["isEventInPastYearAtStart"] == 1) & (full_scalars["isPriorChd"] == 1)]
    gen_km_table(sub_table, events, ref_paths, 'MACE-recent-mi.csv')
    ref_paths = {'MACE': 'literature-KM-data/kmScaledData-FOURIER-Gencer-2020-MACE-RemoteMI.csv'}
    sub_table = full_scalars[(full_scalars["isEventInPastYearAtStart"] == 0) & (full_scalars["isPriorChd"] == 1)]
    gen_km_table(sub_table, events, ref_paths, 'MACE-remote-mi.csv')


def mace_stroke_by_prior_stroke():
    events = {'MACE': ['isMACE3Y', 'timeOfMACE3Y'],
              'Fatal or nonfatal ischemic stroke': ['isIscStroke3Y', 'obsTimeOfIscStroke3Y']}
    ref_paths = {}
    sub_table = full_scalars[full_scalars["isPriorCbvd"] == 1]
    gen_km_table(sub_table, events, ref_paths, 'MACE-stroke-prior-stroke.csv')
    sub_table = full_scalars[full_scalars["isPriorCbvd"] == 0]
    gen_km_table(sub_table, events, ref_paths, 'MACE-stroke-no-prior-stroke.csv')


def mace_male_by_prior_pad():
    events = {'MACE': ['isMACE2_5Y', 'timeOfMACE2_5Y'],
              'MALE': ['isMajorAdverseLimbEvent2_5Y', 'obsTimeOfMajorAdverseLimbEvent2_5Y']}
    ref_paths = {'MACE': 'literature-KM-data/kmScaledData-FOURIER-Bonaca-2018-MACE.csv',
                 'MALE': 'literature-KM-data/kmScaledData-FOURIER-Bonaca-2018-MALE-PAD.csv'}
    sub_table = full_scalars[full_scalars["isPriorPad"] == 1]
    gen_km_table(sub_table, events, ref_paths, 'MACE-MALE-prior-pad.csv')
    ref_paths = {'MACE': 'literature-KM-data/kmScaledData-FOURIER-Bonaca-2018-MACE-NoPAD.csv',
                 'MALE':'literature-KM-data/kmScaledData-FOURIER-Bonaca-2018-MALE-NoPAD.csv'}
    sub_table = full_scalars[full_scalars["isPriorPad"] == 0]
    gen_km_table(sub_table, events, ref_paths, 'MACE-MALE-no-prior-pad.csv')


def mace_by_prior_pad_co_history():
    events = {'MACE': ['isMACE2_5Y', 'timeOfMACE2_5Y']}
    ref_paths = {}
    sub_table = full_scalars[(full_scalars["isPriorPad"] == 1) & (full_scalars["bedInvolvement"] == 1)]
    gen_km_table(sub_table, events, ref_paths, 'MACE-prior-pad-only.csv')
    ref_paths = {}
    sub_table = full_scalars[(full_scalars["isPriorPad"] == 1) & (full_scalars["bedInvolvement"] > 1)]
    gen_km_table(sub_table, events, ref_paths, 'MACE-prior-pad-with-mi-or-stroke.csv')


def mace_by_ldlc_70():
    events = {'MACE': ['isMACE2_5Y', 'timeOfMACE2_5Y']}
    ref_paths = {'MACE': 'literature-KM-data/kmScaledData-FOURIER-Giugliano-2017-MACE-LDLcBelow70.csv'}
    sub_table = full_scalars[full_scalars["FOURIER-placebo_baselineLDLcUncorrectedMass"] < 70]
    gen_km_table(sub_table, events, ref_paths, 'MACE-ldlc-below-70.csv')
    ref_paths = {}
    sub_table = full_scalars[full_scalars["FOURIER-placebo_baselineLDLcUncorrectedMass"] >= 70]
    gen_km_table(sub_table, events, ref_paths, 'MACE-ldlc-above-70.csv')


def mace_by_max_statin():
    events = {'MACE': ['isMACE2_5Y', 'timeOfMACE2_5Y']}
    ref_paths = {'MACE': 'literature-KM-data/kmScaledData-FOURIER-Giugliano-2017-MACE-MaxStatins.csv'}
    sub_table = full_scalars[full_scalars["FOURIER-placebo_isAva80Rva40"] == 1]
    gen_km_table(sub_table, events, ref_paths, 'MACE-max-statin.csv')
    ref_paths = {}
    sub_table = full_scalars[full_scalars["FOURIER-placebo_isAva80Rva40"] == 0]
    gen_km_table(sub_table, events, ref_paths, 'MACE-submax-statin.csv')


def mace_by_lpa_quartiles():
    events = {'MACE': ['isMACE3Y', 'timeOfMACE3Y']}
    ref_paths = {}
    sub_table = full_scalars[full_scalars["FOURIER-placebo_baselineLpaSubstance"] < 13]
    gen_km_table(sub_table, events, ref_paths, 'MACE-lpa-below-13.csv')
    sub_table = full_scalars[(full_scalars["FOURIER-placebo_baselineLpaSubstance"] < 37) & (full_scalars["FOURIER-placebo_baselineLpaSubstance"] >= 13)]
    gen_km_table(sub_table, events, ref_paths, 'MACE-lpa-13-to-37.csv')
    sub_table = full_scalars[(full_scalars["FOURIER-placebo_baselineLpaSubstance"] < 165) & (full_scalars["FOURIER-placebo_baselineLpaSubstance"] >= 37)]
    gen_km_table(sub_table, events, ref_paths, 'MACE-lpa-37-to-165.csv')
    sub_table = full_scalars[full_scalars["FOURIER-placebo_baselineLpaSubstance"] >= 165]
    gen_km_table(sub_table, events, ref_paths, 'MACE-lpa-above-165.csv')


def mace_by_lpa_37_120():
    events = {'MACE': ['isMACE3Y', 'timeOfMACE3Y']}
    ref_paths = {}
    sub_table = full_scalars[full_scalars["FOURIER-placebo_baselineLpaSubstance"] <= 37]
    gen_km_table(sub_table, events, ref_paths, 'MACE-lpa-below-37.csv')
    sub_table = full_scalars[full_scalars["FOURIER-placebo_baselineLpaSubstance"] > 37]
    gen_km_table(sub_table, events, ref_paths, 'MACE-lpa-above-37.csv')
    sub_table = full_scalars[full_scalars["FOURIER-placebo_baselineLpaSubstance"] <= 120]
    gen_km_table(sub_table, events, ref_paths, 'MACE-lpa-below-120.csv')
    sub_table = full_scalars[full_scalars["FOURIER-placebo_baselineLpaSubstance"] > 120]
    gen_km_table(sub_table, events, ref_paths, 'MACE-lpa-above-120.csv')


gen_data_FOURIER()
