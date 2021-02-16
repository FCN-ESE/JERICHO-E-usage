# -*- coding: utf-8 -*-
"""* * @Author: Lars Nolting, Christina Kockel, and Aaron Praktiknjo  * @Date: 2020-10-19  *"""

import pandas as pd
import numpy as np
import itertools
from statistics import mean 
import sys
import os


dirname = os.path.dirname(__file__)

# Abbreviations according to VDI 4655:
# Typical days 
# Workday W, Sunday S
# Serene S Clouded C
# Transition T
# Transition: TWS=Transition_Workday_Serene TWC=Transition_Workday_Clouded TSS=Transition_Sunday_Serene TSC=Transition_Sunday_Clouded
# Summer: SWX=Summer_Workday_None Unterscheidung SSX=Summer_Sunday_None Unterscheidung
# Winter: WWS=Winter_Workday_Serene WWC=Winter_Workday_Clouded WSS=Winter_Sunday_Serene WSC=Winter_Sunday_Clouded

def create_randomness(NUTS3_regions):
	
	p = [0.25, 0.5, 0.25] # distribution for shift in hours

	daily_energy_demand_scaling_helper = {}
	day_shift_helper = {}

	for NUTS3 in NUTS3_regions:

		daily_energy_demand_scaling_helper[NUTS3] = [np.random.normal(loc=1, scale=(0.15/3)*1) for i in range(365)]
		daily_energy_demand_scaling_helper[NUTS3] = [value for value in daily_energy_demand_scaling_helper[NUTS3] for i in range(24)]

		day_shift_helper[NUTS3] = [np.random.choice(np.arange(-1, 2, 1), p=p) for i in range(365)]
		day_shift_helper[NUTS3] = [value for value in day_shift_helper[NUTS3] for i in range(24)]

	daily_energy_demand_scaling = pd.DataFrame.from_dict(daily_energy_demand_scaling_helper)
	day_shift = pd.DataFrame.from_dict(day_shift_helper)

	return daily_energy_demand_scaling, day_shift


def change_order_of_hours(df, steps):
    ### steps can be postiive (moving up) or negative (moving down)
    df_copy = df.copy()

    index_list = df_copy.index.tolist()
    index_list = [index - steps for index in index_list]
    index_list = [index - 24 if index > 24 else index for index in index_list]
    index_list = [index + 24  if index < 1 else index for index in index_list]

    df_copy = df_copy.reindex(index_list).reset_index(drop=True)
    df_copy.index = df_copy.index + 1

    return df_copy


def calculate_intermediate_ts():

    # 1. Load required data sets
    # ------------------------------------------------------------------------------------------------------------------------

    print('1. Load required data sets')

    NUTS3_building_data = pd.read_csv(dirname + '/Input Data/Residential/NUTS3_building_data.csv', sep=';', encoding="ISO-8859-1", decimal=',', index_col=[0])
    NUTS3_climate_zone = pd.read_excel(dirname + '/Input Data/Residential/NUTS3_climate_zone.xlsx', engine="openpyxl", sheet_name='NUTS3')
    P_SFH_MFH = pd.read_excel(dirname + '/Input Data/Residential/Typical_building_data.xlsx', engine="openpyxl", sheet_name='P_SFH_MFH', index_col='Building type')
    TD_factors = pd.read_excel(dirname + '/Input Data/Residential/Typical_building_data.xlsx', engine="openpyxl", sheet_name='TD_factors', header=[0,1], index_col=[0])
    TD_SFH_MFH_SHD_Sk_base = pd.read_excel(dirname + '/Input Data/Residential/Typical_building_data.xlsx', engine="openpyxl", sheet_name='TD_SFH_MFH_SHD_Sk', header=[0,1], index_col=[0])
    TD_SFH_MFH_HWD_Sk_base = pd.read_excel(dirname + '/Input Data/Residential/Typical_building_data.xlsx', engine="openpyxl", sheet_name='TD_SFH_MFH_HWD_Sk', header=[0,1], index_col=[0])
    TD_SFH_MFH_OTH_Sk_base = pd.read_excel(dirname + '/Input Data/Residential/Typical_building_data.xlsx', engine="openpyxl", sheet_name='TD_SFH_MFH_OTH_Sk', header=[0,1], index_col=[0])
    public_holidays = pd.read_excel(dirname + '/Input Data/General/Public_Holidays_2019.xlsx', engine="openpyxl", sheet_name='Public_Holidays')
    nuts2_time_series_header = pd.read_excel(dirname + '/Input Data/Residential/NUTS2_time_series_header.xlsx', engine="openpyxl", sheet_name='header', index_col=[0])
    temperature_time_series = pd.read_excel(dirname + '/Input Data/Residential/time_series_climate_zone_2015.xlsx', engine="openpyxl", sheet_name='temperature')
    cloud_coverage_time_series = pd.read_excel(dirname + '/Input Data/Residential/time_series_climate_zone_2015.xlsx', engine="openpyxl", sheet_name='cloud coverage')
    TD_time_series = pd.read_excel(dirname + '/Input Data/Residential/time_series_climate_zone_2015.xlsx', engine="openpyxl", sheet_name='typical day')
    NUTS3_NUTS2 = pd.read_csv(dirname + '/Input Data/General/NUTS_translation.csv', sep=';', encoding="ISO-8859-1", decimal=',')


    # 2. Preprocess data
    # ------------------------------------------------------------------------------------------------------------------------

    print('2. Preprocess data')

    NUTS3_NUTS2 = NUTS3_NUTS2.drop(['NUTS1'], axis=1)
    NUTS2s = sorted(list(set(NUTS3_NUTS2['NUTS2'].values.tolist())))
    NUTS3_regions = sorted(list(set(NUTS3_NUTS2['NUTS3'].values.tolist())))

    daily_energy_demand_scaling, day_shift = create_randomness(NUTS3_regions)

    NUTS3_climate_zone_columns = ['NUTS3', 'climate_zone']
    NUTS3_climate_zone = NUTS3_climate_zone[NUTS3_climate_zone_columns]

    P_SFH_MFH_columns = ['Average number of persons per apartment']
    P_SFH_MFH = P_SFH_MFH[P_SFH_MFH_columns]

    public_holidays = public_holidays['DD_cont'].values.tolist()
    public_holidays = [entry+1 for entry in public_holidays]    # +1 because 2019 the year begins with Tuesday


    # 3. Create shifted variants of the daily scaling factors
    # ------------------------------------------------------------------------------------------------------------------------

    print('3. Create shifted variants of the daily scaling factors')

    TD_SFH_MFH_SHD_Sk = {}
    TD_SFH_MFH_HWD_Sk = {}
    TD_SFH_MFH_OTH_Sk = {}

    TD_SFH_MFH_SHD_Sk['0'] = TD_SFH_MFH_SHD_Sk_base
    TD_SFH_MFH_SHD_Sk['-1'] = change_order_of_hours(TD_SFH_MFH_SHD_Sk['0'], steps=-1)
    TD_SFH_MFH_SHD_Sk['-2'] = change_order_of_hours(TD_SFH_MFH_SHD_Sk['0'], steps=-2)
    TD_SFH_MFH_SHD_Sk['1'] = change_order_of_hours(TD_SFH_MFH_SHD_Sk['0'], steps=1)
    TD_SFH_MFH_SHD_Sk['2'] = change_order_of_hours(TD_SFH_MFH_SHD_Sk['0'], steps=2)

    TD_SFH_MFH_HWD_Sk['0'] = TD_SFH_MFH_HWD_Sk_base
    TD_SFH_MFH_HWD_Sk['-1'] = change_order_of_hours(TD_SFH_MFH_HWD_Sk['0'], steps=-1)
    TD_SFH_MFH_HWD_Sk['-2'] = change_order_of_hours(TD_SFH_MFH_HWD_Sk['0'], steps=-2)
    TD_SFH_MFH_HWD_Sk['1'] = change_order_of_hours(TD_SFH_MFH_HWD_Sk['0'], steps=1)
    TD_SFH_MFH_HWD_Sk['2'] = change_order_of_hours(TD_SFH_MFH_HWD_Sk['0'], steps=2)

    TD_SFH_MFH_OTH_Sk['0'] = TD_SFH_MFH_OTH_Sk_base
    TD_SFH_MFH_OTH_Sk['-1'] = change_order_of_hours(TD_SFH_MFH_OTH_Sk['0'], steps=-1)
    TD_SFH_MFH_OTH_Sk['-2'] = change_order_of_hours(TD_SFH_MFH_OTH_Sk['0'], steps=-2)
    TD_SFH_MFH_OTH_Sk['1'] = change_order_of_hours(TD_SFH_MFH_OTH_Sk['0'], steps=1)
    TD_SFH_MFH_OTH_Sk['2'] = change_order_of_hours(TD_SFH_MFH_OTH_Sk['0'], steps=2)


    # 4. Adjust time series per climate zone
    # ------------------------------------------------------------------------------------------------------------------------

    print('4. Adjust time series per climate zone')

    # 4.1 Classification of days in Workday/Sunday
    ## Add day counter
    temperature_time_series['DD_cont'] = (temperature_time_series.index.values.astype(int)/24).astype(int)+1+1  # +1 because 2019 the year begins with Tuesday
    cloud_coverage_time_series['DD_cont'] = (cloud_coverage_time_series.index.values.astype(int)/24).astype(int)+1+1  # +1 because 2019 the year begins with Tuesday
    ## Add the Workday (WD) column and first define all WD as Workday
    temperature_time_series['WD'] = 'W'
    cloud_coverage_time_series['WD'] = 'W'
    ## Alle Sonntage definieren
    temperature_time_series.loc[temperature_time_series['DD_cont']%7 == 0, 'WD'] = 'S'
    cloud_coverage_time_series.loc[cloud_coverage_time_series['DD_cont']%7 == 0, 'WD'] = 'S'

    # 4.2 Adding holidays as Sundays
    temperature_time_series.loc[temperature_time_series['DD_cont'].isin(public_holidays), 'WD'] = 'S'
    cloud_coverage_time_series.loc[cloud_coverage_time_series['DD_cont'].isin(public_holidays), 'WD'] = 'S'

    # 4.3 Add daily average temperature und define season Summer/Transition/Winter (S/T/W)
    climate_zone_header_t = ['CZ'+str(i)+'_t_[°C]' for i in range(1,16)]
    for header in climate_zone_header_t:
        t_daily_average = list(temperature_time_series[header].groupby(np.arange(len(temperature_time_series[header]))//24).mean())
        temperature_time_series[header + '_daily_average'] = list(itertools.chain.from_iterable([[t]*24 for t in t_daily_average]))

        temperature_time_series[header + '_season'] = 'W'
        temperature_time_series.loc[temperature_time_series[header + '_daily_average']>=5, header + '_season'] = 'T'
        temperature_time_series.loc[temperature_time_series[header + '_daily_average']>15, header + '_season'] = 'S'

    # 4.4 Add daily mean cloud coverage define if serene or covered (S/C) 
    climate_zone_header_N = ['CZ'+str(i)+'_N_[1/8]' for i in range(1,16)]
    for header in climate_zone_header_N:
        N_daily_average = list(cloud_coverage_time_series[header].groupby(np.arange(len(cloud_coverage_time_series[header]))//24).mean())
        cloud_coverage_time_series[header + '_daily_average'] = list(itertools.chain.from_iterable([[t]*24 for t in N_daily_average]))

        cloud_coverage_time_series[header + '_cloud_coverage'] = 'S'
        cloud_coverage_time_series.loc[cloud_coverage_time_series[header + '_daily_average']>=5, header + '_cloud_coverage'] = 'C'

    # 4.5 Classification of days into typical days
    climate_zone_header_TD = ['CZ'+str(i) for i in range(1,16)]
    for header in climate_zone_header_TD:
        season = temperature_time_series[header + '_t_[°C]_season'].values.tolist()
        cloud_coverage = cloud_coverage_time_series[header + '_N_[1/8]_cloud_coverage'].values.tolist()
        WD = temperature_time_series['WD'].values.tolist()
        TD = [season[i] + WD[i] + cloud_coverage[i] if season[i] != "S" else season[i] + WD[i] + 'X' for i in range (len(season))]
        TD_time_series[header] = TD
        temperature_time_series[header] = TD


    # 5. Compile NUTS2 data
    # ------------------------------------------------------------------------------------------------------------------------

    print('5. Compile NUTS2 data')

    nuts2_time_series_shd_kw = nuts2_time_series_header.copy()
    nuts2_time_series_hwd_kw = nuts2_time_series_header.copy()
    nuts2_time_series_oth_kw = nuts2_time_series_header.copy()

    index_counter = 1
    progress_counter = 1
    progress_counter_max = len(NUTS3_building_data.index.tolist())
    for NUTS3, row in NUTS3_building_data.iterrows():

        NUTS2 = NUTS3_NUTS2.loc[NUTS3_NUTS2.NUTS3 == NUTS3, 'NUTS2'].values.tolist()[0]

        # 5.1 Find assigned climate zone for current county
        climate_zone = NUTS3_climate_zone.loc[NUTS3_climate_zone['NUTS3'] == NUTS3, 'climate_zone'].values[0]
        climate_zone_header = 'CZ' + str(climate_zone)

        # 5.2 Adjust factors to current climate zone
        TD_names = ['TWS', 'TWC', 'TSS', 'TSC', 'SWX', 'SSX', 'WWS', 'WWC', 'WSS', 'WSC']
        Number_TD = {TD_name: int(len(TD_time_series.loc[TD_time_series[climate_zone_header] == TD_name, climate_zone_header])/24) for TD_name in TD_names}
        FSH_TD_scaled_SFH = [TD_factors.loc['FSH,TD',('SFH', slice(None))].tolist()[i] * list(Number_TD.values())[i] for i in range(len(list(Number_TD.values())))]
        FSH_TD_scaled_MFH = [TD_factors.loc['FSH,TD',('MFH', slice(None))].tolist()[i] * list(Number_TD.values())[i] for i in range(len(list(Number_TD.values())))]
        FSH_TD_scaled_normalized_SFH = [1/(sum(FSH_TD_scaled_SFH)) * FSH_TD_scaled_SFH[i]/list(Number_TD.values())[i] for i in range(len(FSH_TD_scaled_SFH))]
        FSH_TD_scaled_normalized_MFH = [1/(sum(FSH_TD_scaled_MFH)) * FSH_TD_scaled_MFH[i]/list(Number_TD.values())[i] for i in range(len(FSH_TD_scaled_MFH))]
        TD_factors.loc['FSH,TD,scaled,normalized'] = FSH_TD_scaled_normalized_SFH + FSH_TD_scaled_normalized_MFH
        TD_factors.loc['TD_t_daily_average'] = [mean(temperature_time_series[temperature_time_series[climate_zone_header] == TD_name][climate_zone_header+'_t_[°C]_daily_average']) for TD_name in TD_names] + [0]*10
        
        # 5.3 Calculate daily energy demand depending on typical day 
        SHD_day_SFH = [row['Number SFH'] * row['SHD_SFH'] * TD_factors.at['FSH,TD,scaled,normalized', ('SFH', TD)] for TD in TD_time_series[climate_zone_header]]
        SHD_day_SFH = [SHD_day_SFH[i] * daily_energy_demand_scaling.iloc[i][NUTS3] for i in range(len(SHD_day_SFH))]
        SHD_day_MFH = [row['Number MFH'] * row['SHD_MFH'] * TD_factors.at['FSH,TD,scaled,normalized', ('MFH', TD)] for TD in TD_time_series[climate_zone_header]]
        SHD_day_MFH = [SHD_day_MFH[i] * daily_energy_demand_scaling.iloc[i][NUTS3] for i in range(len(SHD_day_MFH))]

        HWD_day_SFH = [row['Number SFH'] * row['HWD_SFH'] * (1/365 + row['Appartment per SFH'] * P_SFH_MFH.at['SFH', 'Average number of persons per apartment'] * TD_factors.at['FHW,TD', ('SFH', TD)]) for TD in TD_time_series[climate_zone_header]]
        HWD_day_SFH = [HWD_day_SFH[i] * daily_energy_demand_scaling.iloc[i][NUTS3] for i in range(len(HWD_day_SFH))]
        HWD_day_SFH = [HWD_day_SFH[i] if HWD_day_SFH[i] >= 0 else 0 for i in range(len(HWD_day_SFH))]
        HWD_day_MFH = [row['Number MFH'] * row['HWD_MFH'] * (1/365 + row['Appartment per MFH'] * TD_factors.at['FHW,TD', ('MFH', TD)]) for TD in TD_time_series[climate_zone_header]]
        HWD_day_MFH = [HWD_day_MFH[i] if HWD_day_MFH[i] >= 0 else 0 for i in range(len(HWD_day_MFH))]
        HWD_day_MFH = [HWD_day_MFH[i] * daily_energy_demand_scaling.iloc[i][NUTS3] for i in range(len(HWD_day_MFH))]

        OTH_day_SFH = [row['Number SFH'] * row['OTH_SFH'] * (1/365 + row['Appartment per SFH'] * P_SFH_MFH.at['SFH', 'Average number of persons per apartment'] * TD_factors.at['FOTH,TD', ('SFH', TD)]) for TD in TD_time_series[climate_zone_header]]
        OTH_day_SFH = [OTH_day_SFH[i] * daily_energy_demand_scaling.iloc[i][NUTS3] for i in range(len(OTH_day_SFH))]
        OTH_day_SFH = [OTH_day_SFH[i] if OTH_day_SFH[i] >= 0 else 0 for i in range(len(OTH_day_SFH))]
        OTH_day_MFH = [row['Number MFH'] * row['OTH_MFH'] * (1/365 + row['Appartment per MFH'] * TD_factors.at['FOTH,TD', ('MFH', TD)]) for TD in TD_time_series[climate_zone_header]]    
        OTH_day_MFH = [OTH_day_MFH[i] * daily_energy_demand_scaling.iloc[i][NUTS3] for i in range(len(OTH_day_MFH))]
        OTH_day_MFH = [OTH_day_MFH[i] if OTH_day_MFH[i] >= 0 else 0 for i in range(len(OTH_day_MFH))]

        # 5.4 Calculate hourly energy demand
        SHD_hourly_SFH = [SHD_day_SFH[i] * TD_SFH_MFH_SHD_Sk[str(day_shift.iloc[i][NUTS3])].at[TD_time_series['HH'][i], ('SFH', TD_time_series[climate_zone_header][i])] * (273.15+temperature_time_series[climate_zone_header+'_t_[°C]_daily_average'][i])/(273.15+TD_factors.at['TD_t_daily_average', ('SFH', TD_time_series[climate_zone_header][i])]) for i in range(len(SHD_day_SFH))]
        SHD_hourly_MFH = [SHD_day_MFH[i] * TD_SFH_MFH_SHD_Sk[str(day_shift.iloc[i][NUTS3])].at[TD_time_series['HH'][i], ('MFH', TD_time_series[climate_zone_header][i])] * (273.15+temperature_time_series[climate_zone_header+'_t_[°C]_daily_average'][i])/(273.15+TD_factors.at['TD_t_daily_average', ('SFH', TD_time_series[climate_zone_header][i])]) for i in range(len(SHD_day_MFH))]
        HWD_hourly_SFH = [HWD_day_SFH[i] * TD_SFH_MFH_HWD_Sk[str(day_shift.iloc[i][NUTS3])].at[TD_time_series['HH'][i], ('SFH', TD_time_series[climate_zone_header][i])] for i in range(len(HWD_day_SFH))]
        HWD_hourly_MFH = [HWD_day_MFH[i] * TD_SFH_MFH_HWD_Sk[str(day_shift.iloc[i][NUTS3])].at[TD_time_series['HH'][i], ('MFH', TD_time_series[climate_zone_header][i])] for i in range(len(HWD_day_MFH))]
        OTH_hourly_SFH = [OTH_day_SFH[i] * TD_SFH_MFH_OTH_Sk[str(day_shift.iloc[i][NUTS3])].at[TD_time_series['HH'][i], ('SFH', TD_time_series[climate_zone_header][i])] for i in range(len(OTH_day_SFH))]
        OTH_hourly_MFH = [OTH_day_MFH[i] * TD_SFH_MFH_OTH_Sk[str(day_shift.iloc[i][NUTS3])].at[TD_time_series['HH'][i], ('MFH', TD_time_series[climate_zone_header][i])] for i in range(len(OTH_day_MFH))]

        # 5.5 Add time series to NUTS2 regions
        nuts2_time_series_shd_kw[NUTS2] = nuts2_time_series_shd_kw[NUTS2] + SHD_hourly_SFH
        nuts2_time_series_shd_kw[NUTS2] = nuts2_time_series_shd_kw[NUTS2] + SHD_hourly_MFH
        nuts2_time_series_hwd_kw[NUTS2] = nuts2_time_series_hwd_kw[NUTS2] + HWD_hourly_SFH
        nuts2_time_series_hwd_kw[NUTS2] = nuts2_time_series_hwd_kw[NUTS2] + HWD_hourly_MFH
        nuts2_time_series_oth_kw[NUTS2] = nuts2_time_series_oth_kw[NUTS2] + OTH_hourly_SFH
        nuts2_time_series_oth_kw[NUTS2] = nuts2_time_series_oth_kw[NUTS2] + OTH_hourly_MFH

        # print progress
        if progress_counter / progress_counter_max > 0.1:
            print(str(round(100*index_counter/progress_counter_max,0)) + '% done')
            progress_counter = 1
        index_counter += 1
        progress_counter += 1

        
    # 6. Save intermediate results
    # ------------------------------------------------------------------------------------------------------------------------

    print('6. Save intermediate results')

    nuts2_time_series_shd_kw.to_csv(dirname + '/Input Data/Residential/nuts2_hourly_res_Space Heat_kw_intermediate.csv', sep=';', encoding="ISO-8859-1", decimal=',')
    nuts2_time_series_hwd_kw.to_csv(dirname + '/Input Data/Residential/nuts2_hourly_res_Hot Water_kw_intermediate.csv', sep=';', encoding="ISO-8859-1", decimal=',')
    nuts2_time_series_oth_kw.to_csv(dirname + '/Input Data/Residential/nuts2_hourly_res_oth_kw_intermediate.csv', sep=';', encoding="ISO-8859-1", decimal=',')


def split_demand_ts():

    end_use_translation = ['Space Heat','Hot Water','Process Heat-Direct','Space Cooling','Process Cooling','Mechanical','Information','Light']


    # 1. Load required data sets
    # ------------------------------------------------------------------------------------------------------------------------

    print('1. Load required data sets')

    nuts2_hourly_oth_kw = pd.read_csv(dirname + '/Input Data/Residential/nuts2_hourly_res_oth_kw_intermediate.csv', sep=';', encoding="ISO-8859-1", decimal=',', index_col=[0])
    ageb_end_use = pd.read_excel(dirname + '/Input Data/General/ageb_end_use_2019.xlsx', engine="openpyxl", sheet_name='res', skiprows=[0,1,2], index_col=[0])
    end_use_con_eff = pd.read_excel(dirname + '/Input Data/General/end_use_conversion_efficiencies.xlsx', engine="openpyxl", sheet_name='Conversion Efficiencies')


    # 2. Preprocess data
    # ------------------------------------------------------------------------------------------------------------------------

    print('2. Preprocess data')

    end_use_con_eff = end_use_con_eff.loc[end_use_con_eff.Sector == 'Residential']

    pj_2_kwh = (1e12 / 3600)

    # 3. Split other residential demand in end use types
    # ------------------------------------------------------------------------------------------------------------------------

    print('3. Split other residential demand in end use types')

    # reduce efficiency data to required data
    end_use_con_eff_oth_elec = end_use_con_eff.loc[(end_use_con_eff.Fuel == 'Electricity') & ~(end_use_con_eff['End Use'].isin(['Space Heat', 'Hot Water']))]

    # normalize time series
    nuts2_hourly_oth_kw_shares = nuts2_hourly_oth_kw.div(nuts2_hourly_oth_kw.sum().sum(), axis=1)

    # reduce final energy data to required data
    ageb_end_use_elec = ageb_end_use.loc['Electricity']

    end_use_results_kw = {}
    for end_use in end_use_translation:
        if end_use not in ['Space Heat', 'Hot Water']:

            # get efficiency of conversion from final to useful energy
            efficiency = end_use_con_eff_oth_elec.loc[end_use_con_eff_oth_elec['End Use'] == end_use, 'Efficiency'].values[0]

            # get final energy
            final_energy = ageb_end_use_elec[end_use]

            # convert from PJ to kWh
            final_energy = final_energy * pj_2_kwh

            # calculate useful energy
            useful_energy = final_energy * efficiency
            
            # create time series of useful energy demand
            end_use_results_kw[end_use] = nuts2_hourly_oth_kw_shares.copy().multiply(useful_energy)


    # 4. Add gas demand for process heat
    # ------------------------------------------------------------------------------------------------------------------------

    print('4. Add gas demand for process heat')

    efficiency = end_use_con_eff_oth_elec.loc[end_use_con_eff_oth_elec['End Use'] == 'Process Heat-Direct', 'Efficiency'].values[0]

    final_energy = ageb_end_use.at['Gas','Process Heat-Direct']

    final_energy = final_energy * pj_2_kwh

    useful_energy = final_energy * efficiency

    end_use_results_kw['Process Heat-Direct'] = end_use_results_kw['Process Heat-Direct'].add(nuts2_hourly_oth_kw_shares.copy().multiply(useful_energy), fill_value=0)


    # 4. Add mineral oil demand for meachanical energy
    # ------------------------------------------------------------------------------------------------------------------------

    print('4. Add mineral oil demand for meachanical energy')

    efficiency = end_use_con_eff_oth_elec.loc[end_use_con_eff_oth_elec['End Use'] == 'Mechanical', 'Efficiency'].values[0]

    final_energy = ageb_end_use.at['Liquid fuel','Mechanical']

    final_energy = final_energy * pj_2_kwh

    useful_energy = final_energy * efficiency

    end_use_results_kw['Mechanical'] = end_use_results_kw['Mechanical'].add(nuts2_hourly_oth_kw_shares.copy().multiply(useful_energy), fill_value=0)


    # 5. Save results
    # ------------------------------------------------------------------------------------------------------------------------

    print('5. Save results')

    for key, value in end_use_results_kw.items():
        value.to_csv(dirname + '/Input Data/Residential/nuts2_hourly_res_' + key + '_kw_intermediate.csv', sep=';', encoding="ISO-8859-1", decimal=',')


def scale_ts():

    end_use_translation = ['Space Heat','Hot Water','Process Heat-Direct','Space Cooling','Process Cooling','Mechanical','Information','Light']

    final_energy_translation = ['Liquid fuel','Gas','Electricity','Heat','Coal','Solid biomass & Waste']

    # 1. Load required data sets
    # ------------------------------------------------------------------------------------------------------------------------

    print('1. Load required data sets')

    intermediate_data_kw = {}
    for end_use in end_use_translation:
        intermediate_data_kw[end_use] = pd.read_csv(dirname + '/Input Data/Residential/nuts2_hourly_res_' + end_use + '_kw_intermediate.csv', sep=';', encoding="ISO-8859-1", decimal=',', index_col=[0])

    ageb_end_use = pd.read_excel(dirname + '/Input Data/General/ageb_end_use_2019.xlsx', engine="openpyxl", sheet_name='res', skiprows=[0,1,2], index_col=[0])
    end_use_con_eff = pd.read_excel(dirname + '/Input Data/General/end_use_conversion_efficiencies.xlsx', engine="openpyxl", sheet_name='Conversion Efficiencies')


    # 2. Preprocess data
    # ------------------------------------------------------------------------------------------------------------------------

    print('2. Preprocess data')

    end_use_con_eff = end_use_con_eff.loc[end_use_con_eff.Sector == 'Residential']

    pj_2_kwh = (1e12 / 3600)

    # 3. Validate and scale data
    # ------------------------------------------------------------------------------------------------------------------------

    print('3. Validate and scale data')

    for end_use in end_use_translation:

        total_useful_energy_ageb_kwh = 0.0
        total_useful_energy_intermediate_kwh = intermediate_data_kw[end_use].sum().sum()

        for final_energy_carrier in ageb_end_use.index.tolist():
            if ageb_end_use.at[final_energy_carrier, end_use] > 0.0 and final_energy_carrier in final_energy_translation:

                # get efficiency of conversion from final to useful energy
                efficiency = end_use_con_eff.loc[(end_use_con_eff['Fuel'] == final_energy_carrier) & (end_use_con_eff['End Use'] == end_use), 'Efficiency'].values[0]
            
                # get final energy
                final_energy = ageb_end_use.at[final_energy_carrier, end_use]

                # convert from PJ to kWh
                final_energy = final_energy * pj_2_kwh

                # calculate useful energy
                useful_energy = final_energy * efficiency
                
                # create time series of useful energy demand
                total_useful_energy_ageb_kwh += useful_energy

        # scale data
        intermediate_data_kw[end_use] = intermediate_data_kw[end_use].multiply(total_useful_energy_ageb_kwh/total_useful_energy_intermediate_kwh)

    # 4. Save results
    # ------------------------------------------------------------------------------------------------------------------------

    print('4. Save results')

    for end_use in end_use_translation:
        intermediate_data_kw[end_use].to_csv(dirname + '/../Final Data/Residential/nuts2_hourly_res_' + end_use + '_kw.csv', sep=';', encoding="ISO-8859-1", decimal=',')


if __name__ == "__main__":

    print('Calculate intermediate time series')
    print()
    calculate_intermediate_ts()

    print('Split unspecified demand time series')
    print()
    split_demand_ts()

    print('Scale time series')
    print()
    scale_ts()