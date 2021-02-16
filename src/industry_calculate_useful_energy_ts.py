# -*- coding: utf-8 -*-
"""* * @Author: Lars Nolting, Christina Kockel, and Aaron Praktiknjo  * @Date: 2020-10-19  *"""

import pandas as pd
import numpy as np
import sys
import os


dirname = os.path.dirname(__file__)

def calculate_ts():

    industrial_subsectors = ['Food', 'Glass, Ceramics, Stones', 'Automotive', 'Chemical', 'Paper', 'Mechanical Engineering', 'Iron, Steel', 'Other Industry']

    end_use_translation = ['Space Heat','Hot Water','Process Heat-Direct','Space Cooling','Process Cooling','Mechanical','Information','Light']

    final_energy_translation = ['Liquid fuel','Gas','Electricity','Heat','Coal','Solid biomass & Waste']


    # 1. Load required data sets
    # ------------------------------------------------------------------------------------------------------------------------

    print('1. Load required data sets')

    sectoral_regional_factors = pd.read_excel(dirname + '/Input Data/Industry/sectoral_regional_factors_NUTS2.xlsx', engine="openpyxl", sheet_name='sectoral_regional_factors_NUTS2', index_col=[0])

    industrial_standard_load_profiles = pd.read_excel(dirname + '/Input Data/Industry/industrial_standard_load_profiles.xlsx', engine="openpyxl", sheet_name='industrial_standard_load_profil', index_col=[0]).reset_index(drop=True)

    temperature_scaling = pd.read_excel(dirname + '/Input Data/General/temperature_scaling_NUTS2.xlsx', engine="openpyxl", sheet_name='temperature_scaling_NUTS2', index_col=[0]).reset_index(drop=True)

    ageb_end_use = pd.read_excel(dirname + '/Input Data/General/ageb_end_use_2019.xlsx', engine="openpyxl", sheet_name='ind', skiprows=[0,1,2], index_col=[0])

    end_use_con_eff = pd.read_excel(dirname + '/Input Data/General/end_use_conversion_efficiencies.xlsx', engine="openpyxl", sheet_name='Conversion Efficiencies')

    NUTS2_NUTS2 = pd.read_csv(dirname + '/Input Data/General/NUTS_translation.csv', sep=';', encoding="ISO-8859-1", decimal=',')


    # 2. Preprocess data
    # ------------------------------------------------------------------------------------------------------------------------

    print('2. Preprocess data')

    end_use_con_eff = end_use_con_eff.loc[end_use_con_eff.Sector == 'Industry']

    NUTS2_NUTS2 = NUTS2_NUTS2.drop(['NUTS1'], axis=1)
    NUTS2_regions = sorted(list(set(NUTS2_NUTS2['NUTS2'].values.tolist())))

    pj_2_kwh = (1e12 / 3600)


    # 3. Calculate annual final energy per end use and sector for NUTS 2 regions
    # ------------------------------------------------------------------------------------------------------------------------

    print('3. Calculate annual final energy per end use and sector for NUTS 2 regions')

    annual_FE_per_end_use_and_sector_NUTS2 = {final_energy: {sector: pd.DataFrame(0.0, index=NUTS2_regions, columns=end_use_translation) for sector in industrial_subsectors} for final_energy in final_energy_translation}

    for final_energy in final_energy_translation:

        for end_use in end_use_translation:

            total_FE_per_end_use_kwh = ageb_end_use.at[final_energy, end_use]  * pj_2_kwh

            for sector in industrial_subsectors:

                for NUTS2 in NUTS2_regions:

                    share = sectoral_regional_factors.at[NUTS2, sector]

                    annual_FE_per_end_use_and_sector_NUTS2[final_energy][sector].at[NUTS2, end_use] = share*total_FE_per_end_use_kwh


    # 4. Calculate time series of final energy per end use for NUTS 2 regions
    # ------------------------------------------------------------------------------------------------------------------------

    print('4. Calculate time series of final energy per end use for NUTS 2 regions')

    hourly_FE_per_end_use_NUTS2_kw = {final_energy: {end_use: pd.DataFrame(0.0, index=list(range(8760)), columns=NUTS2_regions) for end_use in end_use_translation} for final_energy in final_energy_translation}

    for final_energy in final_energy_translation:

        print('\t', final_energy)

        for end_use in end_use_translation:        
            
            for sector in industrial_subsectors:

                for NUTS2 in NUTS2_regions:

                    total_FE_per_end_use_kwh = annual_FE_per_end_use_and_sector_NUTS2[final_energy][sector].at[NUTS2, end_use]
                    hourly_FE_per_end_use_kw = industrial_standard_load_profiles[sector].multiply(total_FE_per_end_use_kwh)

                    hourly_FE_per_end_use_NUTS2_kw[final_energy][end_use][NUTS2] = hourly_FE_per_end_use_NUTS2_kw[final_energy][end_use][NUTS2] + \
                        hourly_FE_per_end_use_kw


    # 5. Calculate time series of useful energy per end use for NUTS 2 regions
    # ------------------------------------------------------------------------------------------------------------------------

    print('5. Calculate time series of useful energy per end use for NUTS 2 regions')

    hourly_UE_NUTS2_kw = {useful_energy: pd.DataFrame(0.0, index=list(range(8760)), columns=NUTS2_regions) for useful_energy in end_use_translation}

    for final_energy in final_energy_translation:

        for end_use in end_use_translation:     

            if hourly_FE_per_end_use_NUTS2_kw[final_energy][end_use].sum().sum() > 0.0:

                efficiency = end_use_con_eff.loc[(end_use_con_eff['Fuel'] == final_energy) & (end_use_con_eff['End Use'] == end_use), 'Efficiency'].values.tolist()
                efficiency = sum(efficiency) / len(efficiency)

                for NUTS2 in NUTS2_regions:
        
                    hourly_UE_NUTS2_kw[end_use][NUTS2] = hourly_UE_NUTS2_kw[end_use][NUTS2] + hourly_FE_per_end_use_NUTS2_kw[final_energy][end_use][NUTS2].multiply(efficiency)


    # 6. Rescale space heat profiles
    # ------------------------------------------------------------------------------------------------------------------------

    print('6. Rescale space heat profiles')

    total_UE_for_space_heat_before_kw = hourly_UE_NUTS2_kw['Space Heat'].sum().sum()

    for NUTS2 in NUTS2_regions:

        hourly_UE_NUTS2_kw['Space Heat'][NUTS2] = hourly_UE_NUTS2_kw['Space Heat'][NUTS2] * temperature_scaling[NUTS2]

    total_UE_for_space_heat_after_kwh = hourly_UE_NUTS2_kw['Space Heat'].sum().sum()

    scaling_factor = total_UE_for_space_heat_before_kw / total_UE_for_space_heat_after_kwh

    hourly_UE_NUTS2_kw['Space Heat'] = hourly_UE_NUTS2_kw['Space Heat'].multiply(scaling_factor)


    # 7. Save results
    # ------------------------------------------------------------------------------------------------------------------------

    print('7. Save results')

    for end_use in end_use_translation:      
        hourly_UE_NUTS2_kw[end_use].index.name = 'hour'
        hourly_UE_NUTS2_kw[end_use].to_csv(dirname + '/../Final Data/Industry/nuts2_hourly_ind_' + end_use + '_kw.csv', sep=';', encoding="ISO-8859-1", decimal=',')


if __name__ == "__main__":
    calculate_ts()