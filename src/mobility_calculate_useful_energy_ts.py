# -*- coding: utf-8 -*-
"""* * @Author: Lars Nolting, Christina Kockel, and Aaron Praktiknjo  * @Date: 2020-10-19  *"""

from __future__ import division
import pyomo.environ as pyomo
import pyomo.opt as opt
import numpy as np
import pandas as pd
import sys
import os


dirname = os.path.dirname(__file__)

def transform_data():

    vehicle_codes = ['Pas', 'Mot', 'Tru_sm_35', 'Tru_gr_35', 'Tru_semi', 'Oth', 'Bus']


    # 1. Load required data sets
    # ------------------------------------------------------------------------------------------------------------------------

    print('1. Load required data sets')

    root_path = dirname + '/Input Data/Mobility'
    NUTS3_NUTS2 = pd.read_csv(dirname + '/Input Data/General/NUTS_translation.csv', sep=';', encoding="ISO-8859-1", decimal=',')


    # 2. Preprocess data
    # ------------------------------------------------------------------------------------------------------------------------

    print('2. Preprocess data')

    NUTS3_NUTS2 = NUTS3_NUTS2.drop(['NUTS1'], axis=1)
    NUTS2_regions = sorted(list(set(NUTS3_NUTS2['NUTS2'].values.tolist())))


    # 3. Transform data
    # ------------------------------------------------------------------------------------------------------------------------

    print('3. Transform data')

    driving_profiles = {NUTS2: None for NUTS2 in NUTS2_regions}
    for NUTS2 in NUTS2_regions:
        file_path = root_path + '/' + str(NUTS2) + '_yearly_driving_profiles.csv'
        driving_profiles[NUTS2] = pd.read_csv(file_path, sep=';', encoding="ISO-8859-1", decimal=',')

    nuts_profiles = {code: pd.DataFrame(0, index=list(range(8760)), columns=NUTS2_regions) for code in vehicle_codes}
    for code in vehicle_codes:
        for NUTS2 in NUTS2_regions:
            nuts_profiles[code][NUTS2] = driving_profiles[NUTS2][code]


    # 4. Save intermediate results
    # ------------------------------------------------------------------------------------------------------------------------

    print('4. Save intermediate results')
    for code in vehicle_codes:
        nuts_profiles[code].to_csv(dirname + '/Input Data/Mobility/nuts2_hourly_service_' + code + '_km.csv', sep=';', encoding="ISO-8859-1", decimal=',')
        

def validate_and_scale_ts(bio_fuel_shares, benchmark_values_twh):

    vehicle_translation = {
        'Pas': 'passenger_car',
        'Mot': 'motorcycle',
        'Tru_sm_35': 'truck_light',
        'Tru_gr_35': 'truck_heavy',
        'Tru_semi': 'truck_heavy',
        'Oth': 'other',
    }

    fuel_types = ['diesel', 'gasoline', 'gas', 'electricity']


    # 1. Load required data sets
    # ------------------------------------------------------------------------------------------------------------------------

    print('1. Load required data sets')

    service_profiles = {}
    for key in vehicle_translation.keys():
        service_profiles[key] = pd.read_csv(dirname + '/Input Data/Mobility/nuts2_hourly_service_' + key + '_km.csv', sep=';', encoding="ISO-8859-1", decimal=',', index_col=[0])
        service_profiles[key].index.name = 'hour'

    efficiencies = pd.read_excel(dirname + '/Input Data/Mobility/efficiencies.xlsx', engine="openpyxl", sheet_name='Tabelle1', index_col=[0])

    vehicle_fuel_type_shares = pd.read_excel(dirname + '/Input Data/Mobility/vehicle_fuel_type_shares.xlsx', engine="openpyxl", sheet_name='shares', index_col=[0])

    NUTS3_NUTS2 = pd.read_csv(dirname + '/Input Data/General/NUTS_translation.csv', sep=';', encoding="ISO-8859-1", decimal=',')


    # 2. Preprocess data
    # ------------------------------------------------------------------------------------------------------------------------

    print('2. Preprocess data')

    NUTS3_NUTS2 = NUTS3_NUTS2.drop(['NUTS1'], axis=1)
    NUTS2_regions = sorted(list(set(NUTS3_NUTS2['NUTS2'].values.tolist())))


    # 3. Validate  data and calculate scaling factors
    # ------------------------------------------------------------------------------------------------------------------------

    print('3. Validate  data and calculate scaling factors')

    annual_useful_energy_twh = {}
    annual_final_energy_twh = {}

    for fuel_type in fuel_types:

        annual_useful_energy_twh[fuel_type] = 0.0
        annual_final_energy_twh[fuel_type] = 0.0

        for key in vehicle_translation.keys():

            try:
                fuel_share = vehicle_fuel_type_shares.loc[(vehicle_fuel_type_shares.index == key) & (vehicle_fuel_type_shares.fuel_type==fuel_type), 'vehicle_km_share'].values.tolist()[0]

                efficiency_useful_energy = efficiencies.loc[(efficiencies.index == vehicle_translation[key]) & (efficiencies.fuel_type == fuel_type), 'useful_to_service [kWh/100km]'].values.tolist()[0]
                efficiency_final_energy = efficiencies.loc[(efficiencies.index == vehicle_translation[key]) & (efficiencies.fuel_type == fuel_type), 'final_to_service [kWh/100km]'].values.tolist()[0]
                
                useful_energy_twh = service_profiles[key].sum().sum() / 100 * fuel_share * efficiency_useful_energy / 1e9
                final_energy_twh = service_profiles[key].sum().sum() / 100 * fuel_share * efficiency_final_energy / 1e9

                annual_useful_energy_twh[fuel_type] += useful_energy_twh
                annual_final_energy_twh[fuel_type] += final_energy_twh

            except:
                pass

    annual_useful_energy_twh['biomass'] = annual_useful_energy_twh['diesel']*bio_fuel_shares['diesel'] + annual_useful_energy_twh['gasoline']*bio_fuel_shares['gasoline']
    annual_useful_energy_twh['diesel'] *= (1-bio_fuel_shares['diesel'])
    annual_useful_energy_twh['gasoline'] *= (1-bio_fuel_shares['gasoline'])

    annual_final_energy_twh['biomass'] = annual_final_energy_twh['diesel']*bio_fuel_shares['diesel'] + annual_final_energy_twh['gasoline']*bio_fuel_shares['gasoline']
    annual_final_energy_twh['diesel'] *= (1-bio_fuel_shares['diesel'])
    annual_final_energy_twh['gasoline'] *= (1-bio_fuel_shares['gasoline'])


    # 4. find optimal scaling factors
    # ------------------------------------------------------------------------------------------------------------------------

    print('4. find optimal scaling factors')

    model = pyomo.ConcreteModel()

    # define sets of the LP model  
    model.P = pyomo.Set(initialize=vehicle_translation.keys())      # Profiles
    model.F = pyomo.Set(initialize=['diesel', 'gasoline'])          # Fuel types

    # create decision variables
    model.scaling_factor_p = pyomo.Var(model.P, domain=pyomo.NonNegativeReals)                  # scaling factor for profile p
    model.FE_contribution_p_f = pyomo.Var(model.P, model.F, domain=pyomo.NonNegativeReals)      # contribution of a profile to final energy of fuel type in [TWh]
    model.FE_f = pyomo.Var(model.F, domain=pyomo.NonNegativeReals)                              # final energy of fuel type in [TWh]
    model.FE_deviation_pos_f = pyomo.Var(model.F, domain=pyomo.NonNegativeReals)                # positive deviation in final energy of fuel type in [TWh]
    model.FE_deviation_neg_f = pyomo.Var(model.F, domain=pyomo.NonNegativeReals)                # negative deviation in final energy of fuel type in [TWh]
    
    # add constraints
    
    # define contribution
    def define_contribution(model, p, f):            
        if p in vehicle_fuel_type_shares.index.tolist() and \
            f in vehicle_fuel_type_shares.loc[vehicle_fuel_type_shares.index == p, 'fuel_type'].values.tolist() and \
            vehicle_translation[p] in efficiencies.index.tolist() and \
            f in efficiencies.loc[efficiencies.index == vehicle_translation[p], 'fuel_type'].values.tolist():

            fuel_share = vehicle_fuel_type_shares.loc[(vehicle_fuel_type_shares.index == p) & (vehicle_fuel_type_shares.fuel_type == f), 'vehicle_km_share'].values.tolist()[0]

            efficiency_final_energy = efficiencies.loc[(efficiencies.index == vehicle_translation[p]) & (efficiencies.fuel_type == f), 'final_to_service [kWh/100km]'].values.tolist()[0]

            if fuel_share > 0.0:
            
                vehicle_100km = service_profiles[p].sum().sum() / 100

                # print(p,f, fuel_share, efficiency_final_energy, fuel_share * efficiency_final_energy * vehicle_100km)

                return model.scaling_factor_p[p] * fuel_share * efficiency_final_energy * vehicle_100km / 1e9 == model.FE_contribution_p_f[p,f]
            else:
                return model.FE_contribution_p_f[p,f] >= 0
        else:
            return model.FE_contribution_p_f[p,f] >= 0
    model.contribution = pyomo.Constraint(model.P, model.F, rule=define_contribution)

    # sum contributions
    def define_sum_contributions(model, f):
        return sum(model.FE_contribution_p_f[p,f] for p in model.P if vehicle_fuel_type_shares.loc[(vehicle_fuel_type_shares.index == p) & (vehicle_fuel_type_shares.fuel_type == f), 'vehicle_km_share'].values.tolist()[0] > 0.0) * (1-bio_fuel_shares[f])== model.FE_f[f]
    model.sum_contributions = pyomo.Constraint(model.F, rule=define_sum_contributions)

    # define scaling factors 
    def define_scaling_factors(model, f):
        return model.FE_f[f] - benchmark_values_twh[f] == model.FE_deviation_pos_f[f] - model.FE_deviation_neg_f[f] 
    model.scaling_factors = pyomo.Constraint(model.F, rule=define_scaling_factors)

    def define_scaling_factors_larger(model, p):
        return model.scaling_factor_p[p] >= 0.9
    model.scaling_factors_larger = pyomo.Constraint(model.P, rule=define_scaling_factors_larger)

    def define_scaling_factors_smaller(model, p):
        return model.scaling_factor_p[p] <= 1.45
    model.scaling_factors_smaller = pyomo.Constraint(model.P, rule=define_scaling_factors_smaller)

    # define objective function (i.e. dispatch costs)

    def define_objective_function(model):
        return sum(model.FE_deviation_pos_f[f] + model.FE_deviation_neg_f[f] for f in model.F)
    model.Obj = pyomo.Objective(rule=define_objective_function, sense=pyomo.minimize)

    # call solver
    optimizer = opt.SolverFactory('glpk')
    solved_model = optimizer.solve(model, tee=True)

    scaling_factor = {}
    for p in model.P:
        print("Scaling factor of %s is %f" % (p, pyomo.value(model.scaling_factor_p[p])))
        scaling_factor[p] = pyomo.value(model.scaling_factor_p[p])
    for f in model.F:
        print("Deviation in %s is %f TWh" % (f, pyomo.value(model.FE_deviation_pos_f[f] + model.FE_deviation_neg_f[f])))


    # 5. Validate  data and calculate scaling factors
    # ------------------------------------------------------------------------------------------------------------------------

    print('5. Validate  data and calculate scaling factors')

    annual_useful_energy_twh = {}
    annual_final_energy_twh = {}

    for fuel_type in fuel_types:

        annual_useful_energy_twh[fuel_type] = 0.0
        annual_final_energy_twh[fuel_type] = 0.0

        for key in vehicle_translation.keys():

            try:
                fuel_share = vehicle_fuel_type_shares.loc[(vehicle_fuel_type_shares.index == key) & (vehicle_fuel_type_shares.fuel_type==fuel_type), 'vehicle_km_share'].values.tolist()[0]

                efficiency_useful_energy = efficiencies.loc[(efficiencies.index == vehicle_translation[key]) & (efficiencies.fuel_type == fuel_type), 'useful_to_service [kWh/100km]'].values.tolist()[0]
                efficiency_final_energy = efficiencies.loc[(efficiencies.index == vehicle_translation[key]) & (efficiencies.fuel_type == fuel_type), 'final_to_service [kWh/100km]'].values.tolist()[0]
                
                useful_energy_twh = service_profiles[key].sum().sum() / 100 * fuel_share * efficiency_useful_energy / 1e9
                final_energy_twh = scaling_factor[key] * service_profiles[key].sum().sum() / 100 * fuel_share * efficiency_final_energy / 1e9

                annual_useful_energy_twh[fuel_type] += useful_energy_twh
                annual_final_energy_twh[fuel_type] += final_energy_twh

            except:
                pass

    annual_useful_energy_twh['biomass'] = annual_useful_energy_twh['diesel']*bio_fuel_shares['diesel'] + annual_useful_energy_twh['gasoline']*bio_fuel_shares['gasoline']
    annual_useful_energy_twh['diesel'] *= (1-bio_fuel_shares['diesel'])
    annual_useful_energy_twh['gasoline'] *= (1-bio_fuel_shares['gasoline'])

    annual_final_energy_twh['biomass'] = annual_final_energy_twh['diesel']*bio_fuel_shares['diesel'] + annual_final_energy_twh['gasoline']*bio_fuel_shares['gasoline']
    annual_final_energy_twh['diesel'] *= (1-bio_fuel_shares['diesel'])
    annual_final_energy_twh['gasoline'] *= (1-bio_fuel_shares['gasoline'])


    # 6. Scale and save service time series
    # ------------------------------------------------------------------------------------------------------------------------

    print('6. Scale and save service time series')

    for key in vehicle_translation.keys():
        service_profiles[key] = service_profiles[key].multiply(scaling_factor[key])
        service_profiles[key].to_csv(dirname + '/../Final Data/Mobility/nuts2_hourly_mob_service_' + key + '_km.csv', sep=';', encoding="ISO-8859-1", decimal=',')


    # 7. Calculate and save useful energy time series
    # ------------------------------------------------------------------------------------------------------------------------

    print('7. Calculate and save useful energy time series')

    hourly_useful_energy_kw = pd.DataFrame(0.0, index=list(range(8760)), columns=NUTS2_regions)
    hourly_useful_energy_kw.index.name = 'hour'

    for fuel_type in fuel_types:

        for key in vehicle_translation.keys():

            for NUTS2 in NUTS2_regions:

                try:
                    fuel_share = vehicle_fuel_type_shares.loc[(vehicle_fuel_type_shares.index == key) & (vehicle_fuel_type_shares.fuel_type==fuel_type), 'vehicle_km_share'].values.tolist()[0]

                    efficiency_useful_energy = efficiencies.loc[(efficiencies.index == vehicle_translation[key]) & (efficiencies.fuel_type == fuel_type), 'useful_to_service [kWh/100km]'].values.tolist()[0]
                    
                    useful_energy_kwh = service_profiles[key][NUTS2].sum() / 100 * fuel_share * efficiency_useful_energy

                    NUTS2_ts = service_profiles[key][NUTS2].divide(service_profiles[key][NUTS2].sum()).multiply(useful_energy_kwh)

                    hourly_useful_energy_kw[NUTS2] = hourly_useful_energy_kw[NUTS2] + NUTS2_ts

                except:
                    pass

    hourly_useful_energy_kw.to_csv(dirname + '/../Final Data/Mobility/nuts2_hourly_mob_Mechanical_kw.csv', sep=';', encoding="ISO-8859-1", decimal=',')


if __name__ == "__main__":
    
    print('Transform data')
    print()
    transform_data()
    
    benchmark_values_twh = {
        'diesel': 382.5289,
        'gasoline': 203.1603,
        'gas': 5.9308,
        'electricity': 0.2289,
        'biomass': 31.3239
    }

    bio_fuel_shares = {
        'diesel': 0.0584,             
        'gasoline': 0.0429
    }

    print('Validate and scale time series')
    print()
    validate_and_scale_ts(bio_fuel_shares, benchmark_values_twh)