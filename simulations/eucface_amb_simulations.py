#!/usr/bin/env python

""" EucFACE CO2 simulations

Full spin-up and simulations under amb CO conditions for EucFACE
Include spin-up period, post-industrial period (1750 - 2011),
        observed period (2012 - 2019), and predicted period (2020 - 2069)
Include dry and wet climate conditions
Include NOP, MDP and HIP phosphorus treatment

"""

import os
import shutil
import sys
import subprocess

USER = os.getlogin()
sys.path.append('/Users/%s/Documents/Research/Projects/EucFACE_Modeling/GDAY-EucFACE/scripts' % (USER))
import adjust_gday_param_file as ad
import translate_GDAY_output_to_EUCFACE_format as tr


__author__  = "Mingkai Jiang"
__version__ = "2.0 (06.04.2020)"
__email__   = "m.jiang@westernsydney.edu.au"

def main(experiment_id, site, 
         SPIN_UP=True, POST_INDUST=True, 
         OBS_DRY=True, OBS_WET=True,
         PRD_DRY_NOP=True, PRD_WET_NOP=True, 
         PRD_DRY_MDP=True, PRD_WET_MDP=True, 
         PRD_DRY_HIP=True, PRD_WET_HIP=True):

    GDAY_SPIN = "./gday -s -p "
    GDAY = "./gday -p "

    # dir names
    base_param_name = "base_start_with_P"
    base_param_dir = "/Users/%s/Documents/Research/Projects/EucFACE_Modeling/GDAY-EucFACE/GDAY/params" % (USER)
    base_dir = os.path.dirname(os.getcwd())
    param_dir = os.path.join(base_dir, "params")
    met_dir = os.path.join(base_dir, "met_data")
    run_dir = os.path.join(base_dir, "outputs")

    if SPIN_UP == True:

        # copy base files to make two new experiment files
        shutil.copy(os.path.join(base_param_dir, base_param_name + ".cfg"),
                    os.path.join(param_dir, "%s_%s_model_spinup.cfg" % \
                                                (experiment_id, site)))

        # Run model to equilibrium assuming forest, growing C pools from
        # effectively zero
        itag = "%s_%s_model_spinup" % (experiment_id, site)
        otag = "%s_%s_model_spunup" % (experiment_id, site)
        mtag = "%s_met_spinup_daily_50yrs.csv" % (site)
        out_fn = itag + "_equilib.out"
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)

        replace_dict = {
                        # files
                        "out_param_fname": "%s" % (out_param_fname),
                        "cfg_fname": "%s" % (cfg_fname),
                        "met_fname": "%s" % (met_fname),
                        "out_fname": "%s" % (out_fname),

                        # default C:N 25.
                        # Canopy height = 22 m average of 6 plots at UWS, site_description_stuff/EucFACE_Plot_Summary.doc
                        "activesoil": "0.001",
                        "activesoiln": "0.00004",
                        "activesoilp": "0.000002",
                        "age": "90.0",                # EucFACE parameter file
                        "branch": "0.001",
                        "branchn": "0.00004",
                        "branchp": "0.000002",
                        "cstore": "0.0",
                        "nstore": "0.0",
                        "pstore": "0.0",
                        "inorgn": "0.0000",           # 0.00004
                        "inorglabp": "0.0000",        # 0.00004
                        "inorgsorbp": "0.0",
                        "inorgssorbp": "0.0",
                        "inorgoccp": "0.0",
                        "inorgparp": "0.1",
                        "fertilizerp": "0.0",         # Fertilizer P pool
                        "metabsoil": "0.0",
                        "metabsoiln": "0.0",
                        "metabsoilp": "0.0",
                        "metabsurf": "0.0",
                        "metabsurfn": "0.0",
                        "metabsurfp": "0.0",
                        "passivesoil": "0.001",
                        "passivesoiln": "0.0004",
                        "passivesoilp": "0.000002",
                        "prev_sma": "1.0",
                        "root": "0.001",
                        "croot": "0.0",   # don't simulate coarse roots
                        "crootn": "0.0",  # don't simulate coarse roots
                        "crootp": "0.0",  # don't simulate coarse roots
                        "rootn": "0.00004",
                        "rootp": "0.000002",
                        "sapwood": "0.001",
                        "shoot": "0.001",
                        "shootn": "0.00004",
                        "shootp": "0.000002",
                        "slowsoil": "0.001",
                        "slowsoiln": "0.00004",
                        "slowsoilp": "0.000002",
                        "stem": "0.001",
                        "stemn": "0.00004",
                        "stemp": "0.000002",
                        "stemnimm": "0.00004",
                        "stempimm": "0.000002",
                        "stemnmob": "0.0",
                        "stempmob": "0.0",
                        "structsoil": "0.001",
                        "structsoiln": "0.00004",
                        "structsoilp": "0.000002",
                        "structsurf": "0.001",
                        "structsurfn": "0.00004",
                        "structsurfp": "0.0000024",

                        # parameters
                        "resp_coeff": "0.2",      
                        "alpha_j": "0.308",       # Taking the theoretical maximum (from Belinda) 0.385 x 0.8 (leaf absorptance) = 0.308
                        "intercep_frac": "0.15",
                        "max_intercep_lai": "3.0",
                        "latitude": "-33.61",     # EucFACE parameter file
                        "albedo": "0.2",
                        "slamax": "6.34",    # EucFACE parameter file
                        "sla": "5.57",       # EucFACE parameter file
                        "slazero": "5.57",   # EucFACE parameter file
                        "lai_closed": "0.5", # I am effectively turning this feature off by setting it so low
                        "c_alloc_fmax": "0.52",  # EucFACE parameter file
                        "c_alloc_fmin": "0.44",  # EucFACE parameter file
                        "c_alloc_rmax": "0.32",  # EucFACE parameter file
                        "c_alloc_rmin": "0.22",  # EucFACE parameter file
                        "c_alloc_bmax": "0.1",   # guess
                        "c_alloc_bmin": "0.05",  # guess
                        "c_alloc_cmax": "0.0",   # turn off coarse roots!
                        "fretransn": "0.0",#31",     # EucFACE parameter file
                        "fretransp": "0.0",#53",    # EucFACE parameter file
                        "rretrans": "0.0",#3",      # EucFACE parameter file
                        "bretrans": "0.0",#7",      # EucFACE parameter file
                        "wretrans": "0.0",#82",     # EucFACE parameter file
                        "retransmob": "0.0",#82",     # EucFACE parameter file
                        "cretrans": "0.0",
                        "ncwnewz": "0.003",          #New stem ring N:C at zero leaf N:C (mobile)
                        "ncwnew": "0.003",           #New stem ring N:C at critical leaf N:C (mob)
                        "ncwimmz": "0.003",          #Immobile stem N C at zero leaf N C
                        "ncwimm": "0.003",           #Immobile stem N C at critical leaf N C
                        "ncbnewz": "0.003",          #new branch N C at zero leaf N C
                        "ncbnew": "0.003",           #new branch N C at critical leaf N C
                        "nccnewz": "0.003",          #new coarse root N C at zero leaf N C
                        "nccnew": "0.003",           #new coarse root N C at critical leaf N C
                        "ncrfac": "0.8",
                        "ncmaxfyoung": "0.03",
                        "ncmaxfold": "0.03",
                        "ncmaxr": "0.018",
                        "retransmob": "0.0",
                        "fdecay": "0.83",            # foliage decay rate, guess parameter
                        "fdecaydry": "0.83",         # foliage decay rate, guess parameter
                        "rdecay": "0.6",             # EucFACE parameter file
                        "rdecaydry": "0.6",          # as above
                        "crdecay": "0.00",           # turn off coarse roots!
                        "bdecay": "0.1",            # no idea, assuming 25 years
                        "wdecay": "0.1",            # no idea, assuming 25 years
                        "watdecaydry": "0.0",
                        "watdecaywet": "0.1",
                        "ligshoot": "0.15",          # EucFACE parameter file
                        "ligroot": "0.2",            # Same as Medlyn 2016
                        "rateuptake": "1.8",
                        "rateloss": "0.05",           # guess value
                        "topsoil_depth": "450.0",    # Not needed as I have supplied the root zone water and topsoil water available
                        "rooting_depth": "2500.0",   # Not needed as I have supplied the root zone water and topsoil water available
                        "wcapac_root": "300.0",      # [mm] (FC-WP)*rooting_depth. But using 2.0 m, site_description_stuff/EucFACE_Plot_Summary.doc
                        "wcapac_topsoil": "67.5",    # [mm] (FC-WP)*rooting_depth. But using 0.45 m, site_description_stuff/EucFACE_Plot_Summary.doc
                        "ctheta_topsoil": "0.65",     # Derive based on soil type loamy_sand
                        "ntheta_topsoil": "8.0",     # Derive based on soil type loamy_sand
                        "ctheta_root": "0.525",      # Derive based on soil type sandy_clay_loam
                        "ntheta_root": "5.5",        # Derive based on soil type sandy_clay_loam
                        "topsoil_type": "loamy_sand",
                        "rootsoil_type": "sandy_clay_loam",
                        "kp": "0.3",
                        "krp": "0.00001",
                        "dz0v_dh": "0.05",         # Using Value from JULES for TREE PFTs as I don't know what is best. However I have used value from Jarvis, quoted in Jones 1992, pg. 67. Produces a value within the bounds of 3.5-1.1 mol m-2 s-1 Drake, 2010, GCB for canht=17
                        "displace_ratio": "0.75",  # From Jones, pg 67, following Jarvis et al. 1976
                        "z0h_z0m": "1.0",
                        # root exudation
                        "a0rhizo": "0.05",
                        "a1rhizo": "0.6",

                        "g1": "3.04",            # EucFACE parameter
                        "jmaxna": "49.930",      # forcing intercept to zero; if use all species df, 49.743
                        "jmaxpa": "933.90",      # forcing intercept to zero; if use all species df, 842.46 
                        "jmaxnb": "0.0",         # forcing intercept to zero
                        "jmaxpb": "0.0",         # forcing intercept to zero
                        "vcmaxna": "27.707",     # forcing intercept to zero; if use all species df, 27.627
                        "vcmaxpa": "516.83",     # forcing intercept to zero; if use all species df, 468.76
                        "vcmaxnb": "0.0",        # forcing intercept to zero
                        "vcmaxpb": "0.0",        # forcing intercept to zero
                        "jmax": "162.91",         # EucFACE parameter file
                        "vcmax": "92.85",         # EucFACE parameter file
                        "measurement_temp": "25.0", # parameters obtained at 22 not 25 degrees
                        "heighto": "4.826",
                        "htpower": "0.35",
                        "height0": "5.0",
                        "height1": "25.0",
                        "leafsap0": "4000.0",     # "4000.0",
                        "leafsap1": "2700.0",     # 2700
                        "branch0": "5.61",
                        "branch1": "0.346",
                        "croot0": "0.34",
                        "croot1": "0.84",
                        "targ_sens": "0.5",
                        "density": "492.0",                   # EucFACE parameter file
                        "nf_min": "0.005", 
                        "sapturnover": "0.1",                 # guess value for EucFACE
                        "p_atm_deposition": "0.0",            # read in from met data now. 
                        "p_rate_par_weather": "0.00001",        # 
                        "p_rate_release_fertilizer": "1.0",   # 10 - 15 month release rate for slow-release fertilizer
                        "rate_sorb_ssorb": "0.01",           # fitted value, reasonable for EucFACE
                        "rate_ssorb_occ": "0.048",            # fitted value, reasonable for EucFACE
                        # sorption calculation
                        "smax": "0.01",                       # convert to unit in t ha-1, Yang et al., 2016, GRL, Table S2 
                        "ks": "0.006",                        # convert to unit in t ha-1, Yang et al., 2016, GRL, Table S2 
                        # biochemical P mineralization
                        "biochemical_p_constant": "150.0",
                        "max_p_biochemical": "0.001",
                        "crit_n_cost_of_p": "15.0",
                        "actpcmin": "0.01",         # empirical
                        "actpcmax": "0.02",         # empirical
                        "slowpcmin": "0.005",       # empirical
                        "slowpcmax": "0.011111",    # empirical
                        "passpcmin": "0.0051",      # empirical
                        "passpcmax": "0.0051",      # empirical
                        "pcbnew": "0.000286",       # same as sapwood
                        "pcbnewz": "0.000286",      # same as sapwood
                        "pccnew": "0.000286",       # same as sapwood
                        "pccnewz": "0.000286",      # same as sapwood
                        "pcmaxfold": "0.0014",      # EucFACE parameter file
                        "pcmaxfyoung": "0.002",     # EucFACE parameter file
                        "pcmaxr": "0.0006",         # EucFACE parameter file
                        "pcrfac": "0.8",
                        "pcwimm": "0.00013",        # EucFACE parameter file
                        "pcwimmz": "0.00013",       # EucFACE parameter file
                        "pcwnew": "0.000286",       # EucFACE parameter file
                        "pcwnewz": "0.000286",      # EucFACE parameter file
                        "pf_min": "0.0002",
                        # for calculate_p_ssorb_to_sorb function
                        "finesoil": "0.2",         # EucFACE parameter file
                        "phmin": "5.0",
                        "phmax": "14.0",  
                        "soilph": "5.52",          # EucFACE parameter file
                        "phtextmin": "0.000008",
                        "phtextmax": "0.00015",
                        "phtextslope": "0.00004",
                        "psecmnp": "0.000022",
                        # to determine soil SOM PC ratios
                        "pmin0": "0.0",            # set to zero for now
                        "pmincrit": "2.0",
                        "prateloss": "0.05",       # set it to be the same as N rate loss
                        "prateuptake": "0.9",      # Fitted value to obtain balance between uptake N:P ratio and reasonable P labile pool
                        "structcp": "5500.0",
                        "structratp": "0.0",

                        # control
                        "adjust_rtslow": "false",  # priming, off
                        "alloc_model": "fixed",
                        "assim_model": "mate",
                        "calc_sw_params": "true",   #false=use fwp values, true=derive them
                        "deciduous_model": "false",
                        "disturbance": "false",
                        "exudation": "false",
                        "fixed_stem_nc": "true",
                        "fixed_stem_pc": "true",
                        "fixleafnc": "false",
                        "fixleafpc": "false",
                        "grazing": "false",
                        "gs_model": "medlyn",
                        "aci_relationship": "walker",
                        "model_optroot": "false",
                        "modeljm": "1",
                        "ncycle": "true",
                        "pcycle": "true",
                        "nuptake_model": "1",
                        "puptake_model": "1",
                        "triose_p": "false",
                        "output_ascii": "true",
                        "passiveconst": "false",
                        "print_options": "end",
                        "ps_pathway": "c3",
                        "respiration_model": "fixed",
                        "strfloat": "0",
                        "strpfloat": "0",
                        "sw_stress_model": "1",  # Sands and Landsberg
                        "use_eff_nc": "0",
                        "text_effect_p": "1",
                        "water_stress": "true",

        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY_SPIN + cfg_fname)
        
        
    if POST_INDUST == True:

        # copy spunup base files to make two new experiment files
        shutil.copy(os.path.join(param_dir, "%s_%s_model_spunup.cfg" % (experiment_id, site)),
                    os.path.join(param_dir, "%s_%s_model_spunup_adj.cfg" % (experiment_id, site)))

        itag = "%s_%s_model_spunup_adj" % (experiment_id, site)
        otag = "%s_%s_model_indust" % (experiment_id, site)
        mtag = "%s_met_historic_daily_1750_2011.csv" % (site)
        out_fn = "%s_amb_equilib.csv" % (site)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)

        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
     
                         # control
                         "print_options": "end",
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)

    if POST_INDUST == True:

        # copy spunup base files to make two new experiment files
        shutil.copy(os.path.join(param_dir, "%s_%s_model_spunup.cfg" % (experiment_id, site)),
                    os.path.join(param_dir, "%s_%s_model_spunup_adj.cfg" % (experiment_id, site)))

        itag = "%s_%s_model_spunup_adj" % (experiment_id, site)
        otag = "%s_%s_model_indust" % (experiment_id, site)
        mtag = "%s_met_historic_daily_1750_2011.csv" % (site)
        out_fn = "%s_amb_equilib.csv" % (site)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)

        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
     
                         # control
                         "print_options": "daily",
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)
        
        # translate output to EucFACE requested output
        #tr.translate_output(out_fname, met_fname)
    
    # observed (2011-2019) under dry condition: store output 
    if OBS_DRY == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_indust.cfg" % (experiment_id, site)),
                    os.path.join(param_dir, "%s_%s_model_indust_adj.cfg" % (experiment_id, site)))

        itag = "%s_%s_model_indust_adj" % (experiment_id, site)
        otag = "%s_%s_model_DRY_%s_2012_2019" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_DRY_%s_daily_2012_2019.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_DRY_%s_2012_2019.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "daily",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)
        
        # translate output to EucFACE requested output
        #tr.translate_output(out_fname, met_fname)
    
    # observed (2011-2019) under dry condition: store cfg
    if OBS_DRY == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_indust.cfg" % (experiment_id, site)),
                    os.path.join(param_dir, "%s_%s_model_indust_adj.cfg" % (experiment_id, site)))

        itag = "%s_%s_model_indust_adj" % (experiment_id, site)
        otag = "%s_%s_model_DRY_%s_2012_2019" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_DRY_%s_daily_2012_2019.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_DRY_%s_2012_2019.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "end",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)
        
    # predicted (2020-2069) under DRY and NOP condition: store output 
    if PRD_DRY_NOP == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_DRY_%s_2012_2019.cfg" % (experiment_id, site, CO2_treatment)),
                    os.path.join(param_dir, "%s_%s_model_DRY_%s_2012_2019_adj.cfg" % (experiment_id, site, CO2_treatment)))

        itag = "%s_%s_model_DRY_%s_2012_2019_adj" % (experiment_id, site, CO2_treatment)
        otag = "%s_%s_model_DRY_%s_NOP_2020_2069" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_DRY_%s_NOP_daily_2020_2069.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_DRY_%s_NOP_2020_2069.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "daily",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)
        
        # translate output to EucFACE requested output
        #tr.translate_output(out_fname, met_fname)
    
    # predicted (2020-2069) under DRY and NOP condition: store cfg
    if PRD_DRY_NOP == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_DRY_%s_2012_2019.cfg" % (experiment_id, site, CO2_treatment)),
                    os.path.join(param_dir, "%s_%s_model_DRY_%s_2012_2019_adj.cfg" % (experiment_id, site, CO2_treatment)))

        itag = "%s_%s_model_DRY_%s_2012_2019_adj" % (experiment_id, site, CO2_treatment)
        otag = "%s_%s_model_DRY_%s_NOP_2020_2069" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_DRY_%s_NOP_daily_2020_2069.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_DRY_%s_NOP_2020_2069.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "end",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)   
        
    # predicted (2020-2069) under DRY and MDP condition: store output 
    if PRD_DRY_MDP == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_DRY_%s_2012_2019.cfg" % (experiment_id, site, CO2_treatment)),
                    os.path.join(param_dir, "%s_%s_model_DRY_%s_2012_2019_adj.cfg" % (experiment_id, site, CO2_treatment)))

        itag = "%s_%s_model_DRY_%s_2012_2019_adj" % (experiment_id, site, CO2_treatment)
        otag = "%s_%s_model_DRY_%s_MDP_2020_2069" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_DRY_%s_MDP_daily_2020_2069.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_DRY_%s_MDP_2020_2069.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "daily",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)
        
        # translate output to EucFACE requested output
        #tr.translate_output(out_fname, met_fname)
    
    # predicted (2020-2069) under DRY and MDP condition: store cfg
    if PRD_DRY_MDP == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_DRY_%s_2012_2019.cfg" % (experiment_id, site, CO2_treatment)),
                    os.path.join(param_dir, "%s_%s_model_DRY_%s_2012_2019_adj.cfg" % (experiment_id, site, CO2_treatment)))

        itag = "%s_%s_model_DRY_%s_2012_2019_adj" % (experiment_id, site, CO2_treatment)
        otag = "%s_%s_model_DRY_%s_MDP_2020_2069" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_DRY_%s_MDP_daily_2020_2069.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_DRY_%s_MDP_2020_2069.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "end",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)   
        
     # predicted (2020-2069) under DRY and HIP condition: store output 
    if PRD_DRY_HIP == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_DRY_%s_2012_2019.cfg" % (experiment_id, site, CO2_treatment)),
                    os.path.join(param_dir, "%s_%s_model_DRY_%s_2012_2019_adj.cfg" % (experiment_id, site, CO2_treatment)))

        itag = "%s_%s_model_DRY_%s_2012_2019_adj" % (experiment_id, site, CO2_treatment)
        otag = "%s_%s_model_DRY_%s_HIP_2020_2069" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_DRY_%s_HIP_daily_2020_2069.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_DRY_%s_HIP_2020_2069.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "daily",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)
        
        # translate output to EucFACE requested output
        #tr.translate_output(out_fname, met_fname)
    
    # predicted (2020-2069) under DRY and MDP condition: store cfg
    if PRD_DRY_HIP == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_DRY_%s_2012_2019.cfg" % (experiment_id, site, CO2_treatment)),
                    os.path.join(param_dir, "%s_%s_model_DRY_%s_2012_2019_adj.cfg" % (experiment_id, site, CO2_treatment)))

        itag = "%s_%s_model_DRY_%s_2012_2019_adj" % (experiment_id, site, CO2_treatment)
        otag = "%s_%s_model_DRY_%s_HIP_2020_2069" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_DRY_%s_HIP_daily_2020_2069.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_DRY_%s_HIP_2020_2069.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "end",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)   

    # observed (2012-2019) under wet condition: store output 
    if OBS_WET == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_indust.cfg" % (experiment_id, site)),
                    os.path.join(param_dir, "%s_%s_model_indust_adj.cfg" % (experiment_id, site)))

        itag = "%s_%s_model_indust_adj" % (experiment_id, site)
        otag = "%s_%s_model_WET_%s_2012_2019" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_WET_%s_daily_2012_2019.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_WET_%s_2012_2019.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "daily",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)
        
        # translate output to EucFACE requested output
        #tr.translate_output(out_fname, met_fname)
    
    # observed (2012-2019) under wet condition: store cfg
    if OBS_WET == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_indust.cfg" % (experiment_id, site)),
                    os.path.join(param_dir, "%s_%s_model_indust_adj.cfg" % (experiment_id, site)))

        itag = "%s_%s_model_indust_adj" % (experiment_id, site)
        otag = "%s_%s_model_WET_%s_2012_2019" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_WET_%s_daily_2012_2019.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_WET_%s_2012_2019.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "end",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)
        
    # predicted (2020-2069) under WET and NOP condition: store output 
    if PRD_WET_NOP == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_WET_%s_2012_2019.cfg" % (experiment_id, site, CO2_treatment)),
                    os.path.join(param_dir, "%s_%s_model_WET_%s_2012_2019_adj.cfg" % (experiment_id, site, CO2_treatment)))

        itag = "%s_%s_model_WET_%s_2012_2019_adj" % (experiment_id, site, CO2_treatment)
        otag = "%s_%s_model_WET_%s_NOP_2020_2069" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_WET_%s_NOP_daily_2020_2069.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_WET_%s_NOP_2020_2069.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "daily",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)
        
        # translate output to EucFACE requested output
        #tr.translate_output(out_fname, met_fname)
    
    # predicted (2020-2069) under WET and NOP condition: store cfg
    if PRD_WET_NOP == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_WET_%s_2012_2019.cfg" % (experiment_id, site, CO2_treatment)),
                    os.path.join(param_dir, "%s_%s_model_WET_%s_2012_2019_adj.cfg" % (experiment_id, site, CO2_treatment)))

        itag = "%s_%s_model_WET_%s_2012_2019_adj" % (experiment_id, site, CO2_treatment)
        otag = "%s_%s_model_WET_%s_NOP_2020_2069" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_WET_%s_NOP_daily_2020_2069.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_WET_%s_NOP_2020_2069.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "end",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)   
        
    # predicted (2020-2069) under WET and MDP condition: store output 
    if PRD_WET_MDP == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_WET_%s_2012_2019.cfg" % (experiment_id, site, CO2_treatment)),
                    os.path.join(param_dir, "%s_%s_model_WET_%s_2012_2019_adj.cfg" % (experiment_id, site, CO2_treatment)))

        itag = "%s_%s_model_WET_%s_2012_2019_adj" % (experiment_id, site, CO2_treatment)
        otag = "%s_%s_model_WET_%s_MDP_2020_2069" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_WET_%s_MDP_daily_2020_2069.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_WET_%s_MDP_2020_2069.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "daily",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)
        
        # translate output to EucFACE requested output
        #tr.translate_output(out_fname, met_fname)
    
    # predicted (2020-2069) under WET and MDP condition: store cfg
    if PRD_WET_MDP == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_WET_%s_2012_2019.cfg" % (experiment_id, site, CO2_treatment)),
                    os.path.join(param_dir, "%s_%s_model_WET_%s_2012_2019_adj.cfg" % (experiment_id, site, CO2_treatment)))

        itag = "%s_%s_model_WET_%s_2012_2019_adj" % (experiment_id, site, CO2_treatment)
        otag = "%s_%s_model_WET_%s_MDP_2020_2069" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_WET_%s_MDP_daily_2020_2069.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_WET_%s_MDP_2020_2069.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "end",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)   
        
    # predicted (2020-2069) under WET and HIP condition: store output 
    if PRD_WET_HIP == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_WET_%s_2012_2019.cfg" % (experiment_id, site, CO2_treatment)),
                    os.path.join(param_dir, "%s_%s_model_WET_%s_2012_2019_adj.cfg" % (experiment_id, site, CO2_treatment)))

        itag = "%s_%s_model_WET_%s_2012_2019_adj" % (experiment_id, site, CO2_treatment)
        otag = "%s_%s_model_WET_%s_HIP_2020_2069" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_WET_%s_HIP_daily_2020_2069.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_WET_%s_HIP_2020_2069.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "daily",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)
        
        # translate output to EucFACE requested output
        #tr.translate_output(out_fname, met_fname)
    
    # predicted (2020-2069) under WET and MDP condition: store cfg
    if PRD_WET_HIP == True:
        
        # copy last cfg file and make new one
        shutil.copy(os.path.join(param_dir, "%s_%s_model_WET_%s_2012_2019.cfg" % (experiment_id, site, CO2_treatment)),
                    os.path.join(param_dir, "%s_%s_model_WET_%s_2012_2019_adj.cfg" % (experiment_id, site, CO2_treatment)))

        itag = "%s_%s_model_WET_%s_2012_2019_adj" % (experiment_id, site, CO2_treatment)
        otag = "%s_%s_model_WET_%s_HIP_2020_2069" % (experiment_id, site, CO2_treatment)
        mtag = "%s_met_WET_%s_HIP_daily_2020_2069.csv" % (site, CO2_treatment)
        out_fn = "%s_simulated_WET_%s_HIP_2020_2069.csv" % (site, CO2_treatment)
        out_param_fname = os.path.join(param_dir, otag + ".cfg")
        cfg_fname = os.path.join(param_dir, itag + ".cfg")
        met_fname = os.path.join(met_dir, mtag)
        out_fname = os.path.join(run_dir, out_fn)
        replace_dict = {
                         # files
                         "out_param_fname": "%s" % (out_param_fname),
                         "cfg_fname": "%s" % (cfg_fname),
                         "met_fname": "%s" % (met_fname),
                         "out_fname": "%s" % (out_fname),
    
                         # control
                         "print_options": "end",
    
                        }
        ad.adjust_param_file(cfg_fname, replace_dict)
        os.system(GDAY + cfg_fname)  



if __name__ == "__main__":

    experiment_id = "FACE"
    site = "EUC"
    CO2_treatment = "AMB"
    
    main(experiment_id, site, SPIN_UP=True, POST_INDUST=True, 
    OBS_DRY=True, OBS_WET=True, PRD_DRY_NOP=True, PRD_WET_NOP=True,
    PRD_DRY_MDP=True, PRD_WET_MDP=True, PRD_DRY_HIP=True, PRD_WET_HIP=True)
    #
    #main(experiment_id, site, SPIN_UP=True, POST_INDUST=True, 
    #OBS_DRY=True, OBS_WET=False, PRD_DRY_NOP=True, PRD_WET_NOP=False,
    #PRD_DRY_MDP=False, PRD_WET_MDP=False, PRD_DRY_HIP=False, PRD_WET_HIP=False)
