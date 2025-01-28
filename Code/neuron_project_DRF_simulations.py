# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 11:54:40 2021

@author: brndn
"""

from neuron import h 
import numpy as np
import pandas as pd
from neuron_sim_functions import *

# get combined control values from BAF and DRF experiments
def getControlValues_BAFCPZ():
    ih_BAF_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project 2/"
    ih_BAF_filename = "Ih_model_values_BAF312.csv"
    ih_CPZ_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project 2/"
    ih_CPZ_filename = "Ih_model_values_DRF.csv"
    
    ih_BAF_df = pd.read_csv(ih_BAF_path + ih_BAF_filename)
    ih_CPZ_df = pd.read_csv(ih_CPZ_path + ih_CPZ_filename)
    
    BAF_indices = [i for i, elem in enumerate(ih_BAF_df['Name']) if "control_" in elem]
    CPZ_indices = [i for i, elem in enumerate(ih_CPZ_df['Name']) if "control_" in elem]

    features = ["current_density","v_half","slope_factor","tau_fast"]
    controlDict = {}   
    for feature in features:
        control_values = np.concatenate((ih_BAF_df[feature][BAF_indices].values,ih_CPZ_df[feature][CPZ_indices].values),axis=None)       
        controlDict[feature] = np.mean(control_values)
    return controlDict
    
ih_param_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project 2/"
ih_param_filename = "Ih_model_values_DRF.csv"
# CPZ1DRF_3 doesn't fire any spikes with kleak = 0.004, current density too low

worksheet_path = ih_param_path
worksheet_filename = "DRF_TC0_Voltage_Simulations_kleak0.004_combinedcontrols"

cell_names = ["TC[0]","TC[1]","RE[0]","RE[1]"]

norm_ghbar = 1e-5 # default max Ih conductance (S/cm^2)
tStop = 5000 # stop time in ms

h.load_file("stdrun.hoc") # for run control
h.nrn_load_dll("C:/Users/brndn/OneDrive/Desktop/DLGN_NEW/nrn_mech.dll")
h.load_file('C:/Users/brndn/OneDrive/Desktop/DLGN_NEW/Fdelta.oc')
h.tstop = tStop

trec,vrec_tc0,vrec_tc1,vrec_re0,vrec_re1 = h.Vector(),h.Vector(),h.Vector(),h.Vector(),h.Vector()

trec.record(h._ref_t) # record time
vrec_tc0.record(h.TC[0].soma[0](0.5)._ref_v) # record voltage from center
vrec_tc1.record(h.TC[1].soma[0](0.5)._ref_v) # record voltage from center
vrec_re0.record(h.RE[0].soma[0](0.5)._ref_v) # record voltage from center
vrec_re1.record(h.RE[1].soma[0](0.5)._ref_v) # record voltage from center

v_tc0_df = pd.DataFrame() # initialize dataframe for TC[0] voltage traces

ih_param_df = pd.read_csv(ih_param_path + ih_param_filename)

#cd_control,indices = getControlCD_DRF(ih_param_df)
controlDict = getControlValues_BAFCPZ()
cd_control = controlDict["current_density"]

print("cd_control = %.4g (pA/pF)" % cd_control)
#print("")
#print(indices)

for i in range(len(ih_param_df['Name'])):
    
    save_filename = ih_param_df['Name'].iloc[i]   
    print("")
    print(save_filename)
    
   # set new ih params for both TC cells
    h.shift_iar =  75 + ih_param_df['v_half'].iloc[i]
    h.slope_iar = ih_param_df['slope_factor'].iloc[i]
    h.tau130_iar = ih_param_df['tau_fast'].iloc[i]
    
    for iCell in range(int(h.ncells)):    
        h.TC[iCell].kl.gmax = 0.004
        cd = ih_param_df['current_density'].iloc[i] 
        if iCell == 0: print("cd = %.4g (pA/pF)" % cd)
        h.TC[iCell].soma[0].ghbar_iar = norm_ghbar * cd / cd_control
        
    h.frecord_init()
    h.run()
    
    time,v_tc0,v_tc1,v_re0,v_re1 = vecToNumpy(trec,vrec_tc0,vrec_tc1,vrec_re0,vrec_re1)
    
    data_matrix = np.concatenate((time,v_tc0,v_tc1,v_re0,v_re1)).reshape(5,trec.size()).transpose()
    
    inter_IBI_dict,inter_freq_dict,nAP_dict,mean_nAP_dict,time_AP1_dict,intra_IBI_dict,intra_freq_dict = makeIBIdict(data_matrix,cell_names,printnum=1)
    
    plotSim(data_matrix,cell_names,inter_freq_dict,time_AP1_dict,intra_freq_dict,mean_nAP_dict,save_filename)
    
    saveFig(save_filename)
    
    saveData(save_filename,[data_matrix,inter_IBI_dict,inter_freq_dict,nAP_dict,mean_nAP_dict,time_AP1_dict,intra_IBI_dict,intra_freq_dict])
    
    if i == 0: v_tc0_df["time_ms"] = time
    
    v_tc0_df[save_filename] = v_tc0
    
v_tc0_df.to_excel("%s%s.xlsx" % (worksheet_path,worksheet_filename))
