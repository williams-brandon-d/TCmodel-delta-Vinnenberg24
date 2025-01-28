# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 12:22:33 2021

@author: brndn
"""
# imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy import stats
import math
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import random


# functions
def loadData(save_filename):
    data_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project/Data/"
    data_filename = "%s.pckl" % save_filename   
    f = open(data_path+data_filename, 'rb')
    [_,_,inter_freq_dict,_,mean_nAP_dict,time_AP1_dict,_,intra_freq_dict] = pickle.load(f)
    f.close()
    return inter_freq_dict,mean_nAP_dict,time_AP1_dict,intra_freq_dict

def plotBar(group_names,mean_list,xlabel_string,sem_list,dataDict,group_labels,title_string):
    
    rev_label_list = group_labels.copy()
    rev_label_list.reverse()
    rev_mean_list = mean_list.copy()
    rev_mean_list.reverse()
    rev_sem_list = sem_list.copy()
    rev_sem_list.reverse()
    
    colors = ['red','orange','blue','cyan','green','olive','purple','pink']
    color_list = colors[:len(group_names)]
    rev_color_list = color_list.copy()
    rev_color_list.reverse()
    rand_range = 0.1
    
    x_pos = [i for i, _ in enumerate(rev_label_list)]
    
    plt.figure()
    plt.barh(x_pos, rev_mean_list,color=color_list, xerr=rev_sem_list,alpha=0.5)
    plt.xlabel(xlabel_string,fontweight="bold")
    plt.yticks(x_pos, rev_label_list,fontweight="bold")
    plt.xticks(fontweight="bold")
    plt.suptitle("P-Values (Independent sample comparison to control)")
    plt.title(title_string)
    for iGroup in range(len(group_names)):
        group = group_names[iGroup]
        rev_idx = -(iGroup + 1)
        rev_pos = x_pos[rev_idx]
        #s = np.random.normal(0, rand_range/3, len(dataDict[group]))
        x_list = [random.uniform(rev_pos-rand_range,rev_pos+rand_range) for i in range(0, len(dataDict[group]))]
        #x_list = rev_pos*np.ones(len(dataDict[group]))
        plt.scatter(dataDict[group],x_list,color=rev_color_list[iGroup],alpha=0.6)
    plt.show()

def findUniqueGroups(ih_param_df):
    group_name_list = [] # find unique group names
    for iCell in range(len(ih_param_df['Name'])):
        cell_name = ih_param_df['Name'].iloc[iCell]
        underscore_index = cell_name.find('_')
        group_name = cell_name[:underscore_index] # remove cell number from name
#        if  group_name.find('control') > 0:
#            group_name = group_name[:-len('control')] # remove control from name
        group_name_list.append(group_name) 
    unique_names = np.unique(np.array(group_name_list)).tolist()
    unique_names.sort() # sorts normally by alphabetical order
    unique_names.sort(key=str.lower) # sorts by lowercase first
    #unique_names.sort(key=len) # sorts by descending length
    return unique_names

def findNamesInGroup(names,group):
    new_names = [] # find cell names in group
    for iCell in range(len(names)):
        cell_name = names[iCell]
        underscore_index = cell_name.find('_')
        group_name = cell_name[:underscore_index] # remove cell number from name
        new_names.append(group_name) 

    #new_names = [elem[ : -2] for elem in names] # remove last 2 characters
    indices = [i for i, elem in enumerate(new_names) if group == elem] # compare to group
    return [names[i] for i in indices]

def round_to_n(x,n): # round x to n significant digits
    if x == 0: 
        return x 
    else: 
        return round(x, -int( math.floor( math.log10( abs(x) ) ) ) + (n - 1) )
    
def tukey_hsd(group_names , *args ):
    endog = np.hstack(args)
    groups_list = []
    for i in range(len(args)):
        for j in range(len(args[i])):
            groups_list.append(group_names[i])
    groups = np.array(groups_list)
    res = pairwise_tukeyhsd(endog, groups)
    #print (res.pvalues) #print only p-value
    print("")
    print(res) #print result
    
# choose ih parameter file to get names and data
ih_param_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project 2/"
ih_param_filename = "Ih_model_values_DRF.csv"

# choose data groups to analyze and plot
# group_names = ["CPZday1control","CPZday1","CPZday7control","CPZday7","CPZday25control","CPZday25"]
# group_names = ["bl6control","BL6_IFNalpha","BL6_IL1beta"]

# choose a path and filename to save data in an excel spreadsheet
save_path = ih_param_path
save_filename = "DRF Bar Chart Data"

# choose data features to analyze and provide labels for graphs and spreadsheet names
features = ["v_half","current_density","slope_factor","tau_fast"]
feature_labels = [r"$\bf{V_{0.5}}$ (mV)","Current Density (pA/pF)","Slope Factor",r"$\bf{Tau_{fast}}$ (ms)"]

#data_key = "TC[0]" # choose cell data to load from each simulation

sigDigits = 3 # choose significant digits for displaying P-values

alphaValue = 0.05 # choose alpha value for comparison tests

# load parameter data from .csv file
ih_param_df = pd.read_csv(ih_param_path+ih_param_filename)

# find unique group names
#group_names = findUniqueGroups(ih_param_df)
group_names = ['control','CPZ1','CPZ7','CPZ25'];

# intialize dicts and lists
dataDict = {}
for group in group_names: 
    dataDict[group] = {}
    for feature in features: dataDict[group][feature] = [] 

meanDict = {}
for feature in features: meanDict[feature] = []
    
semDict = {}
for feature in features: semDict[feature] = []
        
cellDict = {}

# calculate mean and S.E.M. for each group and feature
for group in group_names:    
    cell_names = findNamesInGroup(ih_param_df['Name'].values.tolist(),group)

    for cell in cell_names:
        #if cell == "CPZ1_3" or cell == "CPZ1DRF_3": continue # skip datapoints due to no spikes
        #cellDict[features[0]],cellDict[features[1]],cellDict[features[2]],cellDict[features[3]] = loadData(cell)
        for feature in features:
            dataDict[group][feature].append( ih_param_df.loc[ih_param_df['Name'] == cell][feature].values[0] )
            # remove empty values
            #dataDict[group][feature] = list(filter(None, dataDict[group][feature]))

    for feature in features:
        meanDict[feature].append( np.mean( dataDict[group][feature] ) )
        semDict[feature].append( stats.sem( dataDict[group][feature] ) )

# save bar graph values in excel spreadsheet
flipDict = {}
dfDict = {}
#writer = pd.ExcelWriter("%s.xlsx" % (save_path+save_filename))
for iFeature in range(len(features)):
    feature = features[iFeature]
    flipDict[feature] = {}
    for group in group_names:
        flipDict[feature][group] = dataDict[group][feature]
    dfDict[feature] = pd.DataFrame({k:pd.Series(v) for k,v in flipDict[feature].items()}) # pad dict with NaN for equals column lengths
#    dfDict[feature].to_excel(writer,sheet_name=feature_labels[iFeature],index=False)
#writer.save()
#writer.close()

# check normality of each group
normDict = {}
for feature in features:
    normDict[feature] = {}
    for group in group_names:
        #dataDict[group][feature] = list(filter(None, dataDict[group][feature]))
        _,normDict[feature][group] = stats.shapiro( dataDict[group][feature] ) # Pvalue > 0.05 = probably normal
        normDict[feature][group] = round_to_n( normDict[feature][group], sigDigits )

# test for statistical difference between control and test data for each day (CPZ data)
# group_names = findUniqueGroups(ih_param_df)
# pValueDict = {}
# for feature in features:
#     pValueDict[feature] = {}
#     for name in unique_names:
#         _,pValueDict[feature][name] = stats.ttest_ind(dataDict[name+"control"][feature],dataDict[name][feature])
#         pValueDict[feature][name] = round_to_n( pValueDict[feature][name], sigDigits )

# test for statistical difference between control every other data set
controlPvalueDict = {}
for feature in features:
    controlPvalueDict[feature] = {}
    for group in group_names:
        group1 = group_names[0]
        group2 = group;
        if group == group_names[0]: continue # skip control group
        if normDict[feature][group1] > alphaValue and normDict[feature][group2] > alphaValue: # paired ttest
            _,controlPvalueDict[feature][group] = stats.ttest_ind(dataDict[group1][feature],dataDict[group2][feature])
        else:
            _,controlPvalueDict[feature][group] = stats.ranksums(np.array(dataDict[group1][feature]),np.array(dataDict[group2][feature]))
        controlPvalueDict[feature][group] = round_to_n( controlPvalueDict[feature][group], sigDigits )
        print("%s T-test Control vs %s: p-value = %.3g" % (feature,group,controlPvalueDict[feature][group]))

# test for statistical difference between age-matched controls
# nGroups = len(group_names)
# nPairs = nGroups/2
# ageDRFpValueDict = {}
# for feature in features:
#     ageDRFpValueDict[feature] = {}
#     for iPair in range(int(nPairs)):
#         group1 = group_names[2*iPair]
#         group2 = group_names[2*iPair+1]
        
#         if normDict[feature][group1] > alphaValue and normDict[feature][group2] > alphaValue: # paired ttest
#             _,ageDRFpValueDict[feature][group1] = stats.ttest_rel(dataDict[group1][feature],dataDict[group2][feature])
#         else:
#             diff = np.array(dataDict[group1][feature]) - np.array(dataDict[group2][feature])
#             _,ageDRFpValueDict[feature][group1] = stats.wilcoxon( diff )
        
#         ageDRFpValueDict[feature][group1] = round_to_n( ageDRFpValueDict[feature][group1], sigDigits )
#     print("\n"+feature+" DRF P-values")
#     print(ageDRFpValueDict[feature])

# anova
anovaPvalueDict = {}
for feature in features:
    _,anovaPvalueDict[feature] = stats.f_oneway(dataDict[group_names[0]][feature],
                                                dataDict[group_names[1]][feature],
                                                dataDict[group_names[2]][feature])
    anovaPvalueDict[feature] = round_to_n(anovaPvalueDict[feature],sigDigits)
    print("\n"+feature)
    tukey_hsd(group_names,dataDict[group_names[0]][feature],
                          dataDict[group_names[1]][feature],
                          dataDict[group_names[2]][feature])

print("Anova P-values")
print(anovaPvalueDict)

# plot bar graphs for each feature
#group_labels = ["Control",r"$\mathbf{IFN_\alpha}$",r"$\mathbf{IL_\beta}$"]
group_labels = group_names
for iFeature in range(len(features)):
    feature = features[iFeature]
    plotBar(group_names,meanDict[feature],feature_labels[iFeature],semDict[feature],flipDict[feature],group_labels,controlPvalueDict[feature])

# compare across different days

# good online stats reference = https://help.xlstat.com/s/article/which-statistical-test-should-you-use?language=en_US

# normality tests
# shapiro-wilk test has greatest statistical power, stats.shapiro()

# equal variance tests
# Levene's test is the most robust to non-normality, stats.levene()
# if normally distributed Bartlett's test can be used, stats.bartlett()

# for comparing two independent samples
# if both normal:
#   if equal variance: ttest, stats.ttest_ind(data1,data2)
#   if not equal variance: Welch's ttest, stats.ttest_ind(data1,data2,equal_var=False)
# if not normal: Mann Whitney U test, stats.mannwhitneyu()

# for comparing two paired samples
# if normal: paired ttest, stats.ttest()
# if not normal: Wilcoxon's signed rank test, stats.wilcoxon()

# for comparing more than 2 independent samples
# if normal: one way ANOVA, stats.f_oneway()
#   for multiple comparisons test
#   if sample size is equal: TukeyHSD, statsmodels.stats.multicomp.pairwise_tukeyhsd
#   if sample size is not equal: Tukey-Kramer, 
# if not normal: Kruskal-Wallis H-test, stats.kruskal()
