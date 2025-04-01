""" 
Automatically start simulations and loop over parameters

"""

import os
import numpy as np
import pandas as pd
import random
import concurrent.futures
from itertools import repeat

# IO Parameters

smoldyn_path = "/usr/local/bin/smoldyn" #exe file to call smoldyn
root_path = "/Users/koerbem/Documents/smoldyn/" 
model_file = "model_2_5Phos_kcd45_klck.txt"
results_folder = "/Users/koerbem/Documents/smoldyn/Model_9_5Phos_scan-k/"
simnum = 50

# Model Parameters
pMHC_mdens=[0,10]                       # mol/um^2 
cd45_exclusion=[0,20,40,60,80,100]      # percent of exclusion inside contact, i.e 1-cd45_exclusion inside contact
k_cd45 = [0.05, 0.01, 0.005, 0.001]     # s-1 um-2
k_lck = [0.05, 0.01, 0.005, 0.001]      # s-1 um-2
randomint = 42                          # answer to everything

# Function Definitions

def simulate(smoldyn_path, root_path, model_file, results_folder, simnum, pMHC_dens, cd45_exclusion, k_cd45, k_lck):
    os.chdir(root_path)
    string = smoldyn_path + ' ' + model_file + " --define python=42 --define results_folder="+results_folder+" --define sim_num={:.0f} --define randomint={:.0f} --define pMHC_mdens={:.0f} --define cd45_exclusion={:.0f} --define k_cd45={:.5f} --define k_lck={:.5f} -tqw".format(simnum, random.randint(1,1e10), pMHC_dens, cd45_exclusion, k_cd45, k_lck)
    print(string)
    os.system(string)
    return 0

# Main

if __name__ == '__main__':
    # Run simulations
    random.seed()
    print("Chosen model: "+root_path+model_file)
    print("Starting Simulation")

    for i_cd45 in k_cd45:
        for i_lck in k_lck:
            for i_exclusion in cd45_exclusion:
                for i_mhc in pMHC_mdens:
                    print('Simulating k_CD45 {} and k_LCK {}'.format(i_cd45, i_lck))
                    with concurrent.futures.ProcessPoolExecutor() as executor:
                        results = executor.map(simulate, repeat(smoldyn_path), repeat(root_path), repeat(model_file), repeat(results_folder), range(1, simnum+1), repeat(i_mhc), repeat(i_exclusion), repeat(i_cd45), repeat(i_lck))
                

    print("\r Well Done!")