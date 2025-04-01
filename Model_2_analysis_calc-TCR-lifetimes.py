# Set input folder, calculate from simulatin output the phosphorylated TCR lifetimes
# Based on simulation output that preserves the serial number of TCRs

import numpy as np
import pandas as pd
from pathlib import Path
from os import listdir
import concurrent.futures
from itertools import repeat

def get_activation_time(p, t=10):
    """ 
    Find the first occurrence of p where p>0 for t consecutive times
    """
    cnt = 0
    time = 0

    p = p>0
    for i, i_p in enumerate(p):
        if i_p and (cnt == 0):
            time = i
            cnt = 1
        elif i_p and (cnt > 0):
            cnt += 1
        elif not i_p:
            cnt = 0
        
        if cnt == t:
            return time
    
    return False

def get_activation_time_lifetime(track, t=10, serialns = [7]):
    """ 
    Find the first occurrence within track where serial_number in serialns for t consecutive times, and track the lifetimes of each state 
    If during the counting of a lifetime the data array ends, no lifetime is reported for this state occurence
    
    t is the minmal lifetime of the serialns
    """

    cnt = 0
    time = 0
    dwell_time = False
    trig_time = False
    trig_inside = False
    lifetimes = []
    dwell_times = []
    inside = []
    track.reset_index(drop=True, inplace=True)

    for i, i_p in track.iterrows():
        is_p = np.any([i_p['species_type']==i for i in serialns])
        if is_p and (cnt == 0):
            time = i
            cnt = 1            
        elif is_p and (cnt > 0):
            cnt += 1
        elif (not is_p) & (cnt > 0):
            lifetimes.append(cnt)
            inside.append(track.iloc[time]['inside'])
            cnt = 0
        
        if (cnt >= t) & (trig_time == False):
            trig_time = time + t 
            trig_inside = track.iloc[time]['inside']       

        if i_p['inside'] & (dwell_time == False):
            dwell_time = i_p['time']
        #elif i_p['inside'] & (dwell_cnt > 0):
        #    dwell_cnt += 1
        elif (not i_p['inside']) & (dwell_time != False):
            dwell_times.append(i_p['time']-dwell_time)
            dwell_time = False

    return trig_time, trig_inside, lifetimes, inside, dwell_times

def calculate_lifetime(file, folder, columnnames_tcrs, columnnames, r_sim, r_c, A_c, A_out, dt, min_lifetime):
    """ 
    min_lifetime in s
    """

    print('Working on file ' + file)
    # TCR equilibrium density
    _,model,CD45_den,pMHC_mdens,k_off,CD45_exclusion,k_CD45,k_lck,fileend= file.split('_')
    repeat, fileend = fileend.split('.')

    # Read population output and calculate equilibrium values
    dat = pd.read_csv(folder/file, header=None)
    dat.columns = columnnames    
    sim_time = np.max(dat['time']) # s
    tcrp3_equ_dens = np.mean(dat.loc[dat['time']>sim_time/2,'tcrp3'])/(r_sim**2 * np.pi)
    CD45_exclusion_measured = np.mean((1-dat.loc[dat['time']>sim_time/2,'cd45i']*A_out/(dat.loc[dat['time']>sim_time/2,'cd45o']*A_c))*100)
    
    # Read individual TCR tracks
    tracks = pd.read_csv(folder/file.replace('.txt', '_tcrs.txt'), header=None)
    tracks.columns = columnnames_tcrs
    tracks['timestep'] = (tracks['timestep']-1) * 10 + 1 # *10 here because output is only written every 10th timestep
    tracks['time'] = tracks['timestep'] * dt
    n_tcrs = np.max(tracks['sn']) # sn is within interval [1,n_tcrs]

    # TCR phosphorylated inside
    tracks['distance'] = np.sqrt(tracks['x']**2 + tracks['y']**2)
    tracks['inside'] = tracks['distance']<r_c
    
    # TCR phosphorylation lifetime 
    lifetimes3b_in = []
    lifetimes3b_out = []
    lifetimes3_in = []
    lifetimes3_out = []
    trigb_times = []
    trig_times = []
    trigb_inside = []
    trig_inside = []
    dwell_times = []
    for sni in range(n_tcrs):
        track = tracks.loc[tracks['sn']==(sni+1)]
        
        # bound tcrs 3P
        tcrp3b_trig_time, tcrp3b_trig_inside, tcrp3b_lifetimes, tcrp3b_inside, tcr_dwell_times = get_activation_time_lifetime(track, t=int(min_lifetime/(dt*10)), serialns = [8])
        trigb_times.append(tcrp3b_trig_time)
        trigb_inside.append(tcrp3b_trig_inside)
        [lifetimes3b_in.append(tcrp3b_lifetimes[i]) for i in range(len(tcrp3b_lifetimes)) if tcrp3b_inside[i]]
        [lifetimes3b_out.append(tcrp3b_lifetimes[i]) for i in range(len(tcrp3b_lifetimes)) if not tcrp3b_inside[i]]
        [dwell_times.append(i) for i in tcr_dwell_times]

        # bound and unbound tcrs 3P
        tcrp3_trig_time, tcrp3_trig_inside, tcrp3_lifetimes, tcrp3_inside, _ = get_activation_time_lifetime(track, t=int(min_lifetime/(dt*10)), serialns = [7,8])
        trig_times.append(tcrp3_trig_time)
        trig_inside.append(tcrp3_trig_inside)
        [lifetimes3_in.append(tcrp3_lifetimes[i]) for i in range(len(tcrp3_lifetimes)) if tcrp3_inside[i]]
        [lifetimes3_out.append(tcrp3_lifetimes[i]) for i in range(len(tcrp3_lifetimes)) if not tcrp3_inside[i]]

    # average tcrp3_cnt
    track_half2 = tracks.loc[(tracks['time']>=sim_time/2)]
    tcrp3_equ_in = np.sum(track_half2['inside'] & ( (track_half2['species_type']==7) | (track_half2['species_type']==8) ))/len(np.unique(track_half2['timestep']))
    tcrp3_equ_out = np.sum(np.logical_not(track_half2['inside']) & ( (track_half2['species_type']==7) | (track_half2['species_type']==8) ))/len(np.unique(track_half2['timestep']))

    results = {
        'model': eval(model),
        'CD45_den':  eval(CD45_den),
        'pMHC_mdens': eval(pMHC_mdens),
        'k_off': eval(k_off),
        'k_CD45': eval(k_CD45),
        'CD45_exclusion': eval(CD45_exclusion),
        'CD45_exclusion_measured': CD45_exclusion_measured,
        'k_lck': eval(k_lck),
        'repeat': eval(repeat),
        'A_c': A_c,
        'A_out': A_out,
        'tcrp3_equ_dens': tcrp3_equ_dens, # p3 but unbound
        'tcrp3_equ_in': tcrp3_equ_in, # only phopshorylated inside the contact
        'tcrp3_equ_out': tcrp3_equ_out, # phosphorylated outside the contact to control if exclusion independent
        'tcrp3_lifetimes_inside': [lifetimes3_in],
        'tcrp3_lifetimes_outside': [lifetimes3_out],
        'tcrp3b_lifetimes_inside': [lifetimes3b_in],
        'tcrp3b_lifetimes_outside': [lifetimes3b_out],
        'tcrp3_trig_times': [trig_times],
        'tcrp3b_trig_times': [trigb_times],
        'tcrp3_trig_inside': [trig_inside],
        'tcrp3b_trig_inside': [trigb_inside],
        'tcr_dwell_times': [dwell_times]
    }
    results = pd.DataFrame(results)
    results.to_json(folder/file.replace('.txt', '_tcr_lifetimes.json'))

def calculate_lifetime_5P(file, folder, columnnames_tcrs, columnnames, r_sim, r_c, A_c, A_out, dt, min_lifetime):
    """ 
    min_lifetime in s
    """

    print('Working on file ' + file)
    # TCR equilibrium density
    _,model,_,CD45_den,pMHC_mdens,k_off,CD45_exclusion,k_CD45,k_lck,fileend= file.split('_')
    repeat, fileend = fileend.split('.')

    # Read population output and calculate equilibrium values
    dat = pd.read_csv(folder/file, header=None)
    dat.columns = columnnames    
    sim_time = np.max(dat['time']) # s
    tcrp5_equ_dens = np.mean(dat.loc[dat['time']>sim_time/2,'tcrp5'])/(r_sim**2 * np.pi)
    CD45_exclusion_measured = np.mean((1-dat.loc[dat['time']>sim_time/2,'cd45i']*A_out/(dat.loc[dat['time']>sim_time/2,'cd45o']*A_c))*100)
    
    # Read individual TCR tracks
    tracks = pd.read_csv(folder/file.replace('.txt', '_tcrs.txt'), header=None)
    tracks.columns = columnnames_tcrs
    tracks['timestep'] = (tracks['timestep']-1) * 10 + 1 # *10 here because output is only written every 10th timestep
    tracks['time'] = tracks['timestep'] * dt
    n_tcrs = np.max(tracks['sn']) # sn is within interval [1,n_tcrs]

    # TCR phosphorylated inside
    tracks['distance'] = np.sqrt(tracks['x']**2 + tracks['y']**2)
    tracks['inside'] = tracks['distance']<r_c
    
    # TCR phosphorylation lifetime 
    lifetimes5b_in = []
    lifetimes5b_out = []
    lifetimes5_in = []
    lifetimes5_out = []
    trigb_times = []
    trig_times = []
    trigb_inside = []
    trig_inside = []
    dwell_times = []
    for sni in range(n_tcrs):
        track = tracks.loc[tracks['sn']==(sni+1)]
        
        # bound tcrs 3P
        tcrp5b_trig_time, tcrp5b_trig_inside, tcrp5b_lifetimes, tcrp5b_inside, tcr_dwell_times = get_activation_time_lifetime(track, t=int(min_lifetime/(dt*10)), serialns = [12])
        trigb_times.append(tcrp5b_trig_time)
        trigb_inside.append(tcrp5b_trig_inside)
        [lifetimes5b_in.append(tcrp5b_lifetimes[i]) for i in range(len(tcrp5b_lifetimes)) if tcrp5b_inside[i]]
        [lifetimes5b_out.append(tcrp5b_lifetimes[i]) for i in range(len(tcrp5b_lifetimes)) if not tcrp5b_inside[i]]
        [dwell_times.append(i) for i in tcr_dwell_times]

        # bound and unbound tcrs 3P
        tcrp5_trig_time, tcrp5_trig_inside, tcrp5_lifetimes, tcrp5_inside, _ = get_activation_time_lifetime(track, t=int(min_lifetime/(dt*10)), serialns = [11,12])
        trig_times.append(tcrp5_trig_time)
        trig_inside.append(tcrp5_trig_inside)
        [lifetimes5_in.append(tcrp5_lifetimes[i]) for i in range(len(tcrp5_lifetimes)) if tcrp5_inside[i]]
        [lifetimes5_out.append(tcrp5_lifetimes[i]) for i in range(len(tcrp5_lifetimes)) if not tcrp5_inside[i]]

    # average tcrp3_cnt
    track_half2 = tracks.loc[(tracks['time']>=sim_time/2)]
    tcrp5_equ_in = np.sum(track_half2['inside'] & ( (track_half2['species_type']==11) | (track_half2['species_type']==12) ))/len(np.unique(track_half2['timestep']))
    tcrp5_equ_out = np.sum(np.logical_not(track_half2['inside']) & ( (track_half2['species_type']==11) | (track_half2['species_type']==12) ))/len(np.unique(track_half2['timestep']))

    results = {
        'model': eval(model),
        'CD45_den':  eval(CD45_den),
        'pMHC_mdens': eval(pMHC_mdens),
        'k_off': eval(k_off),
        'k_CD45': eval(k_CD45),
        'CD45_exclusion': eval(CD45_exclusion),
        'CD45_exclusion_measured': CD45_exclusion_measured,
        'k_lck': eval(k_lck),
        'repeat': eval(repeat),
        'A_c': A_c,
        'A_out': A_out,
        'tcrp5_equ_dens': tcrp5_equ_dens, # p3 but unbound
        'tcrp5_equ_in': tcrp5_equ_in, # only phopshorylated inside the contact
        'tcrp5_equ_out': tcrp5_equ_out, # phosphorylated outside the contact to control if exclusion independent
        'tcrp5_lifetimes_inside': [lifetimes5_in],
        'tcrp5_lifetimes_outside': [lifetimes5_out],
        'tcrp5b_lifetimes_inside': [lifetimes5b_in],
        'tcrp5b_lifetimes_outside': [lifetimes5b_out],
        'tcrp5_trig_times': [trig_times],
        'tcrp5b_trig_times': [trigb_times],
        'tcrp5_trig_inside': [trig_inside],
        'tcrp5b_trig_inside': [trigb_inside],
        'tcr_dwell_times': [dwell_times]
    }
    results = pd.DataFrame(results)
    results.to_json(folder/file.replace('.txt', '_tcr_lifetimes_5P.json'))

if __name__ == '__main__':
     
    # User input
    folder = '/Users/koerbem/Documents/01_Klenerman_PhD/01_Projects/smoldyn/Model_9_5Phos_scan-k/'
    columnnames_tcrs = ['timestep', 'species_type', 'state', 'x', 'y', 'sn']
    #columnnames = ['time', 'tcr', 'tcrb', 'tcrp1', 'tcrp1b', 'tcrp2', 'tcrp2b', 'tcrp3', 'tcrp3b', 'tcrall', 'cd45i', 'cd45o', 'cd45all'] # for model 9 - 3Phos
    columnnames = ['time', 'tcr', 'tcrb', 'tcrp1', 'tcrp1b', 'tcrp2', 'tcrp2b', 'tcrp3', 'tcrp3b', 'tcrp4', 'tcrp4b', 'tcrp5', 'tcrp5b', 'tcrall', 'cd45i', 'cd45o', 'cd45all'] # for model 9 - 5PHos
    r_sim = 1 # µm
    r_c = 0.22 # µm
    A_c = r_c**2 * np.pi
    A_out = r_sim**2 * np.pi - A_c
    dt = 0.002 #s
    n_parallel = 10
    min_lifetime = 0.3 # s

    folder = Path(folder)

    filelist = []
    for ifile in listdir(folder):
        if ('tcrs' in ifile) or ('tcr_lifetimes' in ifile):
            continue
        elif Path(folder/ifile.replace('.txt', '_tcr_lifetimes.json')).exists():
            continue
        elif ('model_9' in ifile) and ('_600_' in ifile):
            filelist.append(ifile)

        # For debugging:
        #calculate_lifetime_5P(ifile,folder,columnnames_tcrs, columnnames, r_sim, r_c, A_c, A_out,dt,min_lifetime)
        
        if len(filelist) == n_parallel:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = executor.map(calculate_lifetime_5P, 
                                       filelist,
                                       repeat(folder),
                                       repeat(columnnames_tcrs), 
                                       repeat(columnnames), 
                                       repeat(r_sim), 
                                       repeat(r_c),    
                                       repeat(A_c), 
                                       repeat(A_out), 
                                       repeat(dt), 
                                       repeat(min_lifetime))
            filelist = []
    
    if len(filelist) > 0 :
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = executor.map(calculate_lifetime_5P, 
                                    filelist,
                                    repeat(folder),
                                    repeat(columnnames_tcrs), 
                                    repeat(columnnames), 
                                    repeat(r_sim), 
                                    repeat(r_c),    
                                    repeat(A_c), 
                                    repeat(A_out), 
                                    repeat(dt),
                                    repeat(min_lifetime))
        filelist = []
    
    print("Done!")
