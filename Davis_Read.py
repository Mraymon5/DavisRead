#%%
import os
import pandas as pd
import numpy as np
import time
import warnings
import matplotlib.pyplot as plt

#%%
def Davis_Read(folder=None, 
               keyword=".txt", 
               file=None, 
               intervalCutoff=50,
               burstCutoff=500,
               minBurstLength=1,
               output_folder=None, 
               ili_folder=None, 
               title=None, 
               plot_raw_interval=False, 
               plot_interval=False, 
               plot_cumulative=False, 
               timer=False, 
               stjohn_file=None):
    if 0:
        folder="/home/ramartin/Documents/Code/Python/Davis/Test_Data"
        keyword=".txt";
        file=None;
        intervalCutoff=50;
        burstCutoff=500;
        minBurstLength=1;
        output_folder="/home/ramartin/Documents/Code/Python/Davis/Test_Result";
        ili_folder=None;
        title=None;
        plot_raw_interval=False;
        plot_interval=False;
        plot_cumulative=False;
        timer=False;
        stjohn_file=None;
    
    if timer:
        start_time = time.time()

    if folder is None and file is None:
        raise ValueError("You must specify a file or a folder")
    
    if folder:
        file_names = [os.path.join(folder, f) for f in os.listdir(folder) if keyword in f]
        directory = pd.DataFrame({'Files': [f for f in os.listdir(folder) if keyword in f], 'Date': None, 'Animal': None})
    elif file:
        file_names = [file]
        directory = pd.DataFrame({'Files': [os.path.basename(file)], 'Date': None, 'Animal': None})
        if folder:
            warnings.warn("File overwrites Folder")
    
    directory['Files'] = directory['Files'].str.replace(".MS8.txt", "", regex=False)
    
    for i in range(len(directory)):
        directory.at[i, 'Date'] = f"{directory.at[i, 'Files'][:2]}/{directory.at[i, 'Files'][2:4]}"
        directory.at[i, 'Animal'] = directory.at[i, 'Files'][4:20]


    for fileN, fileIs in enumerate(file_names):#fileN = 0;fileIs = file_names[fileN]
        lick_export = []

        version = pd.read_csv(fileIs, header=None, sep=None, engine='python', nrows=1)
        
        if stjohn_file:
            temp_lat = pd.read_csv(stjohn_file, skiprows=3, nrows=1, delim_whitespace=True).values.T
            temp_lat = float(temp_lat[5] if temp_lat[6] == "s" else temp_lat[5]) * 1000
            licks_dav = pd.read_csv(stjohn_file, skiprows=9).values.T;
            meta = pd.read_csv(file, skiprows=5, nrows=2, sep=",", index_col=0, skipinitialspace=1)
            meta_data = pd.read_csv(file, skiprows=8, nrows=int(meta.at[2, 1]), sep=",", names=["PRES", "TUBE", "CONCENTRATION", "SOLUTION", "IPI", "LENGTH", "LICKS"], skipinitialspace=1)
            meta_data['Latency'] = temp_lat
        elif version.iat[0, 0] == 1.015:
            meta = pd.read_csv(fileIs, skiprows=4, nrows=2, header=None, index_col=0, sep=",", skipinitialspace=1)
            meta_data = pd.read_csv(fileIs, skiprows=7, nrows=int(meta.iat[1, 0]), sep=",", skipinitialspace=1)
            licks_dav = pd.read_csv(fileIs, skiprows=(8 + int(meta.iat[1, 0])), sep=",", header=None, index_col=0).T;
        else:
            meta = pd.read_csv(fileIs, skiprows=7, nrows=2, header=None, index_col=0, sep=",", skipinitialspace=1)
            meta_data = pd.read_csv(fileIs, skiprows=9, nrows=int(meta.iat[1, 0]), sep=",", skipinitialspace=1)
            licks_dav = pd.DataFrame(index=range(meta.iat[1, 0]), columns=range(max(meta_data['LICKS'])-1))
            for trialN in range(meta.iat[1, 0]):
                licks_dav.iloc[trialN,range(meta_data.loc[trialN,'LICKS']-1)] = pd.read_csv(fileIs, skiprows=(10 + int(meta.iat[1, 0]) + (trialN+1)), sep=",", header=None, index_col=0, skipinitialspace=1,nrows=1);
        
        licks_r = pd.DataFrame(licks_dav).T.reset_index(drop=True)
        licks_r.columns = range(licks_r.shape[1])
        
        NTrials = meta.iat[1, 0]
        
        for trialN in range(NTrials):
            if not pd.isna(licks_r.iloc[0, trialN]):
                licks_r.iloc[0, trialN] += meta_data.at[trialN, 'Latency']
                for lickN in range(1,licks_r.shape[0]):
                    if not pd.isna(licks_r.iloc[lickN, trialN]):
                        licks_r.iloc[lickN, trialN] += licks_r.iloc[lickN-1, trialN]
        
        licks_r = pd.concat([pd.DataFrame(meta_data['Latency']).T, licks_r], ignore_index=True)
        licks_r.iloc[0, licks_r.iloc[0, :] == (meta.iat[0, 0] * 1000)] = np.nan

        for trialN in range(licks_r.shape[1]):
            #burstCutoff = 145
            trial_lab = f"{meta_data.at[trialN, 'CONCENTRATION']} {meta_data.at[trialN, 'SOLUTION']}"
            
            if meta_data.at[trialN, 'LICKS'] == 0:
                lick_data = pd.DataFrame({'RawTime': [0]})
            else:
                lick_data = pd.DataFrame({'RawTime': licks_r.iloc[:, trialN].dropna()})
            lick_data['FilteredTime'] = int(0)
            lick_data['LickTotal'] = 0
            lick_data['BurstN'] = 0
            lick_data['BurstCount'] = 0
            lick_data['RawInterval'] = 0
            lick_data['FilteredInterval'] = 0
            lick_data['DoubleContact'] = 0
            lick_data = pd.concat([pd.DataFrame([[0]*8], columns=lick_data.columns), lick_data], ignore_index=1)
            
            #This Loop Filters licks for double contacts, and keeps track of bursts      
            for lickN in range(1, len(lick_data)):
                lick_data.loc[lickN,'RawInterval'] = lick_data.loc[lickN,'RawTime']-lick_data.loc[lickN-1,'RawTime']
                if lick_data.at[lickN, 'RawTime'] - lick_data.at[lickN-1, 'FilteredTime'] > intervalCutoff:
                    lick_data.at[lickN, 'FilteredTime'] = lick_data.at[lickN, 'RawTime']
                    lick_data.loc[range(lickN,len(lick_data)), 'LickTotal'] += 1
                    if lick_data.at[lickN, 'FilteredTime'] - lick_data.at[lickN-1, 'FilteredTime'] > burstCutoff or lickN == 1:
                        lick_data.loc[range(lickN,len(lick_data)), 'BurstN'] += 1
                        lick_data.at[lickN, 'BurstCount'] = 1
                        pass
                    else:
                        lick_data.at[lickN, 'BurstCount'] = lick_data.at[lickN-1, 'BurstCount']+1
                        pass
                else:
                    lick_data.at[lickN, 'FilteredTime'] = lick_data.at[lickN-1, 'FilteredTime']
                    lick_data.at[lickN, 'DoubleContact'] += lick_data.at[lickN-1, 'DoubleContact']+1
                    lick_data.at[lickN, 'BurstCount'] = lick_data.at[lickN-1, 'BurstCount']
                lick_data.loc[lickN,'FilteredInterval'] = lick_data.loc[lickN,'FilteredTime']-lick_data.loc[lickN-1,'FilteredTime']

            #Summarize Burst data
            NBursts = sum(lick_data['BurstN'].unique() != 0)
            burst_data = pd.DataFrame(  {'Licks': [None]*NBursts, 'Duration': [None]*NBursts, 'Start': [None]*NBursts, 'Stop': [None]*NBursts})            
            for BurstN in range(NBursts):
                BurstWhere = np.where(lick_data['BurstN'] == BurstN+1)
                burst_data['Licks'][BurstN] = lick_data['BurstCount'][np.array(BurstWhere).max()]
                burst_data['Start'][BurstN] = lick_data['FilteredTime'][np.array(BurstWhere).min()]
                burst_data['Stop'][BurstN] = lick_data['FilteredTime'][np.array(BurstWhere).max()]
                burst_data['Duration'][BurstN] = burst_data['Stop'][BurstN] - burst_data['Start'][BurstN]
            burst_data = burst_data.loc[burst_data['Licks'] >= minBurstLength,:]

            #Summarize Pause data
            NPauses = len(burst_data)-1
            pause_data = pd.DataFrame(  {'Duration': [None]*NPauses, 'Start': [None]*NPauses, 'Stop': [None]*NPauses})            
            for BurstN in range(NPauses):
                BurstWhere = np.where(lick_data['BurstN'] == BurstN+1)
                pause_data['Start'][BurstN] = lick_data['FilteredTime'][np.array(BurstWhere).max()]
                pause_data['Stop'][BurstN] = lick_data['FilteredTime'][np.array(BurstWhere).max()+1]
                pause_data['Duration'][BurstN] = pause_data['Stop'][BurstN] - pause_data['Start'][BurstN]

            #Compute MPI and Efficiency, and deal with no-lick trials
            if lick_data['LickTotal'].max() > 0:
                MPI = lick_data['FilteredInterval'][(lick_data['FilteredInterval'] > intervalCutoff) & (lick_data['FilteredInterval'] <= 160)].mean()
                Efficiency = len(lick_data['FilteredInterval'][(lick_data['FilteredInterval'] > intervalCutoff) & (lick_data['FilteredInterval'] <= 160)]) / (lick_data['LickTotal'].max())
                BurstLatency = burst_data['Start'][0]
                DoubleContacts = lick_data['DoubleContact'].sum()

            else:
                MPI = np.nan; Efficiency = np.nan;BurstLatency = np.nan; DoubleContacts = 0
            
            
            # Generate Output
            output = pd.DataFrame({
                'Trial': [trialN+1],#in brackets to fool pandas into thinking it isn't all scalars. Once is enough, I guess!
                'Solution': trial_lab,
                'TrialLength': meta_data['LENGTH'][trialN],
                'pre-IPI': meta_data['IPI  '][trialN],
                'Latency': meta_data['Latency'][trialN],
                'LickTotal': lick_data['LickTotal'].max(),
                'FirstMinLicks': lick_data['LickTotal'][lick_data['FilteredInterval']<60000].max(),
                'BurstTotal': len(burst_data),
                'BurstLicksMean': burst_data['Licks'].mean(),
                'BurstDurationMean': burst_data['Duration'].mean(),
                'BurstLatency': BurstLatency,
                'PauseTotal': len(pause_data),
                'PauseDurationMean': pause_data['Duration'].mean(),
                'MPI': MPI,
                'Efficiency': Efficiency,
                'DoubleContacts': DoubleContacts,
                'Date': directory.at[fileN, 'Date'],
                'Animal': directory.at[fileN, 'Animal'],
                'File': fileIs}
                )
# Save Output            
            lick_export.append(output) #Tag this trial's data onto the existing data

        lick_export_df = pd.concat(lick_export, ignore_index=True)

        if output_folder:
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
                warnings.warn('Specified output folder did not exist, making it')
            output_path = os.path.join(output_folder, f"{directory['Files'][fileN]}_lick_export.csv")
            lick_export_df.to_csv(output_path, index=False)
            print(f"Data exported to {output_path}")
    
        if ili_folder:
            if not os.path.exists(ili_folder):
                os.makedirs(ili_folder)
                warnings.warn('Specified ILI ouput folder did not exist, making it')
            output_path = os.path.join(ili_folder, f"{directory['Files'][fileN]}_ILI_export.csv")
            licks_dav.to_csv(output_path, index=False, header=0)
            print(f"Data exported to {output_path}")
        
        if timer:
            elapsed_time = time.time() - start_time
            print(f"Execution time: {elapsed_time:.2f} seconds")
               
    return lick_export_df

#%% Example usage
#output = Davis_Read(folder="/home/ramartin/Documents/Code/Python/Davis/Test_Data", output_folder="/home/ramartin/Documents/Code/Python/Davis/Test_Result", ili_folder="/home/ramartin/Documents/Code/Python/Davis/Test_Result")
