import numpy as np
import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import style
mpl.rcParams.update(mpl.rcParamsDefault)
import csv

junction = int 


def SETFOS_load(Abs_profile_fp, junction=int, start_string = str):

    # columns = ['AOI', 'WL', 'R_tot', 'T_tot', 'A_tot', 'Rs_tot', 'Ts_tot', 'As_tot', 'Rp_tot', 'Tp_tot', 'Ap_tot', 'Haze', 'Anp', 'ITOtop_abs', 'SnOx_abs', 'PSK_abs', 'HTL_abs', 'ITObot_abs', 'Si_abs']
    columns = ['AOI', 'WL', 'R_tot', 'T_tot', 'A_tot', 'Rs_tot', 'Ts_tot', 'As_tot', 'Rp_tot', 'Tp_tot', 'Ap_tot', 'Haze', 'Anp', 'MgF2_abs', 'ITOtop_abs', 'PSK_abs', 'ITOmid_abs', 'Si_abs', 'Au_abs']

    # columns = ['AOI', 'WL', 'R_tot', 'T_tot', 'A_tot', 'Rs_tot', 'Ts_tot', 'As_tot', 'Rp_tot', 'Tp_tot', 'Ap_tot', 'Haze', 'Anp', 'Glass_abs', 'ITOtop_abs', 'PSK_abs', 'Au_abs']
    data = {c:[] for c in columns}
    data_idx = -1
    current_angle = None
    with open(Abs_profile_fp, 'r') as f:
        line = ''
        while not line.startswith(start_string): #need to seach your document and fill this 
           line = f.readline()
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:
            if len(row) == 0:
                continue
            this_angle = row[0]
            if current_angle != this_angle:
                 data_idx += 1
                 current_angle = this_angle
                 for c in columns:
                     data[c].append([])
            for col, val in zip(columns, row):
    #             data[col][data_idx].append(val)
                try:
                  writeval = float(val)
                except:
                  writeval = val
                data[col][data_idx].append(writeval)   

    for k, v in data.items(): data[k] = [np.array(v_) for v_ in v]

    df = pd.DataFrame(data)

    def remove_list(d):
        return(d[0])
    df['AOI'] = df['AOI'].apply(remove_list)
    
    if junction  == 2:
        df['PA'] = 1 - df['PSK_abs'] - df['Si_abs'] - df['R_tot']
    if junction == 1:
        df['PA'] = 1 - df['PSK_abs'] - df['R_tot']
    return df

def abs_to_pabs(Spectrum_fp, df, junction=int):

    S = pd.read_csv(Spectrum_fp, delimiter=',', header=0)
    WL_1, power_per_nm = S.iloc[:,0], S.iloc[:,1]

    WL = df['WL'][0]
    AM15 = np.interp(WL, WL_1, power_per_nm)

    df['R_pabs'] = np.array
    df['PA_pabs'] = np.array
    df['PSK_pabs'] = np.array
    if junction  ==2:
        df['Si_pabs'] = np.array
    
    df['R_pabs_stack'] = np.array
    df['PA_pabs_stack'] = np.array
    df['PSK_pabs_stack'] = np.array
    if junction  ==2:
        df['Si_pabs_stack'] = np.array
    
    for n in range(len(df)):
        df['R_pabs'][n]=AM15*df['R_tot'][n]
        df['PA_pabs'][n]=AM15*df['PA'][n]
        df['PSK_pabs'][n]=AM15*df['PSK_abs'][n]
        if junction  ==2:
            df['Si_pabs'][n]=AM15*df['Si_abs'][n]
        df['R_pabs_stack'][n] = df['R_pabs'][n]
        df['PA_pabs_stack'][n] = df['R_pabs'][n] + df['PA_pabs'][n]
        df['PSK_pabs_stack'][n] = df['PA_pabs_stack'][n] + df['PSK_pabs'][n]
        if junction  ==2:
            df['Si_pabs_stack'][n] = df['PSK_pabs_stack'][n] + df['Si_pabs'][n]

def pabs_integrate(Spectrum_fp, df, junction=int):
    
    df['R_pabs_value'] = np.float
    df['PA_pabs_value'] = np.float
    df['PSK_pabs_value'] = np.float
    if junction  ==2:
        df['Si_pabs_value'] = np.float
    
    df['R_Jsc'] = np.float
    df['PA_Jsc'] = np.float
    df['PSK_Jsc'] = np.float
    if junction  ==2:
        df['Si_Jsc'] = np.float
    
    df['R_photon'] = np.array
    df['PA_photon'] = np.array
    df['PSK_photon'] = np.array
    if junction  ==2:
        df['Si_photon'] = np.array

    monitor_rez = len(df['WL'][0])
    WL_0 = np.linspace(300, 1500, monitor_rez)
    WL_filltop = np.ones(monitor_rez)
    WL_fillbot = np.zeros(monitor_rez)
    
    for n in range(len(df)):      
        df['R_pabs_value'][n]=np.trapz(df['R_pabs'][n], x=df['WL'][n])
        df['PA_pabs_value'][n]=np.trapz(df['PA_pabs'][n], x=df['WL'][n])
        df['PSK_pabs_value'][n]=np.trapz(df['PSK_pabs'][n], x=df['WL'][n])
        if junction  ==2:
            df['Si_pabs_value'][n]=np.trapz(df['Si_pabs'][n], x=df['WL'][n])
        
        df['R_photon'][n]=df['R_pabs'][n]/(1240/df['WL'][n])
        df['PA_photon'][n]=df['PA_pabs'][n]/(1240/df['WL'][n])
        df['PSK_photon'][n]=df['PSK_pabs'][n]/(1240/df['WL'][n])
        if junction  ==2:
            df['Si_photon'][n]=df['Si_pabs'][n]/(1240/df['WL'][n])

        df['R_Jsc'][n]=np.trapz(df['R_photon'][n], x=df['WL'][n])/10
        df['PA_Jsc'][n]=np.trapz(df['PA_photon'][n], x=df['WL'][n])/10
        df['PSK_Jsc'][n]=np.trapz(df['PSK_photon'][n], x=df['WL'][n])/10
        if junction  ==2:
            df['Si_Jsc'][n]=np.trapz(df['Si_photon'][n], x=df['WL'][n])/10        

def SETFOS_plot(df, index=0, junction=int):
    monitor_rez = len(df['WL'][0])
    WL_0 = np.linspace(300, 1500, monitor_rez)
    WL_filltop = np.ones(monitor_rez)
    WL_fillbot = np.zeros(monitor_rez)
    
    plt.figure(0)
    plt.plot(df['WL'][index], df['R_pabs_stack'][index], label='Reflection', color="linen")
    plt.fill_between(df['WL'][index], WL_fillbot, df['R_pabs_stack'][index], color='whitesmoke')
    plt.plot(df['WL'][index], df['PA_pabs_stack'][index], label='Parasitic', color='gray')
    plt.fill_between(df['WL'][index], df['R_pabs_stack'][index], df['PA_pabs_stack'][index], color='linen')
    plt.plot(df['WL'][index], df['PSK_pabs_stack'][index], label='Top Cell', color='dodgerblue')
    plt.fill_between(df['WL'][index], df['PA_pabs_stack'][index], df['PSK_pabs_stack'][index], color='azure')
    if junction  ==2:
        plt.plot(df['WL'][index], df['Si_pabs_stack'][index], label='Bottom Cell', color='lightcoral')
        plt.fill_between(df['WL'][index], df['PSK_pabs_stack'][index], df['Si_pabs_stack'][index], color='mistyrose')
    plt.legend()
    plt.title('AM$1.5/0˚$ $P_{abs}$ Spectrum ')

    plt.ylim(0, max(df['PSK_pabs_stack'][index])*1.12)
    if junction  ==2:
        plt.ylim(0, max(df['Si_pabs_stack'][index])*1.12)

    if junction == 1:
        plt.xlim(300, 1200)
    if junction  == 2:
        plt.xlim(300,1500)
    plt.ylabel('$W m^{-2} * nm^{-1}$')
    plt.xlabel("λ Wavelength (nm)")
    plt.show()

    plt.figure(1)
    plt.plot(df['WL'][index], 1-df['R_tot'][index], label="Reflection", color="linen") #1-R
    if junction  ==1:
        plt.plot(df['WL'][index], (df['PA'][index]+df['PSK_abs'][index]), label="Parasitic", color="gray", ) #Parasitic Absorption

    if junction  ==2:
        plt.plot(df['WL'][index], (df['PA'][index]+df['PSK_abs'][index]+df['Si_abs'][index]), label="Parasitic", color="gray", ) #Parasitic Absorption

    plt.plot(df['WL'][index], df['PSK_abs'][index], label='Top Cell', color="dodgerblue") #Perovskite Absorption


    if junction  ==2:
        plt.plot(df['WL'][index], df['Si_abs'][index], label='Bottom Cell', color="lightcoral") # Tranmission
    if junction  ==2:
        plt.fill_between(df['WL'][index], (df['PA'][index]+df['PSK_abs'][index]+df['Si_abs'][index]), (df['PSK_abs'][index]+df['Si_abs'][index]), color='lightgray')


    if junction  ==1:
        plt.fill_between(df['WL'][index], (df['PA'][index]+df['PSK_abs'][index]), (df['PSK_abs'][index]), color='lightgray')


    plt.fill_between(df['WL'][index], WL_filltop, 1-df['R_tot'][index],  color='linen')
    plt.fill_between(df['WL'][index], df['PSK_abs'][index], WL_fillbot, color='azure')
    if junction  ==2:
        plt.fill_between(df['WL'][index], df['Si_abs'][index], WL_fillbot, color='mistyrose' )

    plt.title("Total Absorption")
    plt.ylim(0,1)
    plt.yticks(np.arange(0,1.1,.1))

    if junction == 1:
        plt.xlim(300, 1200)
    if junction  ==2:
        plt.xlim(300,1500)
    plt.xlabel("λ Wavelength (nm)")
    plt.ylabel("1-R and EQE (%)")
    plt.legend(bbox_to_anchor=(0,0), loc="lower left")
    plt.show() 
    

    mydict = {'R_Jsc': round(df['R_Jsc'][index],2), 'PA_Jsc': round(df['PA_Jsc'][index],2), 'Absorber_Jsc':round(df['PSK_Jsc'][index],2)}
    if junction  ==2:
        mydict = {'R_Jsc': round(df['R_Jsc'][index],2), 'PA_Jsc': round(df['PA_Jsc'][index],2), 'Top_Jsc':round(df['PSK_Jsc'][index],2), 'Bot_Jsc':round(df['Si_Jsc'][index],2)}
    return mydict