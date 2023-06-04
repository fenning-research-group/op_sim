# Deniz Cakan
# August 21, 2020
import numpy as np
import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import style
mpl.rcParams.update(mpl.rcParamsDefault)

def FDTD_load(R_fp, T_fp, PSK_fp, SI_fp, Spectrum_fp, RayTracing_fp, Gap_cutoff=101, monitor_rez=200, plot=True):
    S = pd.read_csv(Spectrum_fp, delimiter=',', header=0)
    WL_1, power_per_nm = S.iloc[:,0], S.iloc[:,1]
    z= np.zeros(2002)

    R = pd.read_csv(R_fp, delimiter=',', header=None)
    WL_R, R_per_nm = R.iloc[:,0], R.iloc[:,1]

    T = pd.read_csv(T_fp, delimiter=',', header=None)
    WL_2, T_per_nm1 = T.iloc[:,0], T.iloc[:,1]

    WL_0 = np.linspace(300, 1300, monitor_rez)
    WL_filltop = np.ones(monitor_rez)
    WL_fillbot = np.zeros(monitor_rez)

    SI_abs = pd.read_csv(SI_fp, delimiter=',', header=None)
    WL_SI_abs, SI_abs_per_nm = SI_abs.iloc[:,0], SI_abs.iloc[:,1]

    PSK_abs = pd.read_csv(PSK_fp, delimiter=',', header=None)
    WL_PSK_abs, PSK_abs_per_nm = PSK_abs.iloc[:,0], PSK_abs.iloc[:,1]

    # Transferring post E_gap PSK absorption into Si due to fitting error
    test_wl = np.array(WL_PSK_abs)
    test_psk_abs = np.array(PSK_abs_per_nm)
    test_SI_abs = np.array(SI_abs_per_nm)

    for n in range(Gap_cutoff, len(PSK_abs_per_nm)):
        test_SI_abs[n] = test_psk_abs[n] + test_SI_abs[n]

    test_psk_abs[Gap_cutoff:monitor_rez] = 0

    # PA_per_nm = ITO_abs_per_nm + HTL_abs_per_nm + ETL_abs_per_nm

    SI_abs_per_nm = test_SI_abs
    PSK_abs_per_nm = test_psk_abs

    # Adjusting T monitor to exclude silicon
    T_per_nm = T_per_nm1 + SI_abs_per_nm

    PA_per_nm_correction = T_per_nm+PSK_abs_per_nm
    PA = 1-R_per_nm - PA_per_nm_correction

    AM15 = np.interp(WL_PSK_abs, WL_1, power_per_nm)
    R_pabs = R_per_nm * AM15
    PA_pabs = (PA ) * AM15
    PSK_pabs = PSK_abs_per_nm * AM15
    T_pabs = T_per_nm * AM15

    R_pabs_stack = R_pabs
    PA_pabs_stack = R_pabs_stack + PA_pabs
    PSK_pabs_stack = PA_pabs_stack + PSK_pabs
    T_pabs_stack = PSK_pabs_stack + T_pabs

    WL_range = WL_1[40:1041]
    power_per_nm_range = power_per_nm[40:1041]
    Original_AM15=np.trapz(power_per_nm_range, x=WL_range)
    Rounded_AM15=np.trapz(AM15, x=WL_PSK_abs, axis=-1) #about 10% loss
    Loss=(Rounded_AM15-Original_AM15)/Original_AM15*100
    
    R_pabs_value=np.trapz(R_pabs, x=WL_PSK_abs)
    PA_pabs_value=np.trapz(PA_pabs, x=WL_PSK_abs)
    PSK_pabs_value=np.trapz(PSK_pabs, x=WL_PSK_abs)
    T_pabs_value=np.trapz(T_pabs, x=WL_PSK_abs)
    sum_pabs = R_pabs_value + PA_pabs_value + PSK_pabs_value + T_pabs_value
    accuracy_pabs= (sum_pabs-Rounded_AM15) /Rounded_AM15 *100
    
    # remember to manipulate export to be in lamda, absorption format
    A_Si = pd.read_csv(RayTracing_fp, delimiter=',', header=0)
    WL_3, Si_abs_nm1 = A_Si.iloc[:,0], A_Si.iloc[:,5]
    Si_subcell_abs_nm = np.interp(WL_PSK_abs, WL_3, Si_abs_nm1)
    Si_abs_per_nm = Si_subcell_abs_nm * T_per_nm
    TopBot_abs_per_nm = PSK_abs_per_nm + Si_abs_per_nm
    Si_Pabs = Si_abs_per_nm *AM15
    TopBot_abs_per_nm = PSK_abs_per_nm + Si_abs_per_nm
    Total_PA_per_nm = (1-R_per_nm) - TopBot_abs_per_nm
    Total_PA_pabs = Total_PA_per_nm * AM15

    R_photon=np.zeros(monitor_rez)
    PA_photon=np.zeros(monitor_rez)
    PSK_photon=np.zeros(monitor_rez)
    T_photon=np.zeros(monitor_rez)
    Si_photon=np.zeros(monitor_rez)
    Total_PA_photon=np.zeros(monitor_rez)

    for n in range(0, len(PSK_pabs)):
        R_photon[n]=R_pabs[n]/(1240/WL_PSK_abs[n])
        PA_photon[n]=PA_pabs[n]/(1240/WL_PSK_abs[n])
        PSK_photon[n]=PSK_pabs[n]/(1240/WL_PSK_abs[n])
        T_photon[n]=T_pabs[n]/(1240/WL_PSK_abs[n])
        Si_photon[n]=Si_Pabs[n]/(1240/WL_PSK_abs[n])
        Total_PA_photon[n]=Total_PA_pabs[n]/(1240/WL_PSK_abs[n])


    R_Jsc=np.trapz(R_photon, x=WL_PSK_abs)/10
    PA_Jsc=np.trapz(PA_photon, x=WL_PSK_abs)/10
    PSK_Jsc=np.trapz(PSK_photon, x=WL_PSK_abs)/10
    T_Jsc=np.trapz(T_photon, x=WL_PSK_abs)/10
    Si_Jsc=np.trapz(Si_photon, x=WL_PSK_abs)/10
    Total_PA_Jsc=np.trapz(Total_PA_photon, x=WL_PSK_abs)/10

    if plot == True:
        print('Cutoff =',round(WL_PSK_abs[Gap_cutoff],2), 'nm or', round(1240/WL_PSK_abs[Gap_cutoff],2), 'eV')
        plt.figure(4)
        plt.plot(WL_PSK_abs, R_pabs_stack, label='Reflection', color="gray", linestyle='--')
        plt.fill_between(WL_PSK_abs, WL_fillbot, R_pabs_stack, color='whitesmoke')
        plt.plot(WL_PSK_abs, PA_pabs_stack, label='Parasitic', color='gray')
        plt.fill_between(WL_PSK_abs, R_pabs_stack, PA_pabs_stack, color='lightgray')
        plt.plot(WL_PSK_abs, PSK_pabs_stack, label='Perovskite', color='dodgerblue')
        plt.fill_between(WL_PSK_abs, PA_pabs_stack, PSK_pabs_stack, color='azure')
        plt.plot(WL_PSK_abs, T_pabs_stack, label='Transmittance', color='lightcoral')
        plt.fill_between(WL_PSK_abs, PSK_pabs_stack, T_pabs_stack, color='mistyrose')
        plt.legend()
        plt.title('AM$1.5/0˚$ $P_{abs}$ Spectrum ')
        plt.ylim(0, max(T_pabs_stack)*1.12)
        plt.xlim(300,1200)
        plt.ylabel('$W m^{-2} * nm^{-1}$')
        plt.xlabel("λ Wavelength (nm)")
        plt.show()
    
        print('Original_AM15 =',round(Original_AM15,2),'W/m\u00b2')
        print('Rounded_AM15 =',round(Rounded_AM15,2),'W/m\u00b2')
        print('Error =',round(Loss,2),'%')

        # mA/cm\u00b2
        # can calculate PA_pabs by remainder theory 
        print('R_pabs_value =',round(R_pabs_value,2),'W/m\u00b2','or', round(R_pabs_value/Rounded_AM15*100,2),'%' )
        print('PA_pabs_value =',round(PA_pabs_value,2),'W/m\u00b2','or', round(PA_pabs_value/Rounded_AM15*100,2),'%' )
        print('PSK_pabs_value =',round(PSK_pabs_value,2),'W/m\u00b2','or', round(PSK_pabs_value/Rounded_AM15*100,2),'%' )
        print('T_pabs_value =',round(T_pabs_value,2),'W/m\u00b2','or', round(T_pabs_value/Rounded_AM15*100,2),'%' )
        print('Sum_pabs_error =',round(accuracy_pabs,2),'%' )

        plt.figure(1)
        plt.plot(WL_R, 1-R_per_nm, label="Reflection", color="gray", linestyle="--") #1-R
        plt.plot(WL_PSK_abs, PSK_abs_per_nm, label='Perovskite', color="dodgerblue") #Perovskite Absorption
        plt.plot(WL_PSK_abs, TopBot_abs_per_nm, label="Parasitic", color="gray", ) #Parasitic Absorption
        plt.plot(WL_2, Si_abs_per_nm, label='Silicon', color="lightcoral") # Tranmission

        plt.fill_between(WL_PSK_abs, TopBot_abs_per_nm, 1-R_per_nm, color='lightgray')
        plt.fill_between(WL_PSK_abs, WL_filltop, 1-R_per_nm,  color='whitesmoke')
        plt.fill_between(WL_PSK_abs, PSK_abs_per_nm, WL_fillbot, color='azure')
        plt.fill_between(WL_PSK_abs, Si_abs_per_nm, WL_fillbot, color='mistyrose' )

        plt.title("Total Absorption")
        plt.ylim(0,1)
        plt.yticks(np.arange(0,1.1,.1))  
        plt.xlim(300,1200)
        plt.xlabel("λ Wavelength (nm)")
        plt.ylabel("1-R and EQE (%)")
        plt.legend(bbox_to_anchor=(0,0), loc="lower left")

        # plt.legend()
        plt.show()

        # print('T_Jsc =',round(T_Jsc,2),'mA/cm\u00b2')
        # print('Topcell_PA_Jsc =',round(PA_Jsc,2),'mA/cm\u00b2')
        # print('Botomcell_PA_Jsc =',round(Total_PA_Jsc-PA_Jsc,2),'mA/cm\u00b2')
        # print('R_Jsc =',round(R_Jsc,2),'mA/cm\u00b2')
        # print('Total_PA_Jsc =',round(Total_PA_Jsc,2),'mA/cm\u00b2')
        # print('@@@@@','PSK_Jsc =',round(PSK_Jsc,2),'mA/cm\u00b2')
        # print('@@@@@','Si_Jsc =',round(Si_Jsc,2),'mA/cm\u00b2')

    mydict = {'R_Jsc': round(R_Jsc,2), 'PA_Jsc': round(PA_Jsc,2), 'PSK_Jsc':round(PSK_Jsc,2), 'Si_Jsc':round(Si_Jsc,2)}
    # DF = pd.DataFrame({'R_Jsc':[R_Jsc], 'PA_Jsc':Total_PA_Jsc, 'PSK_Jsc':PSK_Jsc, 'Si_Jsc':Si_Jsc})
    return mydict