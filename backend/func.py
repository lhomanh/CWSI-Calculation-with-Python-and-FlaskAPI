from pathlib import Path

from numba import njit
from numpy import abs, arange, array, diag, empty, inf, ones_like, sort, sum, polyfit, nan, sqrt, maximum, minimum, pi, zeros, ones
from numpy.linalg import lstsq
from pandas import read_excel, concat, DataFrame, Index, MultiIndex, date_range
import seaborn as sns
from scipy.signal import savgol_filter as savgol
from scipy.optimize import minimize
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
from pdb import set_trace

### Igno
import warnings
warnings.filterwarnings("ignore")


# from asce import asce
# from cwsi import cwsi
# from irls import irls

## Linh changed import functions
import backend.asce as asce
import backend.cwsi as cwsi
import backend.irls as irls
import numpy as np

# _DELTA_bar = 0.2340377704788737
# _Rnc_bar = 1.5236942756363387
# _a = 3.0353501931710416
# _b = -1.980071078099523
# _ra_bar = 13.339291717465803
# _rc_bar = 50.215548695601186
# _rho_bar = 0.987210057943333

def load(file_name, sheet):
    return read_excel(file_name, sheet_name=sheet)





def load_filename(file_name):

    #### Loading the information
    st = load(file_name, 'sites')
    wx = load(file_name, 'weather')
    idx = load(file_name, 'index')

    ### canopy temperature
    tc = load(file_name, 'canopy temperature') \
            .merge(idx, on='plot', how='right') \
            .groupby(['DOY', 'hour', 'treatment'])['Tc'] \
            .mean() \
            .reset_index()
            

    cc = load(file_name, 'canopy cover') \
            .merge(idx, on='plot', how='right') \
            .groupby(['DOY', 'treatment'])['cc'] \
            .mean() \
            .reset_index()

    ### SAP: provides a measurement of the transpiration of the plant
    sa = load(file_name, 'sap flow') \
            .merge(idx, on='plot', how='right') \
            .groupby(['treatment', 'DOY', 'hour'])['sap'] \
            .mean() \
            .reset_index() \
            .rename(columns={'sap': 'Sa'})

    return tc, cc, wx, st
        
def treatment_selection(tc, cc, wx, treatment = 'SWB_90'):
    
    wx['DOY'] = wx['day']               ### Main file Sc computation has problems without DOY in wx, Linh changed
    
    ### Filter the the calculation in the treatment
    cc_tmp = cc[cc['treatment'] == treatment]
    doy = Index(range(cc_tmp['DOY'].min(), cc_tmp['DOY'].max() + 1), name='DOY')
    cc_tmp = cc_tmp.set_index('DOY').reindex(doy).sort_index().interpolate() \
            .reset_index()
    tc_tmp = tc[tc['treatment'] == treatment]
    tmp = []
    for _, df in tc_tmp.groupby('DOY'):
            if len(df) == 24:
                tmp.append(df)
                
    tc_tmp = concat(tmp).sort_values(['DOY', 'hour'])
    days = set(cc_tmp['DOY']).intersection(tc_tmp['DOY']).intersection(wx['DOY'])
    cc_tmp = cc_tmp[cc_tmp['DOY'].isin(days)]
    tc_tmp = tc_tmp[tc_tmp['DOY'].isin(days)]
    wx_tmp = wx[wx['DOY'].isin(days)]

    ### There is a mask applied as a filter to dataframe wx
    df = tc_tmp.merge(wx_tmp, on=['DOY', 'hour'])
    
    return df

def CWSI_plot(df, wx, st, cc, treatment, start_date, end_date, compute_delta_bar = True, compute_Rnc_bar = True, compute_a_reg = True, compute_b_reg = True, compute_ra_bar = True, compute_rc_bar = True, compute_rho_bar = True, fixed = True):
    steps = end_date - start_date + 1
    days = np.linspace(start_date, end_date, steps)
    days = [int(item) for item in days]

    ### Data for the date 196 not available
    if 196 in days:
        days.remove(196)
    else:
        pass
    
    ## Start the calculation here:
        
    cwsi_days_empirical = []
    cwsi_days_theoretical = []

    for day in days:  
        # print("Calucalte on the day: ", str(day))
        df_oneday = df[df['DOY'] ==day]
        
        
        ### Start the compute here        
        wx['DOY'] = wx['day']               ### Main file Sc computation has problems without DOY in wx, Linh changed
        
        ### Site calculation. One value for the whole process
        st['P'] = P = asce.P(st['elev'][0])         ### atmospheric pressure (kPa)
        st['Lm'] = Lm = asce.Lm(st['lon'][0])     ### See function Lm: longitude of measurement site (positive degrees west of Greenwich England)
        st['gamma'] = gamma = asce.gamma(P)         #### psychrometric constant (kPa*C^-1)
        st['phi'] = phi = asce.phi(st['lat'][0])    ### latitude of measurement site (rad)
        Lz = st['Lz'][0]
        
        # ### Daily parameters: changing on daily basis
        df_oneday['Sc'] = Sc = asce.Sc(df_oneday['DOY'].values)    #### seasonal correction for solar time (h)
        df_oneday['delta'] = delta = asce.delta(df_oneday['DOY'].values)      ### solar declination (rad)
        df_oneday['omegas'] = omegas = asce.omegas(delta, phi)
        df_oneday['dr'] = dr = asce.dr(df_oneday['DOY'].values)
        
        hours = df_oneday['hour']
        Rs_d = df_oneday['Rs'].values
        Ta_d = df_oneday['Ta'].values
        Tc_d = df_oneday['Tc'].values
        RH_d = df_oneday['RH'].values
        u2_d = df_oneday['u2'].values
        
        
        omega_d = np.zeros([len(hours)])
        omega2_d = np.zeros([len(hours)])
        omega1_d = np.zeros([len(hours)])
        beta_d = np.zeros([len(hours)])
        Ra_d = np.zeros([len(hours)])
        Rso_d = np.zeros([len(hours)])
        fcd_d = np.zeros([len(hours)])
        es_d = np.zeros([len(hours)])
        ea_d = np.zeros([len(hours)])
        Rnl_d =  np.zeros([len(hours)])
        Rns_d =  np.zeros([len(hours)])
        Rn_d =  np.zeros([len(hours)])
        Cd_d =  np.zeros([len(hours)])
        G_d =  np.zeros([len(hours)])
        DELTA_d =  np.zeros([len(hours)])
        ETRs_d =  np.zeros([len(hours)])
        
        rho_d = np.zeros([len(hours)])
        Rnc_d = np.zeros([len(hours)])
        VPD_d = np.zeros([len(hours)])
        Tca_d = np.zeros([len(hours)])   #### 'Tc- Ta'
        
        ols_d = np.zeros([len(hours)])
        irls_d = np.zeros([len(hours)])
        
        #### Constant
        Cn = 66
        
        for i, hour in enumerate(hours):
            Sc_d = np.unique(Sc)[0]
            omegas_d = np.unique(omegas)[0]
            dr_d = np.unique(dr)[0]
            delta_d = np.unique(delta)[0]
            omega_d[i] = asce.omega(Lm, Lz, Sc_d,hour,1,1)
            omega2_d[i] = asce.omega2(omega_d[i], omegas_d, 1)
            omega1_d[i] = asce.omega1(omega_d[i], omega2_d[i], omegas_d, 1) 
            beta_d[i] = asce.beta(delta_d, omega_d[i], phi)   
            Ra_d[i] = asce.Ra(dr_d, delta_d, omega1_d[i], omega2_d[i], phi)
            Rso_d[i] = asce.Rso(Ra_d[i], st['elev'][0])
            
            if i == 0:
                fcd_d[i] = asce.fcd(Rs_d[i], Rso_d[i], beta_d[i], 0.055)   #### Initial value of fcd
            else:
                fcd_d[i] = asce.fcd(Rs_d[i], Rso_d[i], beta_d[i], fcd_d[i-1])
                
            es_d[i] = asce.es(Ta_d[i])
            ea_d[i] = asce.ea(RH_d[i], es_d[i])
            Rnl_d[i] = asce.Rnl(Ta_d[i], ea_d[i], fcd_d[i]) 
            Rns_d[i] = asce.Rns(Rs_d[i])
            Rn_d[i] = asce.Rn(Rnl_d[i], Rns_d[i])
            
            Cd_d[i] = asce.Cd(Rn_d[i])
            G_d[i] = asce.G(Rn_d[i])
            DELTA_d[i] = asce.DELTA(Ta_d[i])
            ETRs_d[i] = asce.ETsz(Cd_d[i], 66, DELTA_d[i], G_d[i], Rn_d[i], Ta_d[i], ea_d[i], es_d[i], gamma, u2_d[i])
            
            # CWSI coef: Or start the use of CWSI lib
            rho_d[i] = cwsi.rho(st['P'][0], Ta_d[i])
            Rnc_d[i] = cwsi.Rnc(Rn_d[i])
            
            VPD_d[i] = es_d[i] - ea_d[i]
            Tca_d[i] = Tc_d[i] - Ta_d[i]
                
        x = VPD_d
        y = Tca_d
        x = irls.design(x, 1)
        beta_irls, w = irls.irls(x, y)
        beta_ols, _ = irls.irls(x, y, weight=None)
        ols_d = x @ beta_ols
        irls_d = x @ beta_irls
        
        if compute_a_reg == True:
            a = beta_irls[0]
        else:
            a = 3.0353501931710416
        
        if compute_b_reg == True:
            b = beta_irls[1]
        else:
            b = -1.980071078099523
        
        
        #### The part of theoretical approach
        if compute_delta_bar == True:
            DELTA_bar = DELTA_d.mean()
        else:
            DELTA_bar = 0.2340377704788737
        
        if compute_Rnc_bar == True:
            Rnc_bar = Rnc_d.mean()
        else:
            Rnc_bar = 1.5236942756363387
            
        if compute_rho_bar == True:
            rho_bar = rho_d.mean()
        else:
            rho_bar = 0.987210057943333
        
        
        if compute_ra_bar == True:
            ra_bar = cwsi.ra(DELTA_bar, Rnc_bar, a, b,
                    rho_bar)
        else:
            ra_bar = 13.339291717465803
            
        if compute_rc_bar == True:
            rc_bar = cwsi.rc(DELTA_bar, b, st['gamma'][0], ra_bar)
        else:
            rc_bar = 50.215548695601186
            
        
        
        upper = cwsi.ul(Rnc_d, ra_bar,
                rho_d)
        lower = cwsi.ll(DELTA_d, ea_d,
                es_d, gamma, ra_bar, rc_bar, upper)
        
        cwsi_t = np.zeros([len(hours)]) 
        for i, hour in enumerate(hours):
            cwsi_t[i] = cwsi.cwsi(Ta_d[i], Tc_d[i], lower[i], upper[i])
        
        #### The part of emperical approach
        vpg_d = np.zeros([len(hours)])
        hi_E_d = np.zeros([len(hours)]) 
        low_E_d = np.zeros([len(hours)])  
        cwsi_e = np.zeros([len(hours)]) 
        for i, hour in enumerate(hours):
            vpg_d[i] = es_d[i] - asce.es(Ta_d[i] + a)
            hi_E_d[i] = a + b * vpg_d[i]
            low_E_d[i] = a + b * VPD_d[i]
            cwsi_e[i] = cwsi.cwsi(Ta_d[i], Tc_d[i], low_E_d[i], hi_E_d[i])
        cwsi_days_empirical.append(cwsi_e[13])
        cwsi_days_theoretical.append(cwsi_t[13])
        ### indent to end the calculation loop for days
    
    ### Plot the results here
    if fixed == True:
            plt.figure()
            plt.title("CWSI for treatment " + treatment + ' fixed')
            plt.plot(days, cwsi_days_theoretical, 'blue', label = "Theoretical")
            plt.plot(days, cwsi_days_empirical, 'r', label = "Emperical")
            plt.xlabel('Days')
            plt.ylabel('CWSI coefficients')
            plt.legend()
            plt.savefig("./static/result_2pm_fixed.png")
            # plt.show()
    if fixed == False:
            plt.figure()
            plt.title("CWSI for treatment " + treatment + ' no fix')
            plt.plot(days, cwsi_days_theoretical, 'blue', label = "Theoretical")
            plt.plot(days, cwsi_days_empirical, 'r', label = "Emperical")
            plt.xlabel('Days')
            plt.ylabel('CWSI coefficients')
            plt.legend()
            plt.savefig("./static/result_2pm_no_fix.png")
            # plt.show()
            
def plot_2pm(file_name, treatment, start_day, end_day):    
    tc, cc, wx, st = load_filename(file_name)
    df = treatment_selection(tc, cc, wx, treatment)
    CWSI_plot(df, wx, st, cc, treatment, start_day, end_day)
    CWSI_plot(df, wx, st, cc, treatment, start_day, end_day, compute_delta_bar = False, compute_Rnc_bar = False, compute_a_reg = False, compute_b_reg = False, compute_ra_bar = False, compute_rc_bar = False, compute_rho_bar = False, fixed = False)
    


