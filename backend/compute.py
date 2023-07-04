
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


# from asce import asce
# from cwsi import cwsi
# from irls import irls

import asce
import cwsi
import irls

_DELTA_bar = 0.2340377704788737
_Rnc_bar = 1.5236942756363387
_a = 3.0353501931710416
_b = -1.980071078099523
_ra_bar = 13.339291717465803
_rc_bar = 50.215548695601186
_rho_bar = 0.987210057943333

def load(sheet):
    return read_excel('LIRF_2019.xlsx', sheet_name=sheet)

def data():
    st = load('sites')
    wx = load('weather')
    idx = load('index')
    tc = load('canopy temperature') \
        .merge(idx, on='plot', how='right') \
        .groupby(['DOY', 'hour', 'treatment'])['Tc'] \
        .mean() \
        .reset_index()
        #.apply(lambda x: bw.iterative_mean(x.values)) \
    cc = load('canopy cover') \
        .merge(idx, on='plot', how='right') \
        .groupby(['DOY', 'treatment'])['cc'] \
        .mean() \
        .reset_index()
        #.apply(lambda x: bw.iterative_mean(x.values))
    sa = load('sap flow') \
        .merge(idx, on='plot', how='right') \
        .groupby(['treatment', 'DOY', 'hour'])['sap'] \
        .mean() \
        .reset_index() \
        .rename(columns={'sap': 'Sa'})
        #.apply(lambda x: bw.iterative_mean(x.values))
    return cc, sa, st, tc, wx

def compute(cc, sa, st, tc, wx):
    # refet
    st['Lm'] = Lm = asce.Lm(st['lon'][0])
    st['P'] = P = asce.P(st['elev'][0])
    st['gamma'] = gamma = asce.gamma(P)
    st['phi'] = phi = asce.phi(st['lat'][0])
    wx['Sc'] = Sc = asce.Sc(wx['DOY'].values)
    # wx['omega'] = omega = asce.omega(Lm, st['Lz'][0], Sc,
    #     wx['hour'].values)
    wx['omega'] = omega = asce.omega(Lm, st['Lz'][0], Sc,
        wx['hour'].values, 1,1)         # Linh changed
    wx['delta'] = delta = asce.delta(wx['DOY'].values)
    wx['omegas'] = omegas = asce.omegas(delta, phi)
    # wx['omega2'] = omega2 = asce.omega2(omega, omegas)
    # wx['omega1'] = omega1 = asce.omega1(omega, omega2, omegas)
    wx['omega2'] = omega2 = asce.omega2(omega, omegas, 1)
    wx['omega1'] = omega1 = asce.omega1(omega, omega2, omegas, 1)
    wx['beta'] = beta = asce.beta(delta, omega, phi)
    wx['dr'] = dr = asce.dr(wx['DOY'].values)
    wx['Ra'] = Ra = asce.Ra(dr, delta, omega1, omega2, phi)
    wx['Rso'] = Rso = asce.Rso(Ra, st['elev'][0])
    wx['fcd'] = fcd = asce.fcd(wx['Rs'].values, Rso, beta)
    wx['es'] = es = asce.es(wx['Ta'].values)
    wx['ea'] = ea = asce.ea(wx['RH'].values, es)
    wx['Rnl'] = Rnl = asce.Rnl(wx['Ta'].values, ea, fcd)
    wx['Rns'] = Rns = asce.Rns(wx['Rs'].values)
    wx['Rn'] = Rn = asce.Rn(Rnl, Rns)
    wx['Cd'] = Cd = asce.Cd(Rn)
    wx['G'] = G = asce.G(Rn)
    wx['DELTA'] = DELTA = asce.DELTA(wx['Ta'].values)
    wx['ETrs'] = asce.ETsz(Cd, 66, DELTA, G, Rn, wx['Ta'].values, ea, es,
        gamma, wx['u2'].values)
    # CWSI coef
    wx['rho'] = cwsi.rho(st['P'][0], wx['Ta'].values)
    wx['Rnc'] = cwsi.Rnc(wx['Rn'].values)
    cc_tmp = cc[cc['treatment'] == 'SWB_90']
    doy = Index(range(cc_tmp['DOY'].min(), cc_tmp['DOY'].max() + 1), name='DOY')
    cc_tmp = cc_tmp.set_index('DOY').reindex(doy).sort_index().interpolate() \
        .reset_index()
    tc_tmp = tc[tc['treatment'] == 'SWB_90']
    tmp = []
    for _, df in tc_tmp.groupby('DOY'):
        if len(df) == 24:
            tmp.append(df)
    tc_tmp = concat(tmp).sort_values(['DOY', 'hour'])
    days = set(cc_tmp['DOY']).intersection(tc_tmp['DOY']).intersection(wx['DOY'])
    cc_tmp = cc_tmp[cc_tmp['DOY'].isin(days)]
    tc_tmp = tc_tmp[tc_tmp['DOY'].isin(days)]
    wx_tmp = wx[wx['DOY'].isin(days)]
    mask = cwsi.mask(wx_tmp['Rs'].values, wx_tmp['Rso'].values,
        cc_tmp['cc'].values, wx_tmp['omega'].values, wx_tmp['omegas'].values)
    df = tc_tmp.merge(wx_tmp[mask], on=['DOY', 'hour'])
    df['VPD'] = x = df['es'] - df['ea']
    df['Tc-Ta'] = y = df['Tc'] - df['Ta']
    x = irls.design_p(x.values, 1)
    beta_irls, df['w'] = irls.b_irls(x, y)
    beta_ols, _ = irls.b_irls(x, y, weight=None)
    df['ols'] = x @ beta_ols
    df['irls'] = x @ beta_irls
    cwsi_ab_plot(df, beta_ols, beta_irls)
    st['a'], st['b'] = a, b = beta_irls
    st['DELTA_bar'] = DELTA_bar = df['DELTA'].mean()
    st['Rnc_bar'] = Rnc_bar = df['Rnc'].mean()
    st['rho_bar'] = rho_bar = df['rho'].mean()
    st['ra_bar'] = ra_bar = cwsi.ra(DELTA_bar, Rnc_bar, a, b,
        rho_bar)
    st['rc_bar'] = rc_bar = cwsi.rc(DELTA_bar, b, st['gamma'][0], ra_bar)
    # CWSI values
    wx['CWSI_u'] = upper = cwsi.upper(wx['Rnc'].values, ra_bar,
        wx['rho'].values)
    wx['CWSI_l'] = cwsi.lower(wx['DELTA'].values, wx['ea'].values,
        wx['es'].values, st['gamma'][0], ra_bar, rc_bar, upper)
    df = wx.merge(tc, on=['DOY', 'hour'])
    df['CWSI'] = cwsi_ = cwsi.cwsi(df['Ta'].values, df['Tc'].values,
        df['CWSI_l'].values, df['CWSI_u'].values)
    df['Ks'] = 1 - cwsi_
    ##### Questionable
    #df.loc[df['Ks'] < 0.2, 'Ks'] = 0.2
    #####
    # add sap
    df = df.merge(sa, on=['treatment', 'DOY', 'hour'])
    df['Sp'] = df['Sa'] / df['Ks']
    for trt, data in df.groupby('treatment'):
        trt = trt.split('_')[1]
        kcb = data.groupby('DOY').apply(
            # force regression through the origin
            lambda x: (x['ETrs'] @ x['Sp']) / (x['ETrs'] @ x['ETrs']))
        kcb = kcb.reset_index().rename(columns={0: 'Kcb_{}'.format(trt)})
        df = df.merge(kcb, on=['DOY'])
    for trt in (40, 65, 90):
        df['S_{}'.format(trt)] = df['ETrs'] * df['Kcb_{}'.format(trt)] \
            * df['Ks']
    df = df.sort_values(['treatment', 'DOY', 'hour'])
    df = df.replace([inf, -inf], nan)
    return df, st

def cwsi_ab_plot(df, beta_ols, beta_irls):
    # main plot
    df = df.sort_values('VPD')
    fig = plt.figure(figsize=(6.5, 4))
    plt.plot(df['VPD'].values, df['ols'].values,
        label='OLS: ' + str_eq(beta_ols), c='black', ls=':')
    plt.plot(df['VPD'].values, df['irls'].values,
        label='IRLS: ' + str_eq(beta_irls), c='black', ls='-')
    inc = df['w'] > 0
    s = df['w'][inc] * 50
    plt.scatter(df['VPD'][inc], df['Tc-Ta'][inc], edgecolors='black',
        marker='o', facecolors='none', s=s, label='Included')
    plt.scatter(df['VPD'][~inc], df['Tc-Ta'][~inc], color='black', marker='x',
        label='Excluded', s=16)
    plt.legend(frameon=False, fontsize='small', loc=3)
    plt.ylabel(r'$\mathregular{T_c - T_a}$ ($^\circ$C)')
    plt.xlabel(r'VPD (kPa)')
    plt.tight_layout()
    plt.show()
    #plt.savefig('../documents/two/figures/cwsi_coef.pdf')

def add_date(df):
    date = DataFrame({'date': date_range('2019', '2020', freq='d')[:-1]})
    date['DOY'] = date['date'].dt.dayofyear
    return df.merge(date, on='DOY')

def kcb_plots(df):
    df = add_date(df)
    ylim = max(df['S_90'].dropna().max(), df['S_65'].dropna().max(), df['S_40'].dropna().max())
    ylim = (-0.1 * ylim, 1.1 * ylim)
    c = {
        '90': 'tab:green',
        '65': 'tab:orange',
        '40': 'tab:red'}
    lim = 0
    for n, doy in df.groupby('DOY'):
        lim += 1
        if lim > 10:
            return
        t = 0
        date = doy['date'].iloc[0]
        fig, axs = plt.subplots(2, 3, figsize=(16,12), sharey='row',
            sharex='row', gridspec_kw={'wspace': 0})
        fig.suptitle('{} ({})'.format(date, n))
        #    doy.loc[doy['treatment'] == 'SWB_40', 'cc'][0].round(3),
        #    doy.loc[doy['treatment'] == 'SWB_65', 'cc'][0].round(3),
        #    doy.loc[doy['treatment'] == 'SWB_90', 'cc'][0].round(3)))
        for trt_name, trt in doy.groupby('treatment'):
            n = trt_name.split('_')[1]
            # upper
            ax = axs[0][t]
            ax.set_title(trt_name, color=c[n])
            ax.scatter(trt['hour'], trt['Sa'], c=c[n])
            #ax.plot(trt['hour'], trt['Sp' + n], c=c[n], ls=':')
            ax.plot(trt['hour'], trt['S_90'], c='tab:green', label=None)
            ax.plot(trt['hour'], trt['S_65'], c='tab:orange', label=None)
            ax.plot(trt['hour'], trt['S_40'], c='tab:red', label=None)
            ax.set_ylim(ylim)
            ax.set_ylabel('S (mm)')
            ax.set_xlabel('hour')
            ax.legend(frameon=False, fontsize='small', loc=2)
            # lower
            ax = axs[1][t]
            trt = trt.sort_values('ETrs')
            ax.scatter(trt['ETrs'], trt['Sp'], c=c[n])
            ax.plot(trt['ETrs'], trt['ETrs'] * trt['Kcb_' + n], c=c[n],
                label='Kcb={}'.format(trt['Kcb_' + n].iloc[0].round(3)))
            ax.set_xlabel('ETrs (mm)')
            ax.set_ylabel('Sp (mm)')
            ax.legend(frameon=False, fontsize='small', loc=2)
            t += 1
        plt.show()

def test_plots(wx, cc, gamma, P):
    eq1 = r'$y = {:.2f}$'
    eq2 = r'$y = {:.2f} x + {:.2f}$'
    mask = cwsi.mask(wx['Rs'].values, wx['Rso'].values, cc['cc'].values,
        wx['omega'].values, wx['omegas'].values)
    wx = wx[mask]
    # 0 order polynomial
    DELTA_bar = wx['DELTA'].mean()
    Rnc = cwsi.Rnc(wx['Rn'].values)
    Rnc_bar = Rnc.mean()
    rho = cwsi.rho(P, wx['Ta'].values)
    rho_bar = rho.mean()
    ra_bar = cwsi.ra_bar(DELTA_bar, Rnc_bar, a, b, rho_bar)
    rc_bar = cwsi.rc_bar(DELTA_bar, b, gamma, ra_bar)
    return DELTA_bar, Rnc_bar, ra_bar, rc_bar, rho_bar
    # first order
    #grpr = wx.groupby('DOY')
    #DELTA = grpr['DELTA'].mean()
    #Rnc = grpr.apply(lambda x: (x['Rn'] - x['G']).mean())
    #rho = 3.484 * P / (grpr['Ta'].mean() + 273.15)
    #ra = cwsi.ra_bar(DELTA.values, Rnc.values, a, b, rho.values)
    #rc = cwsi.rc_bar(DELTA.values, b, gamma, ra)
    #rab = irls.yhat(irls.beta(doy, ra), doy)
    #rcb = irls.yhat(irls.beta(doy, rc), doy)
    # first order
    all_doy = wx['DOY']

    DELTA = wx['DELTA'].values
    Rnc = (wx['Rn'] - wx['G']).values
    rho = (3.484 * P / (wx['Ta'] + 273.15)).values
    ra = cwsi.ra_bar(DELTA, Rnc, a, b, rho)
    rc = cwsi.rc_bar(DELTA, b, gamma, ra)
    doy = wx['DOY'].values
    rab, raw = irls.beta(doy, ra)
    rcb, rcw = irls.beta(doy, rc)
    # plots
    fig, (ax1, ax2) = plt.subplots(2, figsize=(6.5, 8), sharex='all',
        gridspec_kw={'hspace': 0, 'wspace': 0})
    # ra
    inc = raw > 0
    s = raw[inc] * 50
    ax1.scatter(doy[inc], ra[inc], edgecolors='black', marker='o',
        facecolors='none', s=s, label='Included')
    ax1.scatter(doy[~inc], ra[~inc], color='black', marker='x',
        label='Excluded', s=16)
    ax1.plot(all_doy, irls.yhat(rab, all_doy), color='black',
        label=eq2.format(*rab[::-1]))
    ax1.axhline(ra_bar, color='black', linestyle=':', label=eq1.format(ra_bar))
    ax1.set_ylabel(r'$\mathregular{r_a\;(s/m)}$')
    ax1.legend(frameon=False, fontsize='small')
    # rc
    inc = rcw > 0
    s = rcw[inc] * 50
    ax2.scatter(doy[inc], rc[inc], edgecolors='black', marker='o',
        facecolors='none', s=s, label='Included')
    ax2.scatter(doy[~inc], rc[~inc], color='black', marker='x',
        label='Excluded', s=16)
    ax2.plot(all_doy, irls.yhat(rcb, all_doy), color='black',
        label=eq2.format(*rcb[::-1]))
    ax2.axhline(rc_bar, color='black', linestyle=':', label=eq1.format(rc_bar))
    ax2.set_ylabel(r'$\mathregular{r_c\;(s/m)}$')
    ax2.set_xlabel(r'DOY')
    ax2.legend(frameon=False, fontsize='small')
    plt.tight_layout()
    plt.savefig('figures/cwsis_coef.pdf')
    plt.show()




def cwsi_data(tc, info, wx, sap, cc):
    # canopy temperature
    tc = tc.groupby(['DOY', 'hour', 'treatment'])['Tc'].apply(lambda x: bw.iterative_mean(x.values))
    tc = tc.reset_index()
    # weather
    Lm = asce.Lm(info['lon'].values[0])
    P = asce.P(info['elev'].iloc[0])
    phi = asce.phi(info['lat'].values[0])
    gamma = asce.gamma(P)
    es = asce.es(wx['Ta'].values)
    ea = asce.ea(wx['RH'].values, es)
    Sc = asce.Sc(wx['DOY'].values)
    omega = asce.omega(Lm, info['Lz'].values[0], Sc, wx['hour'].values)
    delta = asce.delta(wx['DOY'].values)
    omegas = asce.omegas(delta, phi)
    omega2 = asce.omega2(omega, omegas)
    omega1 = asce.omega1(omega, omega2, omegas)
    beta = asce.beta(delta, omega, phi)
    dr = asce.dr(wx['DOY'].values)
    Ra = asce.Ra(dr, delta, omega1, omega2, phi)
    Rso = asce.Rso(Ra, info['elev'].values[0])
    fcd = asce.fcd(wx['Rs'].values, Rso, beta)
    Rnl = asce.Rnl(wx['Ta'].values, ea, fcd)
    Rns = asce.Rns(wx['Rs'].values)
    Rn = asce.Rn(Rnl, Rns)
    Cd = asce.Cd(Rn)
    Cn = 66
    G = asce.G(Rn)
    DELTA = asce.DELTA(wx['Ta'].values)
    ETrs = asce.ETsz(Cd, Cn, DELTA, G, Rn, wx['Ta'].values, ea, es, gamma,
        wx['u2'].values)
    ETrs[ETrs < 0] = 0
    wx['ETrs'] = ETrs
    # theoretical
    Rnc = cwsi.Rnc(Rn)
    rho = cwsi.rho(P, wx['Ta'].values)
    wx['hi-T'] = upper = cwsi.upper(Rnc, ra_bar, rho)
    wx['lo-T'] = cwsi.lower(DELTA, ea, es, gamma, ra_bar, rc_bar, upper)
    # empirical
    wx['VPD'] = vpd = es - ea
    vpg = es - asce.es(wx['Ta'].values + a)
    wx['hi-E'] = a + b * vpg
    wx['lo-E'] = a + b * vpd
    # sap flow
    tmp = []
    for keys, df in sap.groupby(['treatment', 'DOY']):
        df = df.groupby('hour')['sap'].apply(
            lambda x: bw.iterative_mean(x.values))
        df = df.reset_index()
        if len(df) == 24:
            df['treatment'], df['DOY'] = keys
            tmp.append(df)
    sap = concat(tmp)
    # cover
    cc = cc.groupby(['DOY', 'treatment'])['cc'].apply(lambda x: bw.iterative_mean(x.values))
    ix = cc.index.levels[0]
    ix = Index(range(ix.min(), ix.max() + 1), name='DOY')
    tmp = []
    for trt, data in cc.groupby(['treatment']):
        data = data.reset_index(level=1)
        idx = Index(range(data.index.min(), data.index.max() + 1), name='DOY')
        data = data.reindex(ix).sort_index().interpolate().reset_index()
        data['treatment'] = trt
        tmp.append(data)
    cc = concat(tmp)
    # combine
    df = wx.merge(tc, on=['DOY', 'hour'])
    df = df.merge(cc, on=['treatment', 'DOY'])
    df = df.merge(sap, on=['treatment', 'DOY', 'hour'])
    df = df.sort_values(['treatment', 'DOY', 'hour'])
    # indices
    df['CWSI-E'] = cwsi.cwsi(df['Ta'].values, df['Tc'].values,
        df['lo-E'].values, df['hi-E'].values)
    df['CWSI'] = cwsi.cwsi(df['Ta'].values, df['Tc'].values,
        df['lo-T'].values, df['hi-T'].values)
    #df.loc[df['VPD'] < 1.5, 'CWSI-T'] = 0
    #df.loc[df['VPD'] < 1.5, 'CWSI-E'] = 0
    df.loc[df['CWSI'] > 0.9, 'CWSI'] = 0.9
    #df.loc[df['CWSI-E'] > 0.9, 'CWSI-E'] = 0.9
    df.loc[df['hi-T'] < 0, 'CWSI'] = 0
    df['Ks'] = 1 - df['CWSI']
    #df['Ks-E'] = 1 - df['CWSI-E']
    for trt, data in df.groupby('treatment'):
        data = data.copy()
        n = trt.split('_')[1]
        c_ETp = 'ETp_' + n
        c_Kcb = 'Kcb_' + n
        data[c_ETp] = data['sap'] / data['Ks']
        for doy, tmp in data.groupby('DOY'):
            x = ones((len(tmp), 2))
            x[:, 1] = tmp['ETrs'].values
            Kcb, _ = irls.fit(x, tmp[c_ETp].values, weight=None)
            #Kcb = tmp[c_ETp].sum() / tmp['ETrs'].sum()
            data.loc[data['DOY'] ==  doy, c_Kcb] = Kcb[1]
        df = df.merge(data[['DOY', 'hour', c_Kcb, c_ETp]], on=['DOY', 'hour'])
        df['Tc_' + n] = df['ETrs'] * df[c_Kcb] * df['Ks']
    return df


def cwsi_plots(df):
    sap_lim = df['sap'].max()
    et_lim = df['ETrs'].max()
    lim = max(sap_lim, et_lim)
    df = add_date(df)
    for doy, data in df.groupby('DOY'):
        date = data['date'].iloc[0]
        #fig, ax1 = plt.subplots()
        #ax1.plot(data['hour'], data['Rs'], label='Rs', color='tab:red')
        fig, axs = plt.subplots(3, figsize=(6.5, 9), sharex='all',
            gridspec_kw={'hspace': 0, 'wspace': 0})
        fig.suptitle('DOY {}: {}'.format(doy, data[data['treatment'] == 'SWB_90']['cc'].iloc[0]))
        c = ('tab:red', 'tab:orange', 'tab:green')
        for n, trt in enumerate(['SWB_40', 'SWB_65', 'SWB_90']):
            ax1 = axs[n]
            ax2 = ax1.twinx()
            tmp = data[data['treatment'] == trt]
            if not tmp.empty:
                print(str(date.date()) + ' ' + trt)
                print('=' * 17)
                print(tmp[['hour', 'Tc', 'Ta', 'hi-E', 'hi-T', 'lo-E', 'lo-T', 'CWSI-E', 'CWSI-T']].set_index('hour'))
                print()
            ax1.set_ylim(-0.1, 1.1)
            ax1.set_ylabel('CWSI')
            ax1.plot(tmp['hour'], tmp['CWSI-T'], c=c[n])
            ax2.plot(tmp['hour'], tmp['sap'], c='tab:blue')
            ax2.plot(tmp['hour'], tmp['ETrs'], c='tab:blue', ls=':')
            ax2.plot(tmp['hour'], tmp['ETp-T'], c='tab:purple')
            ax2.set_ylabel('mm')
            ax2.set_ylim(-lim / 10, lim + lim / 10)
        #plt.legend()
        plt.show()

def cwsi_plot_days(cc, tc, wx, bi, mi):
    eq = r'$y = {:.2f} x + {:.2f}$'
    df = tc.merge(wx, on=['DOY', 'hour'])
    df['flag'] = cwsi.flags(wx['Rs'].values, wx['Rso'].values, cc['cc'].values,
        wx['omega'].values, wx['omegas'].values, th_cc=0)
    df['VPD'] = df['es'] - df['ea']
    df['dif'] = df['Tc'] - df['Ta']
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 3.5), sharex='all',
        sharey='all', gridspec_kw={'hspace': 0, 'wspace': 0})
    letters = ('a', 'b')
    for n, doy in enumerate((184, 210)):
        ax = axs[n]
        data = df[df['DOY'] == doy]
        day = data[data['flag']]
        night = data[~data['flag']]
        for i in range(len(data) - 1):
            x1, y1 = data.iloc[i][['VPD', 'dif']]
            x2, y2 = data.iloc[i + 1][['VPD', 'dif']]
            dx, dy = x2 - x1, y2 - y1
            ax.arrow(x1, y1, dx, dy, ls=':')
        ax.plot(data['VPD'], data['VPD'] * mi + bi, color='black',
            label=eq.format(mi, bi))
        ax.scatter(day['VPD'], day['dif'], edgecolors='black', marker='o',
            facecolors='none', label='Daytime')
        ax.scatter(night['VPD'], night['dif'], color='black', marker='x',
            label='Non-Daytime')
        ax.annotate('({})'.format(letters[n]), (0.05, 0.04),
            xycoords='axes fraction', ha='left', va='bottom')
        if n == 0:
            ax.set_ylabel(r'$\mathregular{T_c - T_a}$ ($^\circ$C)')
            ax.set_xlabel(r'VPD (kPa)', ha='center', position=(1,0))
        if n == 1:
            ax.legend(frameon=False, fontsize='small', loc=1)
    plt.tight_layout()
    #plt.show()
    plt.savefig('../documents/two/figures/cwsi_days.pdf')

def str_eq(b):
    res = '{:.2f}'.format(b[0])
    n = 1
    for p in b[1:]:
        res += '{:+.2f}x'.format(p)
        if n > 1:
            res += '^{}'.format(n)
        n += 1
    return '$y={}$'.format(res)
