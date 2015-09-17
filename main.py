'''
------------------------------------------------------------------------
Last updated 9/17/2015

This program runs the steady state solver as well as the time path
iteration solution for the model with S-period lived agents, exogenous
labor, and two industries and two goods.

This Python script calls the following other file(s) with the associated
functions:
    ssfuncs.py
        feasible
        SS
------------------------------------------------------------------------
'''

# Import packages
import numpy as np
import scipy.optimize as opt
import ssfuncs as ssf
reload(ssf)

'''
------------------------------------------------------------------------
Declare parameters
------------------------------------------------------------------------
S           = integer in [3,80], number of periods an individual lives
alpha       = scalar in (0,1), expenditure share on good 1
c1til       = scalar > 0, minimum consumption of good 1
c2til       = scalar > 0, minimum consumption of good 2
cmtilvec    = [2,] vector, minimum consumption values for all goods
beta_ann    = scalar in [0,1), discount factor for one year
beta        = scalar in [0,1), discount factor for each model period
sigma       = scalar > 0, coefficient of relative risk aversion
nvec        = [S,] vector, exogenous labor supply n_{s,t}
A1          = scalar > 0, total factor productivity in industry 1
A2          = scalar > 0, total factor productivity in industry 2
Avec        = [2,] vector, total factor productivity values for all
              industries
gam1        = scalar in (0,1), capital share of income in industry 1
gam2        = scalar in (0,1), capital share of income in industry 2
gamvec      = [2,] vector, capital shares of income for all industries
eps1        = scalar in (0,1), elasticity of substitution between
              capital and labor in industry 1
eps2        = scalar in (0,1), elasticity of substitution between
              capital and labor in industry 2
epsvec      = [2,] vector, elasticities of substitution between capital
              and labor for all industries
del1_ann    = scalar in [0,1], one-year depreciation rate of capital in
              industry 1
del1        = scalar in [0,1], model period depreciation rate of capital
              in industry 1
del2_ann    = scalar in [0,1], one-year depreciation rate of capital in
              industry 2
del2        = scalar in [0,1], model period depreciation rate of capital
              in industry 2
delvec      = [2,] vector, model period depreciation rates for all
              industries
ss_tol      = scalar > 0, tolerance level for steady-state fsolve
ss_graphs   = boolean, =True if want graphs of steady-state objects
------------------------------------------------------------------------
'''
# Household parameters
S = int(80)
alpha = 0.4
c1til = 0.6
c2til = 0.6
cmtilvec = np.array([c1til, c2til])
beta_annual = 0.96
beta = beta_annual ** (80 / S)
sigma = 3.0
nvec = np.zeros(S)
nvec[:int(round(2 * S / 3))] = 1.
nvec[int(round(2 * S / 3)):] = 0.9
# Firm parameters
A1 = 1.
A2 = 1.2
Avec = np.array([A1, A2])
gam1 = 0.15
gam2 = 0.2
gamvec = np.array([gam1, gam2])
eps1 = 0.6
eps2 = 0.6
epsvec = np.array([eps1, eps2])
del1_ann = .04
del1 = 1 - ((1-del1_ann) ** (80 / S))
del2_ann = .05
del2 = 1 - ((1-del2_ann) ** (80 / S))
delvec = np.array([del1, del2])
# SS parameters
ss_tol = 1e-13
ss_graphs = True

'''
------------------------------------------------------------------------
Compute the steady state
------------------------------------------------------------------------
rbar_init    = scalar > 1, initial guess for steady-state model period
               interest rate
wbar_init    = scalar > 1, initial guess for steady-state real wage
rwbar_init   = [2,] vector, initial guesses for steady-state r and w
feas_params  = length 5 tuple, parameters for feasible function:
               (S, alpha, beta, sigma, ss_tol)
b_guess      = [S-1,] vector, initial guess for savings to use in fsolve
               in c11ssf.get_cbess
GoodGuess    = boolean, =True if initial steady-state guess is feasible
r_cstr_ss    = boolean, =True if initial r + delvec <= 0
w_cstr_ss    = boolean, =True if initial w <= 0
c_cstr_ss    = [S,] boolean vector, =True if c_s<=0 for initial r and w
cm_cstr_ss   = [2, S] boolean matrix, =True if c_{m,s}<=0 for initial r
               and w
K1K2_cstr_ss = boolean, =True if K1+K2<=0 for initial r and w
ss_params    = length 5 tuple, parameters for SS function:
               (S, alpha, beta, sigma, ss_tol)
r_ss         = scalar, steady-state interest rate
w_ss         = scalar > 0, steady-state wage
pm_ss        = [2,] vector, steady-state prices in each industry
p_ss         = scalar > 0, steady-state composite good price
b_ss         = [S-1,] vector, steady-state savings
c_ss         = [S,] vector, steady-state composite consumption
cm_ss        = [2,S] matrix, steady-state consumption of each good
eul_ss       = [S-1,] vector, steady-state Euler errors
Cm_ss        = [2,] vector, total demand for goods from each industry
Ym_ss        = [2,] vector, steady-state total output for each industry
Km_ss        = [2,] vector, steady-state capital demand for each industry
Lm_ss        = [2,] vector, steady-state labor demand for each industry
MCK_err_ss   = scalar, steady-state capital market clearing error
MCL_err_ss   = scalar, steady-state labor market clearing error
ss_time      = scalar, number of seconds to compute SS solution
rcmdiff_ss   = [2,] vector, steady-state difference in goods market
               clearing (resource constraint) in each industry
------------------------------------------------------------------------
'''

# Make sure initial guess of r and w is feasible
rbar_init = ((1 + 0.04) ** (80 / S)) - 1
wbar_init = 1.
rwbar_init = np.array([rbar_init, wbar_init])
feas_params = (S, alpha, beta, sigma, ss_tol)
b_guess = np.zeros(S-1)
b_guess[:int(round(2 * S / 3))] = \
    (np.linspace(0.003, 0.3, int(round(2 * S / 3))))
b_guess[int(round(2 * S / 3)):] = \
    (np.linspace(0.3, 0.003, S - 1 - int(round(2 * S / 3))))
GoodGuess, r_cstr_ss, w_cstr_ss, c_cstr_ss, cm_cstr_ss, K1K2_cstr_ss \
    = ssf.feasible(feas_params, rwbar_init, b_guess, cmtilvec, Avec,
    gamvec, epsvec, delvec, nvec)

if r_cstr_ss == True and w_cstr_ss == True:
    print 'Initial guess is not feasible because both r + delvec, w <= 0.'
elif r_cstr_ss == True and w_cstr_ss == False:
    print 'Initial guess is not feasible because r + delvec <= 0.'
elif r_cstr_ss == False and w_cstr_ss == True:
    print 'Initial guess is not feasible because w <= 0.'
elif (r_cstr_ss == False and w_cstr_ss == False and c_cstr_ss.max() == 1
  and K1K2_cstr_ss == False):
    print 'Initial guess is not feasible because c_s<=0 for some s.'
elif (r_cstr_ss == False and w_cstr_ss == False and c_cstr_ss.max() == 1
  and K1K2_cstr_ss == True):
    print 'Initial guess is not feasible because c_s<=0 for some s and K1+K2<=0.'
elif (r_cstr_ss == False and w_cstr_ss == False and c_cstr_ss.max() == 0
  and cm_cstr_ss.max() == 1 and K1K2_cstr_ss == False):
    print 'Initial guess is not feasible because c_{m,s}<=0 for some m and s.'
elif (r_cstr_ss == False and w_cstr_ss == False and c_cstr_ss.max() == 0
  and cm_cstr_ss.max() == 1 and K1K2_cstr_ss == True):
    print 'Initial guess is not feasible because c_{m,s}<=0 for some m and s and K1+K2<=0.'
elif (r_cstr_ss == False and w_cstr_ss == False and c_cstr_ss.max() == 0
  and cm_cstr_ss.max() == 0 and K1K2_cstr_ss == True):
    print 'Initial guess is not feasible because K1+K2<=0.'
elif GoodGuess == True:
    print 'Initial guess is feasible.'

    # Compute steady state
    print 'BEGIN STEADY STATE COMPUTATION'
    ss_params = (S, alpha, beta, sigma, ss_tol)
    (r_ss, w_ss, pm_ss, p_ss, b_ss, c_ss, cm_ss, eul_ss, Cm_ss, Ym_ss,
        Km_ss, Lm_ss, MCK_err_ss, MCL_err_ss, ss_time) = \
        ssf.SS(ss_params, rwbar_init, b_guess, cmtilvec, Avec, gamvec,
        epsvec, delvec, nvec, ss_graphs)

    # Print diagnostics
    print 'The maximum absolute steady-state Euler error is: ', \
        np.absolute(eul_ss).max()
    print 'The capital and labor market clearing errors are: ', \
        (MCK_err_ss, MCL_err_ss)
    print 'The steady-state distribution of capital is:'
    print b_ss
    print 'The steady-state distribution of composite consumption is:'
    print c_ss
    print 'The steady-state distribution of goods consumption is:'
    print cm_ss
    print 'The steady-state interest rate and wage:'
    print np.array([r_ss, w_ss])
    print 'Steady-state industry prices and composite price are:'
    print pm_ss, p_ss
    print 'Aggregate output, capital stock and consumption for each industry are:'
    print np.array([[Ym_ss], [Km_ss], [Cm_ss]])
    rcmdiff_ss = Ym_ss - Cm_ss - delvec * Km_ss
    print 'The difference Ym_ss - Cm_ss - delta_m * Km_ss is: ', rcmdiff_ss

    # Print SS computation time
    if ss_time < 60: # seconds
        secs = round(ss_time, 3)
        print 'SS computation time: ', secs, ' sec'
    elif ss_time >= 60 and ss_time < 3600: # minutes
        mins = int(ss_time / 60)
        secs = round(((ss_time / 60) - mins) * 60, 1)
        print 'SS computation time: ', mins, ' min, ', secs, ' sec'
    elif ss_time >= 3600 and ss_time < 86400: # hours
        hrs = int(ss_time / 3600)
        mins = int(((ss_time / 3600) - hrs) * 60)
        secs = round(((ss_time / 60) - mins) * 60, 1)
        print 'SS computation time: ', hrs, ' hrs, ', mins, ' min, ', secs, ' sec'
    elif ss_time >= 86400: # days
        days = int(ss_time / 86400)
        hrs = int(((ss_time / 86400) - days) * 24)
        mins = int(((ss_time / 3600) - hrs) * 60)
        secs = round(((ss_time / 60) - mins) * 60, 1)
        print 'SS computation time: ', days, ' days,', hrs, ' hrs, ', mins, ' min, ', secs, ' sec'


