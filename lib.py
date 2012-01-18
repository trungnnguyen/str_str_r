"""
Library module
slope
__slope__
__windowed_rate_values__
__uniform_strain_range__
polynomial_fitting
polynomial_fitting_hr
polynomial_fitting_hr_multi
polynomial_fitting_hr_recursive
"""

def slope(x, y, upy, lowy):
    """
    Returns slope of the given stress and strain
    within the range indicated
    """
    import numpy as np
    # getting proper lower and upper index w.r.t. upy and lowy
    i = 0
    while True:
        if y[i]>lowy:
            ind1 = i
            break
        else : i = i + 1
    i = 0
    while True:
        if y[i]>upy:
            ind2 = i
            break
        else: i = i + 1
        
    # print 'Indices =', ind1, ind2
    # fitting and obtain the slope

    if ind1 == ind2 :
        ind2 = ind1 + 10
    z = np.polyfit(x[ind1:ind2], y[ind1:ind2],1)

    """
    try:
        z = np.polyfit(x[ind1:ind2],y[ind1:ind2],1)
    except TypeError:
        print x[ind1:ind2], y[ind1:ind2]
    """
    
    # returns the slope
    return z[0]

def __slope__(x, y):
    """
    Obtain the slope of y w.r.t. x
    """
    import numpy as np
    z = np.polyfit(x,y,1)
    return z[0]

def __windowed_rate_values__(
    fout, le, engie, sig, sig_engi, te,
    we, time, E_pl, delt,wrk):
    """
    Calculates windowsed rate values.

    Arguments
    ---------
      fout     : 'open' object for saving data, it should be clean new file.
      le       : logarithmic strain (+)
      engie    : engineering strain (+)
      sig      : True stress
      sig_engi : engineering stress
      te       : thickness strain   (-)
      we       : width strain       (-)
      time     : time
      E_pl     : plastic strain
      delt     : 1/2 size of the window
      wrk      : plastic work

    Dependent
    ---------
     1) numpy as np 2) __slope__

    Returns
     1) R-value 2) Strain rate 3) Hardeing rate 4) Engineering HR
    """
    import numpy as np
    
    ## new fout file
    ##----------------------------------------------------
    fout.writelines(
        '%13s %13s %13s %13s '%(
            'Strain','E_pl','Stress','Plastic_Work'))
    fout.writelines(
        '%13s %13s %13s %13s '%(
            'time','pl_SR','dSig/dEps','R-value'))
    fout.writelines(
        '%13s %13s %13s \n'%(
            'EngiStress','E_pl_trans', 'ER_pl_trans'))
    fout.writelines(
        '%13s %13s %13s %13s '%(
            'NA','NA','MPa','MJ/m^3'))
    fout.writelines(
        '%13s %13s %13s %13s '%(
            's','/s','MPa','NA'))
    fout.writelines(
        '%13s %13s %13s \n'%(
            'NA','NA','NA'))

    fname = fout.name
    fout.close()
    f = open(fname, 'r')
    header = f.read(); f.close()
    ##----------------------------------------------------

    ## array initialization is done by cloning le.
    hr = np.zeros(len(le))
    hr_engi = hr.copy()
    Rv = hr.copy()
    er_t = hr.copy()
    er_w = hr.copy()
    sr = hr.copy()
    srp = hr.copy()
    ##

    # window size: delt * 2
    # array : [i - delt : i + delt + 1]
    for i in range(len(le)):
        #only if the window is within the range
        if i-delt>=0 and i+delt+1<=len(le):
            
            HR      = __slope__(
                x=   le[i-delt:i+delt], y=     sig[i-delt:i+delt])
            HR_ENGI = __slope__(
                x=engie[i-delt:i+delt], y=sig_engi[i-delt:i+delt])
            RV      = __slope__(
                x=   te[i-delt:i+delt], y=      we[i-delt:i+delt])
            ER_T    = __slope__(
                x= time[i-delt:i+delt], y=      te[i-delt:i+delt])
            ER_W    = __slope__(
                x= time[i-delt:i+delt], y=      we[i-delt:i+delt])
            SR      = __slope__(
                x= time[i-delt:i+delt], y=      le[i-delt:i+delt])
            SRP     = __slope__(
                x= time[i-delt:i+delt], y=    E_pl[i-delt:i+delt])

            # Assign the slopes to the arrays
            hr[i] = HR
            hr_engi[i] = HR_ENGI
            Rv[i] = RV
            er_t[i] = ER_T
            er_w[i] = ER_W
            sr[i] = SR
            srp[i] = SRP
        else: # if out of the window,
            ## masking values
            hr[i] = np.NAN
            hr_engi[i] = np.NAN
            Rv[i] = np.NAN
            er_t[i] = np.NAN
            er_w[i] = np.NAN
            sr[i] = np.NAN  # strain rate
            srp[i] = np.NAN # plastic strain rate

    # transposed data
    data = np.array(
        [le, E_pl, sig, wrk, time,
         sr, hr, Rv, sig_engi, we, er_w]).T
    # no header allowed when writing data on files.
    np.savetxt(fname.split('.')[0]+'.st0', data,fmt='%13.6e') 
    cont = open(fname.split('.')[0]+'.st0', 'r').read()
    f = open(fname.split('.')[0]+'.st0', 'w')
    # prepend the header to the data-saved file.
    f.write(header);f.write(cont);f.close()

    return  Rv, sr, hr, hr_engi

def __uniform_strain_range__(force):
    """
    Provide the trim_index based on the relevant
    judgement upon the stress-strain curve.

    Arguments
    =========
    force   : Either force or engineering stress

    Dependents
    ==========
    numpy as np,

    Returns
    =======
    maxind
    """
    import numpy as np
    if type(force).__name__=='list':
        force = np.array(force)
    ## simply picking up the maximum force
    maxind = np.where(max(force)==force)[0][0]
    return maxind
    
def polynomial_fitting(strain, stress, lower_window, upper_window):
    """
    Within the windowed strain range, perform an polynomial
    fitting of the flow curve and return fitted flow curve,
    hardening rate, i.e. d(sig)/d(eps), parameter A and exponent.

    ## polynomial fitting related variables ----------------------
       sig = A eps^n
       ln(sig) = n ln(eps) + ln(A)
       -> y = nx + a
    ## -----------------------------------------------------------

    """
    import numpy as np
    for i in range(len(strain)):
        if strain[i]>lower_window: si=i;break
        pass
    for i in range(len(strain)):
        if strain[i]>upper_window: sf=i;break
        pass
    # sig = A eps^n : fitting polynomial function
    # ln(sig) = n ln(eps) + ln(A)
    # -> y = nx + A
    # z = [n,A]
    x = np.log(strain[si:sf])    
    y = np.log(stress[si:sf])
    z = np.polyfit(x,y,1) # linear fitting
    n = z[0]; A = z[1]
    A = np.exp(A)
    #print '%10s: %5.2f'%('exponent',n)
    #print '%10s: %5.2f'%('A',A)

    ## polynomial
    p = np.poly1d(z) # y = nx + A; ln(sig) = n ln(eps) + ln(A)
    ln_sig = p(np.log(strain))
    sig_fitted = np.exp(ln_sig)
    ## based on the polynomial's exponent and constant A,
    # d(sig)/d(eps) = A n eps ^(n-1)
    hp = A * n * strain ** (n-1)
    return sig_fitted, hp, A, n 


def polynomial_fitting_hr_recursive(
    hr, strain, lower_window=0.03, upper_window=0.10,
    order=3, maxiter=10,
    ):
    """
    Recursive polynomial fitting for hardening rate
    and return the optimized poly1d

    Arguments
    =========
     hr: hardening rate
     strain: 

    Dependents:
    ===========
     polynomial_fitting_hr_multi, 
     
    """
    ## initial guess
    pfhm = polynomial_fitting_hr_multi ## abbreviation
    hrp0, z = pfhm(
        hr=hr, strain=strain,
        lower_window=lower_window,
        upper_window=upper_window, order=order
        )

    ## hrp0: fitted hr for the given straina
    e0 = lower_window
    e1 = upper_window

    ## recursive calls.
    # calculates the strain where hrp0 = 0
    # -> that point is going to be the next upper_window
    # calculates the change of upper_window
    # break if iteration reach ...
    # break if change of upper_windo is negligible
    # <(0.0001 strain for instance.)
    ##

    ## recursive calls.
    count = 0
    while True:
        # calculates the strain where hrp0 = 0
        """the condition is based on the presumption that
        initial hrp0 is decreasing with respect to strain...

        hrp0: ln(hr) = N * ln(eps) + Bprime
        hr = exp (N * ln(eps) + Bprime)
        """
        if hrp0(e0) < hrp0(e1):
            raise IOError, 'unexpected behavior of hrp0'

        ## find the root. hrp0 = 0
        roots = hrp0.roots
        #print np.iscomplex(roots)
        if any(np.iscomplex(roots)==False):
            index = np.where(np.iscomplex(roots)==False)[0]
            if len(index)==1:
                #print 'single root'
                root = roots[index[0]]
            else: raise IOError,'multiple roots'

        #print root
        #raw_input()
        
        # the next upper_window
        
        e1_new = root # next upper_window
        most_recent_root = root
        del root
        
        # calculates the change of upper_window        
        delt = abs(e1 - e1_new)
        if delt<0.0001:
            print count
            break
        
        e1 = e1_new; del e1_new

        # new iteration
        hrp0, z = pfhm(
            hr=hr, strain=strain,
            lower_window=e0, upper_window=e1, order=order)
        count = count + 1
        # break if too many iterations
        if count> maxiter:
            print 'max iteration %i counts reached'%maxiter
            break

    print most_recent_root.real
    return hrp0, z, e1 #poly1d and polyfit: p = np.poly1d(z) 

def polynomial_fitting_hr(
    hr, strain, lower_window=None, upper_window=None):
    """
    Within the windwed strain range, perform an polynomial
    fitting of the hardening rate curve - d(sig)/d(eps)
    """
    import numpy as np
    if lower_window==None:si = 0
    else:
        for i in range(len(strain)):
            if strain[i]>lower_window:si=i;break

    if upper_window==None:sf = -1
    else:
        for i in range(len(strain)):
            if strain[i]>upper_window:sf=i;break

    hr = hr[si:sf]; eps = strain[si:sf]
    #d(sig)/d(eps) = A n eps^(n-1)
    # y = B x^(N)  :y = hr, x = eps, B = A * n, N = n-1
    # ln(y) = N ln(x) + ln(B)
    # Y = N X + B`: Y = ln(y), X = ln(x), B` = ln(B),N = n-1

    #Y = ln(y) = ln(hr)
    #X = ln(x) = ln(eps)
    Y = np.log(hr)
    X = np.log(eps)
    Z = np.polyfit(X,Y,1)
    
    p = np.poly1d(Z) #polynomial
    
    N = Z[0]
    Bprime = Z[1]
    
    B = np.exp(Bprime)
    n = N + 1
    A = B / n
    #print '%10s: %5.2f'%('exponent',n)
    #print '%10s: %5.2f'%('A',A)

    hr_fitted = p(np.log(strain))
    return hr_fitted, p, Bprime, N # p: ln(hr) = N * ln(eps) + Bprime

def polynomial_fitting_hr_multi(
    hr, strain, lower_window=None,
    upper_window=None, order=3):
    """
    Within the windwed (index wise) strain range,
    perform an n-order polynomial fitting of
    the hardening rate curve - d(sig)/d(eps)

    Returns
    =======
    p, z

    Arguments
    =========
     hr                   : hardening rate
     strain               : strain
     lower_window  = None : starting index for the strain array
     upper_window  = None : ending index for the strain array
     order         = 3    :
    """
    import numpy as np
    if lower_window==None: si = 0
    else:
        for i in range(len(strain)):
            if strain[i]>lower_window:si=i;break
            pass
        pass
    if upper_window==None: sf = -1
    else:
        for i in range(len(strain)):
            if strain[i]>upper_window:sf=i;break
            pass
        pass

    hr = hr[si:sf]
    eps = strain[si:sf]
    z = np.polyfit(eps, hr, order)
    p = np.poly1d(z)

    return p,z
