##
## 2011 Nov Steglich's Stresss-Strain curve calculation
##

# data file reader
def single_reader(filename='POSTE31_TE_0%F8_09.dat'):
    """
    2011 Nov
    """
    f = open(filename, 'r')
    cont = f.read()

    ld = cont.split('\r\n')
    #print ld; raw_input()

    data = []
    ld = ld[1:] # Remove header
    for i in range(len(ld)):
        temp = ld[i].split('\t')
        if len(temp)> 3: data.append(ld[i].split('\t'))
        else: pass
        pass
    return data

def multi_reader():
    import glob
    data = []
    files = glob.glob('*.dat')
    print 'files:\n', files
    for f in files:
        data.append(single_reader(f))
        pass
    return data

# data file analysis
def analysis(single_data, delt=100):
    """
    Single data
    """
    from str_str_r import __slope__ # slope
    import numpy as np
    from str_str_r import slope  # modulus slope

    # strain, standard force, test time, reduct. in wid, Lankford, eng.strain, eng. Stress
    # what is necessary: eng. eps, eng. sig, time, Lankford,
    if type(single_data).__name__=='ndarray':
        pass
    elif type(single_data).__name__=='list':
        single_data = np.array(single_data)
        pass
    else: raise IOError, 'Unexpected type'

    ## transpose
    single_data = single_data.T
    # column -> Var
    force = single_data[1]  # force
    time =  single_data[2]
    R    =  single_data[4]
    eps  =  single_data[5]  # engi
    sig  =  single_data[6]  # engi

    # ',' -> '.'
    force = comma2dot(force)
    time  = comma2dot(time)
    R     = comma2dot(R)
    eps   = comma2dot(eps)
    sig   = comma2dot(sig)
    # empty string trimmer
    force = empty_string_trimmer(force)
    time  = empty_string_trimmer(time)
    R     = empty_string_trimmer(R)
    eps   = empty_string_trimmer(eps)
    sig   = empty_string_trimmer(sig)
    
    tiny = 10**-5
    mx = max(force)
    for i in range(len(force)):
        if abs(force[i]-mx)<tiny:
            maxind = i
            break
        pass
    
    #maxind = np.where(max(force)==force)[0][0]
    force = np.array(force[:maxind])
    time  = np.array(time[:maxind])
    R     = np.array(R[:maxind])
    eps   = np.array(eps[:maxind])
    sig   = np.array(sig[:maxind])
    

    sig = sig * (eps + 1)  # log stress
    eps = np.log(eps + 1)  # log strain

    ## modulus calculation
    print eps, len(eps)
    #raw_input()
    print sig, len(sig)
    #raw_input()
    mod = slope(x=eps, y=sig, upy = 80, lowy = 20)

    print 'modulus: %5.2f [GPa]'%(mod/1000.)
    eps_pl = []
    for i in range(len(eps)):
        temp = eps[i] - sig[i]/mod
        if temp<0: eps_pl.append(0.)
        else: eps_pl.append(temp)
        pass    

    ## strain rate calculation
    epsr = []
    for i in range(len(eps)):
        if i-delt>=0 and i+delt+1<=len(eps):
            sr = __slope__(x=time[i-delt:i+delt], y=eps[i-delt:i+delt])
            epsr.append(sr)
            pass
        else: epsr.append(np.NAN)
        pass

    ## plastic region trimmer
    for i in range(len(eps_pl)):
        if eps_pl[i]>0.002:
            maxi = i; break
            pass
        pass
    eps_pl = eps_pl[maxi:]
    sig    =    sig[maxi:]
    R      =      R[maxi:]
    epsr = np.array(epsr)

    return eps_pl, sig, R  #plastic strain, True Stress, R-value
    
def comma2dot(list):
    """
    Some European countries use ',' instead of '.'
    when separating subzero values from above zero ones.
    ex.) 1.23 -> 1,23
    """
    dum = []
    for i in range(len(list)):
        dum.append(list[i].replace(',','.'))
        pass
    return dum

def empty_string_trimmer(x):
    """Empty string element trimmer"""
    import numpy as np
    dum = []
    for i in range(len(x)):
        try: float(x[i])
        except: dum.append(np.NAN)
        else: dum.append(float(x[i]))
        pass
    return dum

def main(ifig=1, jump=1):
    """
    """
    import matplotlib.pyplot as plt
    import numpy as np

    data = multi_reader()
    eps = []
    sig = []
    R   = []

    # Plastic strain, Stress, R-value
    for i in range(len(data)):
        e, s, rv = analysis(data[i])
        eps.append(e); sig.append(s); R.append(rv)
        pass

    ## save to files
    for i in range(len(data)):
        fn = 'str_%s.txt'%str(i).zfill(2)
        print fn
        np.savetxt(fn, np.array([eps[i],sig[i]]).T)
        fn = 'R_%s.txt'%str(i).zfill(2)
        print fn
        np.savetxt(fn, np.array([eps[i],R[i]]).T)
        pass

    ## save 0~2 pct to files (for fitting tau0s)
    mxi = []
    for i in range(len(data)):
        ## find the index
        for j in range(len(eps[i])):
            if eps[i][j]>0.02:
                mxi.append(j)
                break
            pass
        fn = 'str_%s_pre.txt'%str(i).zfill(2)
        print fn
        np.savetxt(fn, np.array([eps[i][:j],sig[i][:j]]).T)
        fn = 'R_%s_pre.txt'%str(i).zfill(2)
        print fn
        np.savetxt(fn, np.array([eps[i][:j],R[i][:j]]).T)
        pass
    
    ## plot
    figs = []; axr=[]
    figr = []; axs=[]
    for i in range(len(data)):
        figs.append(plt.figure(ifig+i*2));   figs[i].clf()
        figr.append(plt.figure(ifig+i*2+1)); figr[i].clf()
        axs.append(figs[i].add_subplot(111))
        axr.append(figr[i].add_subplot(111))
        axs[i].plot(eps[i][::jump], sig[i][::jump], marker='x', ls='None')
        axr[i].plot(eps[i][::jump],   R[i][::jump], marker='+', ls='None')
        pass
    pass
