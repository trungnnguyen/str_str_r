"""
in-plane variations of selected physical entities
"""
def inpv(
    path=None, ext='*.csv', ifig=92, order=3,
    epsl=0.05, epsu=0.10, delt=10, mod=None,
    mode='NIST',
    ):
    """
    Post-process the in-plane variations of selected physical entities
    Post-processes a series of in-plane tensile tests.
    Files from glob.glob('*.csv') are under process.
    The file name should be 'a_b_c.csv' with a, b and c are 
    referring to material label, angle in 3 digits and
    duplicate number, respectively.

    Basically, it runs nist_column_post_process on the files
    whose extension is matched.

    Arguments
    =========
      path  = None (a name of directory containing the files)
      ext   = '*.csv'
      ifig  = 92
      order = 3
      epsl  = 0.05
      epsu  = 0.10
      delt  = 10
      mod   = None
      mode  = 'NIST', or 'MTS'

    Variables
    =========
      files             : glob.glob(ext)
      angles            : angles are obtained from the filename...
      dang              : angles[1] - angles[0]
      ang_off           : dang / 5.
      FilesForthisAngle : Contains Files of the current angle
      r15accl           : R-value array list
      ys1               : ys_Teps, YS at the offset using total strain
      ys2               : ys_Peps, YS at the offest using plastic strain
      r_lowel           : R value at the lower strain
      r_upel            : R value at the upper strain
      siglowel          : Flow stress at the lower strain
      sigupel           : Flow stress at the upper strain
      r_lowel           : R value at the lower strain
      r_upel            : R value at the upper strain

      ys1               : list containing ys_Teps
      ys2               : list containning ys_Peps
      ys1nl             : list containing ys_Teps
      ys1n              : mean of ys1nl
      ys2nl             : list containing ys_Peps
      ys2n              : mean of ys2nl

      R15A              : mean of r15accl
      IR                : mean of InstRMeanl
      IRSLOP            : mean of InstRslopel
      YS1               : mean of ys1
      YS2               : mean of ys2
      SIGL              : mean of siglowel
      SIGU              : mean of sigupel
      R_LOWE            : mean of r_lowel

      * master variables as below
      ys1_mast          : list containing YS1
      ys2_mast          : list containing YS2
      ysu_mast          : list containing SIGU
      ysl_mast          : list containing SIGL

    Dependents
    ==========
      os, glob, matplotlib.pyplot as plt,
      nist_column_post_process
      mts_column_post_process

    Outputs
    =======
      Matplotlib figures
      figR15a : Accummulative R-value vs angles
      figIR   : Instantaneous R-value vs angles
           (red: InstRMeanl.mean()), 
      figIRS  : Instantaneous R-value / time vs angles to see
                if it decreases or increases.
      figYS1  : ys_Teps, ys_Peps, sig_lowe vs angle
      figYS2  : ys_Teps, ys_Peps, sig_upe  vs angle
      figYSn  : ys1_mast / ys1n with errorbar
    """
    import os, glob
    import matplotlib.pyplot as plt
    from strr import nist_column_post_process, mts_column_post_process

    # subject files whose extension is given as 'ext'
    if path==None: files = glob.glob(ext)
    elif type(path).__name__=='str':
        files = glob.glob(
            '%s%s%s'%(path,os.sep,ext))

    if len(files)==0: raise IOError, 'No files are available.'
    # sorting angles
    # In the below, it is assumed that the files have
    # 'a_b_c.csv' a: is label for material, b: 3-digit angle, c:duplicate
    print 'files as below'
    for i in files: print i
    print '\n'
    for i in range(len(files[-1])/2): print '-',
    print '\n'

    angles = []
    for i in range(len(files)):
        temp = int(files[i].split('_')[1])
        if temp in angles: pass
        else: angles.append(temp)

    angles.sort() # sorting angle

    #mapping... into integers.
    angles = map(int, angles)
    print 'angles as below:'
    print angles

    #figR15  = plt.figure(ifig)
    figR15a = plt.figure(ifig+1)
    figIR   = plt.figure(ifig+2)
    figIRS  = plt.figure(ifig+3)
    figYS1  = plt.figure(ifig+4)
    figYS2  = plt.figure(ifig+5)
    figYSn  = plt.figure(ifig+6)

    #figR15.clf()
    figR15a.clf()
    figIR.clf()
    figIRS.clf()
    figYS1.clf()
    figYS2.clf()
    figYSn.clf()

    #axR15  = figR15.add_subplot(111)
    axR15a = figR15a.add_subplot(111)
    axIR   = figIR.add_subplot(111)
    axIRS  = figIRS.add_subplot(111)
    axYS1  = figYS1.add_subplot(111)
    axYS2  = figYS2.add_subplot(111)
    axYSn  = figYSn.add_subplot(111)

    ## 
    dang = abs(angles[1] - angles[0])
    ang_off = dang / 5.

    ys1_mast = []; ys1_err = []
    ys2_mast = []; ys2_err = []
    ysl_mast = []; ysl_err = []
    ysu_mast = []; ysu_err = []
    ys1nl, ys2nl, yslnl, ysunl  = [], [], [], []

    ## angles include available angles 
    for i in range(len(angles)):
        print 'current angle: %i'%angles[i]
        FilesForThisAngle = []
        for j in range(len(files)):
            if int(files[j].split('_')[1])==angles[i]:
                FilesForThisAngle.append(files[j])
        print '     %i samples tested'%len(
            FilesForThisAngle)

        r15accl, InstRMeanl, InstRslopel = [], [], []
        ys1, ys2 = [], []
        r_lowel, r_upel = [], []
        siglowel, sigupel = [], []

        ## Yield stress to be used for
        # Normalization w.r.t result along 0
        for j in range(len(FilesForThisAngle)):
            ## getting one value parameters from nist_column_post_process.
            if mode.upper()=='NIST':
                rst = nist_column_post_process(
                    datafile=FilesForThisAngle[j], 
                    echo=False, delt=delt, modulus=mod)
            elif mode.upper()=='MTS':
                rst = mts_column_post_process(
                    datafile=FilesForThisAngle[j],
                    #t,L0,W0,area,itime, 
                    delt=delt, mod=mod,
                    )
            else: raise IOError

            # r15,
            sig_lowe, sig_upe, r_lowe, r_upe = rst[0:4]
            r15acc, InstRMean, InstRSTDV, InstRslope = rst[4:8]
            ys_Teps, ys_Peps = rst[8:] #Yield stress

            # axR15.plot(angles[i],r15, 'o', mfc='None', mec='black')
            axR15a.plot(
                angles[i], r15acc,'o', mfc='None', mec='black')
            axIR.errorbar(
                angles[i],InstRMean,yerr=InstRSTDV,color='grey',
                fmt='o', mfc='None', mec='black', ls='--')
            axIR.plot(angles[i] - ang_off, r_lowe,
                      'o', mfc='None', mec='blue')
            axIR.plot(angles[i] + ang_off, r_upe,
                      'o', mfc='None', mec='green')
            # axIR.plot(angles[i],InstRMean,
            #  'o', mfc='None', mec='black')
            axIRS.plot(angles[i],InstRslope,
                       'o', mfc='None', mec='black')
            # yield stress(0.2pct offset from Total strain),
            # flow stress at lowe strain(blue)
            # and flow stress at upe strain(green).
            if i==0 and j==0:
                axYS1.plot(
                    angles[i] - ang_off,ys_Teps,
                    'o', mfc='None', mec='black',
                    label=r'$\sigma^{YS}_{\varpepsilon_{tot}}$')
                axYS1.plot(
                    angles[i], sig_lowe,
                    'o', mfc='None', mec='blue',
                    label=r'$\sigma^{low}_{\varepsilon_{tot}}$')
                axYS1.plot(
                    angles[i] + ang_off,sig_upe,
                    'o', mfc='None', mec='green')
                axYS2.plot(angles[i] - ang_off,ys_Teps,
                           'o', mfc='None', mec='black')            
                axYS2.plot(angles[i], ys_Peps,
                           'o', mfc='None', mec='black')
                axYS2.plot(angles[i] + ang_off,sig_upe,
                           'o', mfc='None', mec='green')
            else:
                axYS1.plot(angles[i] - ang_off,ys_Teps,
                           'o', mfc='None', mec='black')
                axYS1.plot(angles[i], sig_lowe,
                           'o', mfc='None', mec='blue')
                axYS1.plot(angles[i] + ang_off,sig_upe,
                           'o', mfc='None', mec='green')
                axYS2.plot(angles[i] - ang_off,ys_Teps,
                           'o', mfc='None', mec='black')            
                axYS2.plot(angles[i], ys_Peps,
                           'o', mfc='None', mec='black')
                axYS2.plot(angles[i] + ang_off,sig_upe,
                           'o', mfc='None', mec='green')            
            ## list
            # r15l.append(r15)
            r15accl.append(r15acc)
            InstRMeanl.append(InstRMean)
            InstRslopel.append(InstRslope)
            ys1.append(ys_Teps)
            ys2.append(ys_Peps)
            if abs(angles[i]-0.)<0.01: #
                ys1nl.append(ys_Teps)
                ys2nl.append(ys_Peps)
                yslnl.append(sig_lowe)
                ysunl.append(sig_upe)
                pass

            # yield stress and r value at the low
            # and up prescribed strains
            siglowel.append(sig_lowe)
            sigupel.append(sig_upe)
            r_lowel.append(r_lowe)
            r_upel.append(r_upe)
            pass # each sample for an angle

        # mean values for the current angle: angle[i]
        #R15   = np.array(r15l).mean()
        R15A   = np.array(r15accl).mean()
        IR     = np.array(InstRMeanl).mean()
        IRSLOP = np.array(InstRslopel).mean()
        YS1    = np.array(ys1).mean()
        YS2    = np.array(ys2).mean()
        SIGL   = np.array(siglowel).mean()
        SIGU   = np.array(sigupel).mean()
        R_LOWE = np.array(r_lowel).mean()
        R_UPE  = np.array(r_upel).mean()

        # master yield stresses
        ys1_mast.append(YS1)
        ys2_mast.append(YS2)
        ysu_mast.append(SIGU)
        ysl_mast.append(SIGL)

        # standard deviation
        if len(FilesForThisAngle)>1 :
            #axR15.errorbar(
            # angles[i],R15,fmt='+', yerr=np.std(r15l),
            #  color='red')
            """
            index i on every angle
            index j on files on i-th angle
            """
            label=None
            axR15a.errorbar(
                angles[i],R15A,fmt='+', yerr=np.std(r15accl),
                color='red')

            ## Instantaneous R-value at -----------------------------#
            if i==0 and j==0:
                axIRp1, = axIR.errorbar(
                    angles[i], IR, fmt='+', yerr=np.std(InstRMeanl),
                    color='red', label=r'$\bar{R^{inst}}$')
                axIRp2, = axIR.errorbar(
                    angles[i]-ang_off,R_LOWE,fmt='+',
                    yerr=np.std(r_lowel), color='blue',
                    label=r'$R^{inst}_{low}$')
                axIRp3, = axIR.errorbar(
                    angles[i]+ang_off,R_UPE, fmt='+',
                    yerr=np.std(r_upel), color='green',
                    label=r'$R^{inst}_{up}$')
                axIR.legend()
            else:
                axIR.errorbar(
                    angles[i], IR, fmt='+', yerr=np.std(InstRMeanl),
                    color='red', label=label)
                axIR.errorbar(
                    angles[i]-ang_off,R_LOWE,fmt='+',
                    yerr=np.std(r_lowel), color='blue')
                axIR.errorbar(
                    angles[i]+ang_off,R_UPE, fmt='+',
                    yerr=np.std(r_upel), color='green')
            ## ------------------------------------------------------#

            axIRS.errorbar(
                angles[i],IRSLOP,fmt='+',yerr=np.std(InstRslopel),
                color='red')
            axYS1.errorbar(
                angles[i],YS1,fmt='+',yerr=np.std(ys1),
                color='red')
            axYS2.errorbar(
                angles[i],YS2,fmt='+', yerr=np.std(ys2),
                color='red')

            # save standard deviation to the master error lists
            # for yield stress
            ys1_err.append(np.std(ys1))
            ys2_err.append(np.std(ys2))
            ysl_err.append(np.std(siglowel))
            ysu_err.append(np.std(sigupel))
        else:
            #axR15.plot( angles[i],R15,'+',    color='red')
            axR15a.plot( angles[i],R15A,'+',   color='red')

            # Instantaneous R-value -----------------------------#
            axIR.plot( angles[i],IR,'+',     color='red')
            axIR.plot( angles[i]-ang_off, R_LOWE, color='blue')
            axIR.plot( angles[i]+ang_off, R_UPE,  color='green')
            #----------------------------------------------------#

            axIRS.plot( angles[i],IRSLOP,'+', color='red')
            axYS1.plot( angles[i],YS1,'+',    color='red')
            axYS2.plot( angles[i],YS2,'+',    color='red')
            ys1_err.append(np.NAN)
            ys2_err.append(np.NAN)            
            ysl_err.append(np.NAN)
            ysu_err.append(np.NAN)

    ## ticks
    d = abs(angles[1]-angles[0])
    ang0 = min(angles) - d
    if max(angles)>=90.:# -->more than 90 degree
        ang1 = max(angles) + 2 * d
    else: ang1 = max(angles) + d
    if max(angles)>=090.: fin = max(angles) + d
    else: fin = max(angles)    

    ## normalized in-plane yield stress variation.
    ys1n = np.array(ys1nl).mean()
    ys2n = np.array(ys2nl).mean()
    ysln = np.array(yslnl).mean()
    ysun = np.array(ysunl).mean()

    angles = np.array(angles) # numpy the angles list

    axYSn.errorbar(
        angles - 1.5*d/8., ys1_mast/ys1n,
        yerr=ys1_err/ys1n, fmt='o', label=r'$\sigma^{YS}_1$')
    axYSn.errorbar(
        angles - 0.5*d/8., ys2_mast/ys2n,
        yerr=ys2_err/ys2n, fmt='o', label=r'$\sigma^{YS}_2$')
    axYSn.errorbar(
        angles + 0.5*d/8., ysl_mast/ysln,
        yerr=ysl_err/ysln, fmt='o', label=r'$\bar{\sigma^{YS}}_l$')
    axYSn.errorbar(
        angles + 1.5*d/8., ysu_mast/ysun,
        yerr=ysu_err/ysun, fmt='o', label=r'$\bar{\sigma^{YS}}_u$')

    leg = axYSn.legend(loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    axYSn.set_xlim(ang0, ang1)
    axYSn.set_xticks(np.arange(0, fin*1.01, 30.))
    figYSn.savefig('YS_norm.pdf')
    axYSn.set_ylim(0.,)
    figYSn.savefig('YS_norm_rescale.pdf')

    #axR15.set_xlim(ang0, ang1)
    axR15a.set_xlim(ang0, ang1)
    axIR.set_xlim(ang0,   ang1)
    axIRS.set_xlim(ang0,  ang1)
    axYS1.set_xlim(ang0,  ang1)
    axYS2.set_xlim(ang0,  ang1)

    #axR15.set_xticks(np.arange(0,max(angles)*1.01,30.))
    axR15a.set_xticks(np.arange(0,fin*1.01, 30.))
    axIR.set_xticks(np.arange(0,  fin*1.01, 30.))
    axIRS.set_xticks(np.arange(0, fin*1.01, 30.))
    axYS1.set_xticks(np.arange(0, fin*1.01, 30.))
    axYS2.set_xticks(np.arange(0, fin*1.01, 30.))

    #axR15.set_ylabel('R15')
    axR15a.set_ylabel('R15a', dict(fontsize=20))
    axIR.set_ylabel(r'$\bar{R}^{inst}$', dict(fontsize=20))
    axIRS.set_ylabel(r'$\dot{R}^{inst}$', dict(fontsize=20))
    axYS1.set_ylabel(r'$\sigma^{YS}$ offset from $\varepsilon^{tot}$',
                     dict(fontsize=20))
    axYS2.set_ylabel(r'$\sigma^{YS}$ offset from $\varepsilon^{pl}$',
                     dict(fontsize=20))

    #axR15.set_xlabel(r'$\theta$')
    axR15a.set_xlabel(r'$\theta$', dict(fontsize=20))
    axIR.set_xlabel(r'$\theta$', dict(fontsize=20))
    axYS1.set_xlabel(r'$\theta$', dict(fontsize=20))
    axYS2.set_xlabel(r'$\theta$', dict(fontsize=20))

    #figR15.savefig('figR15.pdf')
    figR15a.savefig('figR15a.pdf')
    figIR.savefig('figIR.pdf')
    axIR.set_ylim(0.,); figIR.savefig('figIR_rescale.pdf')

    figIRS.savefig('figIRS.pdf')
    figYS1.savefig('figYS1.pdf')
    axYS1.set_ylim(0,)
    figYS1.savefig('figYS1_rescale.pdf')
    figYS2.savefig('figYS2.pdf')
    axYS2.set_ylim(0.,)
    figYS2.savefig('figYS2_rescale.pdf')
    pass
    # sorting multiple tensile bars for an angle

