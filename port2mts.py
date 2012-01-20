"""
Make a series of raw data file into 'mts' mode compatible files.
"""
def __main__(
    filename=None, nheader=None, label=None, angle=None,
    itime=None, iext=None, iwid=None, ifor=None, 
    iforce_stress=None, # 'force'
    iforce_unit=None, # kN
    axial_strain_or_extension=None, #'extension'
    width_strain_or_extension=None, #'extension'
    T0=None, L0=None, W0=None, delimiter=None, fmt='%12.6e'):
    """
    Make a single file into mts mode compatible
    Output file name has the extension of 'mts'.
    time, force, width extension, t extension,

    when itime==-1: time not available
    when iwid ==-1: wid not available
         (make width change column as zero-array)
    """
    import numpy as np
    import os

    ## argument check-ups
    if filename==None: raise IOError, 'File name is missing'

    # if anyting is missing show the upper part of the file.
    if None in [delimiter, nheader, angle, label,
                iext, iwid, iforce_stress, itime]:
        print '\n\n##############################################'
        with open(filename) as myfile:
            for i in range(8):
                print myfile.next(),
        print '##############################################\n\n'

    if delimiter==None: delimiter = raw_input('\n\ndelimiter? >>')
    if nheader==None: nheader = __nheaderFinder__(filename, delimiter)
    if angle==None: angle = int(raw_input('Angle of this file >>'))
    if label==None: label = raw_input('Label of this file >>')
    if itime==None: itime = int(raw_input('itime >> '))
    if iext==None: iext = int(raw_input('iext >> '))
    if iwid==None: iwid = int(raw_input('iwid >> '))
    if ifor==None: ifor = int(raw_input('ifor >> '))
    if iforce_stress==None: iforce_stress = raw_input('force or stress >>')
    if iforce_unit==None: iforce_unit = raw_input('force unit [kN, N] .. >>')
    if axial_strain_or_extension==None: axial_strain_or_extension = raw_input(
        "Axial 'extension' or 'strain' >>")
    if width_strain_or_extension==None: width_strain_or_extension = raw_input(
        "Width 'extension' or 'strain' >>")    
    if T0==None: T0 = float(raw_input('Thickness? >>'))
    if L0==None: L0 = float(raw_input('Gauge length? >>'))
    if W0==None: W0 = float(raw_input('Sample width? >>'))
    ##    

    time, axialextension, changeinwidth, axialforce = column_info(
        ### 
        T0, L0, W0,
        #
        itime, iext, iwid, ifor,
        #
        iforce_stress, iforce_unit, nheader,
        #
        filename, delimiter,
        #
        axial_strain_or_extension,
        width_strain_or_extension,
        )

    idp = 0
    while True:
        fout = '%s_%s_%s.mtsform'%(label, str(angle).zfill(3), idp)
        if not(os.path.isfile(fout)): break
        idp = idp + 1

    fout = open(fout, 'w')
    fileformat = '%s %s %s %s\n'%(fmt, fmt, fmt, fmt)
    fout.write('%f %f %f\n'%(W0, T0, L0)) ## width, thickness, gauge length
    for i in range(len(time)):
        fout.write(fileformat%(
                time[i], axialextension[i],
                changeinwidth[i], axialforce[i]))

def column_info(
    T0, L0, W0,
    itime, iext, iwid, ifor,
    iforce_stress, iforce_unit, nheader, 
    filename, delimiter,
    axial_strain_or_extension,
    width_strain_or_extension,
    ):
    """
    Change engineering stress to force(kN)
    Make all columns to the fullest.
    engineering strain -> extension[mm]
    engineering stress -> force[kN]
    engineering width strain -> change in width direction (extension [mm])
    """
    import numpy as np
    columns = np.loadtxt(filename, skiprows=nheader, delimiter=delimiter).T
    ndat = len(columns[0])
    # time
    if itime==-1:
        acqrate = float(raw_input(
                'Data acquisition rate in [Hz] >>'))
        time = []
        t = 0.
        for i in range(ndat):
            time.append(t)
            t = t + 1./acqrate
        time = np.array(time)
    else: time = columns[itime]

    # axial extension
    ext = columns[iext]
    if axial_strain_or_extension.upper()=='STRAIN':
        ext = strain_to_ext(ext, L0=L0)
    elif axial_strain_or_extension.upper()=='EXTENSION': pass
    else: raise IOError

    # width strain
    if iwid==-1:
        w = []
        for i in range(ndat): w.append(np.nan)
        wid = np.array(w)
    else:
        wid = columns[iwid]
        if width_strain_or_extension.upper()=='EXTENSION':pass
        elif width_strain_or_extension.upper()=='STRAIN':
            wid = strain_to_ext(wid, W0)
        else: raise IOError

    # force
    force = columns[ifor]
    if iforce_stress.upper()=='FORCE':
        if iforce_unit.upper()=='KN': coef = 1.0
        else:
            print 'The input force unit is not kN instead is %s'%(
                iforce_unit)
            coef = raw_input(
                'Type a coeffcient to be multiplied to be in [kN]')
    elif iforce_stress.upper()=='STRESS':
        force = stress_to_force(force, T0, W0, stress_unit)
    return time, ext, wid, force

def strain_to_ext(strain, L0=None):
    """ L0 should be in [mm] unit. """
    
    import numpy as np
    if type(strain).__name__=='list':
        strain = np.array(strain)
    return strain * L0

def stress_to_force(force, T0, W0, stress_unit):
    """
    T0 and W0 are in [mm]
    """
    import numpy as np
    if type(force).__name__=='list':
        force = np.array(force)
        pass
    # if stress_unit == 'Mpa'
    """
    1. MPa = 10^6          N / m^2
           = 10^6 x 10^-3 kN / (10^3mm)^2
           = 10^3 x 10^-6 kN / mm^2
           = 10^-3        kN/mm^2

    2. MPa = 10^-3        kN/mm^2
           = 10^-3 x 10^3  N/mm^2
           = 1.            N/mm^2
    
    3. MPa = 10^6 Pa = 10^-3 kN/mm^2
           -> Pa = 10^-9 kN/mm^2
    """
    if stress_unit.upper()=='MPA': force = force * (10.**-3)
    elif stress_unit.upper()=='PA': force = force * (10.**-9)
    else: raise IOError
    return force * T0 * W0

def __nheaderFinder__(filename, delimiter):
    f = open(filename, 'r')
    lines = f.readlines()
    for i in range(len(lines)):
        try: map(float, lines[i].split(delimiter))
        except: pass
        else:
            nheader = i
            break

    f.close()
    return nheader
    


def ex(ext='fss/*.txt'):
    from glob import glob
    files = glob(ext)
    for i in range(len(files)):
        __main__(filenames=files[i])
