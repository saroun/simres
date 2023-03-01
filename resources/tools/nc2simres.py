import NCrystal as nc
from scipy.special import bernoulli
import scipy.integrate as integ
import datetime
import numpy as np

# constatnts
_kb = 1.380658e-23 # Boltzman const. [J/K]
_kb_meV = 0.0861739 # Boltzman const. [meV/K]
_hsqov2m = 81.8042531 # h^2/2m_n [meV.A^2]
_h = 6.6260755e-34 # Planck const. [J.s]
_hovm = 3.9560346e3 # h/m_n [A.m/s]
_u = 1.6605402e-27 # mass unit [kg]
_csph = _h**2/_u/_kb*1e20 # h^2/u/kb in [Ang^2.K]

def _ksi(a):
    """Calculate integral of x*coth(x) for x=[0, a]."""
    _alim = 0.05
    def B0(a):
        # approx. for small a
        c = np.array([1.0, 1.0/9, 1.0/225, 2.0/6615])
        r = a
        for i in [1,2,3]:
            r += c[i]*a**(2*i+1)
        return r
    
    def xcoth(x):
        ep = np.exp(x)
        em = np.exp(-x)
        tot = (ep+em)/(ep-em)
        return tot*x  
    if (a<=_alim):
        res = B0(a)
    else:
        res0 = B0(_alim)
        res1 = integ.quad(xcoth, _alim, a)[0]
        res = res0 + res1
    return res

def _BT(T, TD, A):
        """Calculate temperature factor BT in [A^2]."""
        tt = T/TD
        # (J.s)^2 / kg /  (J/K)/K = J/kg.s^2 = m^2
        res = 3*_csph/A*tt**2/TD*_ksi(0.5/tt)
        # res = 3*_h**2/(A*_u*_kb)*tt**2/TD*_ksi(0.5/tt)*1e20 # to A^2
        return res 


def cal_RTDS(x):
    """Calculate the value of R(x=TD/T) function for single-phonon scattering.
    
    See e.g. A. Freund, Nucl. Instr. Meth. 213 (1983) 495.
    """
    b = bernoulli(30)
    if x <= 6.0:
        R = 0.0
        Ifact = 1.0
        Xn = 1.0/x
        for i in range(len(b)):
            R = R + b[i]*Xn/(Ifact*(i+2.5))
            Xn = Xn*x
            Ifact = Ifact*(i+1)
    else:
        R = 3.3/np.sqrt(x**7)
    return R


def dumpInfo(info, fname):
    """Print selected values fromn NCrystal Info object."""
    print("Info on {}".format(fname))
    if info.hasComposition():
        print("Composition:")
        cc = info.getComposition()
        for c in cc:
            #print(c[1].__dict__)
            print("\t{:g}, {}".format(* c))
            #print(c[1].captureXS())


    if info.hasTemperature():
        print("Temperature: {:g}".format(info.getTemperature()))

    if info.hasDensity():
        print("Density: {:g}".format(info.getDensity()))

    if info.hasAtomInfo():
        print("Atoms:")
        aa = info.getAtomInfo()
        for a in aa: 
            print("\t{}".format(a))

    if info.hasStructureInfo():
        print("Structure:")
        aa = info.getStructureInfo()
        print("\t{}".format(aa))

    print("Compound parameters:")
    print("SIGA={:g}".format(info.getXSectAbsorption()))
    print("SIGF={:g}".format(info.getXSectFree()))
    print("n_density={:g}".format(info.getNumberDensity()))
    hkl = info.hklList()
    for h in hkl:
        if (h[4]>0.5):
            print(h)


def genHKLlist(info, BT, dmin=0.0):
    """Generate list of reflections for SIMRES table.

    The list contains 6 columns:

    dhkl[A], multiplicity, Fhkl/cell_volume [fm/A^3], h, k, l 

    Fhkl excludes DW factor, which can be calculated from provided BT parameter.
    Fhkl is obtained from NCrystal hkl list, which contains F2 = Fhkl**2 including DW factor.
    Therefore, Fhkl in the table is calculated as sqrt(F2)*exp(BT/dhkl**2).

    """
    hkl = info.hklList()
    ss = info.getStructureInfo()
    volc = ss['volume']
    table = []
    for h in hkl:
        if (h[4]>=dmin):
            Fhkl = 10*np.sqrt(h[5])*np.exp(BT/h[4]**2)/volc
            row = [h[4], h[3], Fhkl, abs(h[0]), abs(h[1]), abs(h[2])]
            table.append(row)
    return table


def genData(inpfile, temperature=None, dmin=0.0):
    info = nc.createInfo(inpfile)
    if temperature:
        T = temperature
    else:
        T = info.getTemperature()
    print("Generating data for {}, T={:g}".format(inpfile, T))
    data = {}
    data['input'] = inpfile
    if not info.hasComposition():
        raise Exception('{}: No composition info available.'.format(inpfile))
    if not info.hasAtomInfo():
        raise Exception('{}: No atom info available.'.format(inpfile))
    if not info.hasStructureInfo():
        raise Exception('{}: No structure info available.'.format(inpfile))
    SIGA = 0.0
    SIGI = 0.0
    SIGC = 0.0
    SIGF = 0.0
    SIGS = 0.0
    BT = 0.0
    BMPH = 0.0
    suma = 0.0
    cc = info.getComposition()
    aa = info.getAtomInfo()
    ss = info.getStructureInfo()
    #print(ss)
    #print(aa)
    volc = ss['volume']
    n_at = ss['n_atoms']
    ncomp = len(cc)
    if (ncomp != len(aa)):
        raise Exception("Unequalt number of components and atoms in {}".format(inpfile))
    for i in range(ncomp):
        comp = cc[i]
        atom = aa[i]
        p = comp[0]
        suma += p
        c = comp[1]
        A = c.averageMassAMU()
        SIGA += p*c.captureXS()
        SIGI += p*c.incoherentXS()
        SIGC += p*c.coherentXS()
        stot= c.incoherentXS() + c.coherentXS()
        SIGF += p*stot*(A/(A+1))**2
        TD = atom.debyeTemperature
        MSD = atom.meanSquaredDisplacement
        R = cal_RTDS(TD/T)
        SIGS += p*3*stot/A*np.sqrt(_kb_meV*TD/_hsqov2m)*R
        BT_part = _BT(T, TD, A)
        #print("BT={:g}, T={:g}, TD={:g}, A={:g}".format(BT_part,T, TD, A))
        BT += p*BT_part
        C2 = 4.27*np.exp(A/61)
        BMPH += p*4*BT_part*C2*_hsqov2m*1e-3
        #print('bc={:g}'.format(c.coherentScatLenFM()))
    data = {}
    params = {}
    # note: convert cross-sections to 1/mm
    params['SIGA'] = 0.1*SIGA/suma*n_at/volc/1.798 # convert to 1 Ang
    params['SIGI'] = 0.1*SIGI/suma*n_at/volc
    params['SIGC'] = 0.1*SIGC/suma*n_at/volc
    params['SIGF'] = 0.1*SIGF/suma*n_at/volc
    params['SIGS'] = 0.1*SIGS/suma*n_at/volc
    params['BT'] = BT/suma
    params['BMPH'] = BMPH/suma
    data['params'] = params
    data['hkl'] = genHKLlist(info, params['BT'], dmin=dmin)
    data['source'] = inpfile
    return data


def mergeData(datalist, names=None):
    result = {}
    params = {}
    # note: convert cross-sections to 1/mm
    params['SIGA'] = 0.0
    params['SIGI'] = 0.0
    params['SIGC'] = 0.0
    params['SIGF'] = 0.0
    params['SIGS'] = 0.0
    params['BT'] = 0.0
    params['BMPH'] = 0.0
    result['params'] = params
    result['hkl'] = []
    result['source'] = []

    # labels for phases
    if names and len(names)==len(datalist):
        lbls = names
    else:
        lbls = []
        for i in range(len(datalist)):
            lbl = 'phase{:d}'.format(i+1)
            lbls.append(lbl)
    for i in range(len(datalist)):
        d = datalist[i]
        f = d[0]
        dset = d[1]
        for p in params.keys():
            result['params'][p] += f*dset['params'][p] 
    # now we have the averaged Debye-Waller coefficient 
    BT_ave = result['params']['BT']
    for i in range(len(datalist)):
        d = datalist[i]
        f = d[0]
        dset = d[1]
        hkl = dset['hkl'].copy()
        # multiply Fhkl by sqrt(fraction)*exp(-W+W_ave)
        BT = dset['params']['BT']
        #print('BT={:g}, BT_ave={:g}'.format(BT, BT_ave))
        fmult = np.sqrt(f)
        for h in hkl:
            W_ave = BT_ave/h[0]**2
            W = BT/h[0]**2
            fac = fmult*np.exp(-W+W_ave)
            h[2] *= fmult*np.exp(-W+W_ave)
            # append phase label to each hkl
            h.append(lbls[i])
        result['hkl'].extend(hkl)
        result['source'].append([lbls[i], dset['source'], f*100])
    # sort relections by dhkl
    result['hkl'].sort(key = lambda x: x[0], reverse=True) 
    return result



def gen_SIMRES_table(inpfile, outfile, dmin=0.5, names=None):

    isList = False
    if isinstance(inpfile,list):
        datalist = []
        for inp in inpfile:
            f = inp[0]
            data = genData(inp[1], temperature=300.0, dmin=dmin)
            datalist.append([f,data])
        data = mergeData(datalist, names=names)
        nphases = len(inpfile)
        isList = True
    else:
        data = genData(inpfile, temperature=300.0, dmin=dmin)
        nphases = 1
    hdr = []
    hdr.extend(["polycrystal lookup table"])
    datestr = "Date: {}".format(str(datetime.datetime.now()))
    hdr.extend([datestr])
    hdr.extend(["Generated from NCrystal (https://github.com/mctools/ncrystal)"])
    hdr.extend(["Reference: X.-X. Cai and T. Kittelmann, Comp. Phys. Comm. 246 (2020) 106851"])
    if isList:
        for src in data['source']:
            hdr.extend(["Phase {}: {}, {:g}%".format(* src)])
        colhdr = "dhkl[A], multiplicity, FHKL/cell_volume [fm/A^3], h, k, l, phase"
    else:
        hdr.extend(["Source file: {}".format(inpfile)])
        colhdr = "dhkl[A], multiplicity, FHKL/cell_volume [fm/A^3], h, k, l"
    params = []
    pdata = data['params']
    for key in pdata.keys():
        params.extend(["{}={:g}".format(key, pdata[key])])
    f = open(outfile, "w")
    fmtc = "# {}\n"
    fmt = "{}\n"
    for h in hdr:
        f.write(fmtc.format(h))
    fmt = "{}\n"
    for p in params:
        f.write(fmt.format(p))
    f.write(fmtc.format(colhdr))
    fmtd = 2*"{:g}\t" + 1*"{:g}\t" + 2*"{:g}\t" + "{:g}"
    if isList:
        fmtd += "\t{}"
    fmtd += "\n"
    for hkl in data['hkl']:
        f.write(fmtd.format(* hkl))
    f.close()
    print("{} converted in {}".format(inpfile, outfile))


# fname = "Fe_sg229_Iron-alpha.ncmat"
#fname = "MgO_sg225_Periclase.ncmat"
#info = nc.createInfo(fname)    
#dumpInfo(info, fname)
#data = genData(fname, temperature=None, dmin=1)
#print(data)


#%% Convert selected data to SIMRES table format
inp = [[0.5,"Fe_sg229_Iron-alpha.ncmat"], [0.5,"Fe_sg225_Iron-gamma.ncmat"]]
gen_SIMRES_table(inp, "duplex_steel.dat", dmin=0.2, names=['alpha', 'gamma'])

gen_SIMRES_table("Fe_sg229_Iron-alpha.ncmat", "Fe_alpha.dat", dmin=0.2)
gen_SIMRES_table("Fe_sg225_Iron-gamma.ncmat", "Fe_gamma.dat", dmin=0.2)
gen_SIMRES_table("Cu_sg225.ncmat", "Cu.dat", dmin=0.2)
gen_SIMRES_table("Al_sg225.ncmat", "Al.dat", dmin=0.2)
gen_SIMRES_table("Mg_sg194.ncmat", "Mg.dat", dmin=0.2)
gen_SIMRES_table("Ni_sg225.ncmat", "Ni.dat", dmin=0.2)
gen_SIMRES_table("SiO2-alpha_sg154_AlphaQuartz.ncmat", "a-SiO2.dat", dmin=0.2)
gen_SIMRES_table("Ti_sg194.ncmat", "Ti.dat", dmin=0.2)

