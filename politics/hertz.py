"""module has spectrum information:

BANDS is a dictionary of spectrum standards

nu2band(f) is a function converting a frequency to band(s)

"""
## \namespace geo.politics.hertz The EM Spectrum Standard
try:
    from numpy import inf
except ImportError:
    from sys import float_info
    inf = float_info.max

from .interval import Interval as I
    
## List of Bands, from 0Hz to MeV (http://en.wikipedia.org/wiki/Radio_spectrum ):\n __NATO__: A, B, C, D, E, F, G, H, I, J, K, L (Lx), M\n __ITU__: ELF, SLF, ULF, VLF, LF, MF, VHF, UHF, SHF, EHF\n __IEEE__: HF, VHF, UHF, L, S, C, X, Ku, K, Ka, Q, V, W\n __RSGB__: L, S, C, X, Ku, K, Ka, Q, U, V, E, W, F, C\n __CIE__: FIR, LWIR, MWIR, SWIT, NIR\n _Visible_: red, orange, yellow, green, blue\n UV: UVA, (NUV), UVB, (MUV), UVC, (FUV), VUV, LUV, SUV, EUV, X-rays, and \f$\gamma\f$ rays.
BANDS = dict(misc=dict(subHz=I(0., 3.),
                       THz=I(300.e9, 3000.e9)),
             NATO=dict(A=I(3., 250e6),
                       B=I(250.e6, 500.e6),
                       C=I(500.e6, 1.e9),
                       D=I(1e9, 2e9),
                       E=I(2e9, 3e9),
                       F=I(3e9, 4e9),
                       G=I(4e9, 6e9),
                       H=I(6e9, 8e9),
                       I=I(8e9, 10e9),
                       J=I(10e9, 20e9),
                       K=I(20e9, 40e9),
                       L=I(40e9, 60e9),
                       M=I(60e9, 100e9)),
             ITU=dict(ELF=I(3, 30.),
                      SLF=I(30., 300.),
                      ULF=I(300., 3000.),
                      VLF=I(3.e3, 30.e3),
                      LF=I(30.e3, .300e3),
                      MF=I(300.e3, 3000.e3),
                      HF=I(3.e6, 30.e6),
                      VHF=I(30.e6, 300.e6),
                      UHF=I(0.3e9, 3.e9),
                      SHF=I(3.e9, 30.e9),
                      EHF=I(30.e9, 300.e9)),
             IEEE=dict(HF=I(3.e6, 30.e6),
                       VHF=I(30.e6, 300.e6),
                       UHF=I(0.3e9, 3.e9),
                       L=I(1.e9, 2.e9),
                       S=I(2.e9, 4.e9),
                       C=I(4.e9, 8.e9),
                       X=I(8.e9, 12.e9),
                       Ku=I(12.e9, 18.e9),
                       K=I(18.e9, 27.e9),
                       Ka=I(26.5e9, 40.e9),
                       Q=I(33.e9, 50.e9),
                       V=I(50.e9, 75.e9),
                       W=I(75.e9, 110.e9),
                       mm=I(110.e9, 300.e9)),
             RSGB=dict(L=I(1.e9, 2.e9),
                       S=I(2.e9, 4.e9),
                       C=I(4.e9, 8.e9),
                       X=I(8.e9, 12.e9),
                       Ku=I(12.e9, 18.e9),
                       K=I(18.e9, 27.e9),
                       Ka=I(26.5e9, 40.e9),
                       Q=I(30.e9, 50.e9),
                       U=I(40.e9, 60.e9),
                       V=I(50.e9, 75.e9),
                       E=I(60.e9, 90.e9),
                       W=I(75.e9, 110.e9),
                       F=I(90.e9, 140.e9),
                       D=I(110.e9, 170.e9)),
             CIE=dict(FIR=I(300.e9, 20.e12),
                      LWIR=I(20.e12, 37.474e12),
                      MWIR=I(37.474e12, 100.e12),
                      SWIR=I(100.e12, 214.13747e12),
                      NIR=I(214.13747e12, 400.e12)),
             visible=dict(red=I(400.e12, 484.e12),
                          orange=I(484.e12, 508.5e12),
                          yellow=I(508.e12, 526.e12),
                          green=I(526.e12, 606.e12),
                          blue=I(606.8e12, 668.e12),
                          violet=I(668.e12, 789.e12)),
             UVx=dict(UVA=(750.e12, 937.e12),
                      UVB=(937.e12, 1071.e12),
                      UVC=(1071.e12, 2998.e12)),
             UV=dict(UVA=I(750.e12, 937.e12),
                     NUV=I(750.e12, 999.3e12),
                     UVB=I(937.e12, 1071.e12),
                     MUV=I(999.3e12, 1499.e12),
                     UVC=I(1071.e12, 2998.e12),
                     FUV=I(1499.e12, 2457.e12),
                     VUV=I(1499.e12, 30.e15),
                     LUV=I(2998.e12, 3406.7e12),
                     SUV=I(1998.6e12, 30.e15),
                     EUV=I(2478.e12, 30.e15)),
             Xray=dict(soft=I(30.e15, 3000.e15),
                       hard=I(3.e18, 30.e18)),
             GammaRay=dict(PhotoElectricEffect=I(1.e19, 7.6e19*1.5),
                           ComptonScattering=I(7.6e19*1.5, 1.5e22),
                           PairProduction=I(1.5e21,  1.519e43),
                           PlankScale=I(1.519e43, inf)))


## Dictionary of \f$L_n\f$ bands \n
# ( http://en.wikipedia.org/wiki/GPS_signals#Frequencies_used_by_GPS )
LBAND = {1: 1575.42,
         2: 1227.60,
         3: 1381.05,
         4: 1379.913,
         5: 1176.45}

_deltaL = 1712000.0 / 2

         
BANDS['IEEE']['Lx'] = {key:I(value-_deltaL, value+_deltaL) for key, value in
                       dict(LA=1452.960e6,
                            LB=1454.672e6,
                            LC=1456.384e6,
                            LD=1458.096e6,
                            LE=1459.808e6,
                            LF=1461.520e6,
                            LG=1463.232e6,
                            LH=1464.944e6,
                            LI=1466.656e6,
                            LJ=1468.368e6,
                            LK=1470.080e6,
                            LL=1471.792e6,
                            LM=1473.504e6,
                            LN=1475.216e6,
                            LO=1476.928e6,
                            LP=1478.640e6,
                            LQ=1480.352e6,
                            LR=1482.064e6,
                            LS=1483.776e6,
                            LT=1485.488e6,
                            LU=1487.200e6,
                            LV=1488.912e6,
                            LW=1490.624e6).iteritems()}





## Convert Frequencty to Band Name(s)
# \param f Frequency (in Hz)
# \returns dict Of names of bands
def nu2band(f):
    """dict of bands = nu2band(f) """
    from itertools import ifilter
    # not pretty- there must be a better way.
    return {item[0][0]:item[0][1] for item in
            ifilter(bool,
                    [[(s, k) for k, v in b.items() if v[0] <= f < v[1]]
                     for s, b in BANDS.items()])}



