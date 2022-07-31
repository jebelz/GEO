"""Ellipsoids do a lot. They provide a foundation for coordinate classes, they
hold transformation functions, affine computations, and distance computations,
latitude conversions,
"""
## \namespace geo.detic.newton.almanac
## <a href="http://en.wikipedia.org/wiki/Earth_ellipsoid">Ellipsoids of
## Yore</a>, and now.

__all__ = ('ALMANAC',)


_AIRY_1830 = """
Ellipsoid [AIRY 1830]
  Semi-Major Axis (Equatorial Radius)..[6377563.396]
  Semi-Minor Axis (Polar Radius).......[6356256.909]
  Flattening...........................[0.00334085067870327]
  Flattening Inverse...................[299.32496126649505]
  First Eccentricity...................[0.0816733743281685]
  First Eccentricity Squared...........[0.006670540074149084]
  Second Eccentricity..................[0.08194714751155018]
  Second Eccentricity Squared..........[0.006715334985279728]
"""


_AUSTRALIAN_1965 = """
Ellipsoid [AUSTRALIAN 1965]
  Semi-Major Axis (Equatorial Radius)..[6378160.0]
  Semi-Minor Axis (Polar Radius).......[6356774.719]
  Flattening...........................[0.003352891899858333]
  Flattening Inverse...................[298.249997276158]
  First Eccentricity...................[0.08182018036905428]
  First Eccentricity Squared...........[0.006694541915624534]
  Second Eccentricity..................[0.08209543749645103]
  Second Eccentricity Squared..........[0.006739660857733726]"""

_BESSEL_1841 = """
Ellipsoid [BESSEL 1841]
  Semi-Major Axis (Equatorial Radius)..[6377397.155]
  Semi-Minor Axis (Polar Radius).......[6356078.963]
  Flattening...........................[0.0033427731536659813]
  Flattening Inverse...................[299.15281535132334]
  First Eccentricity...................[0.0816968308747341]
  First Eccentricity Squared...........[0.006674372174974933]
  Second Eccentricity..................[0.08197084080074579]
  Second Eccentricity Squared..........[0.006719218741581313]"""

_BESSEL_1841_NAMIBIA = """
Ellipsoid [BESSEL 1841 NAMIBIA]
  Semi-Major Axis (Equatorial Radius)..[6377483.865]
  Semi-Minor Axis (Polar Radius).......[6356165.383]
  Flattening...........................[0.003342773176894559]
  Flattening Inverse...................[299.152813272542]
  First Eccentricity...................[0.0816968311581122]
  First Eccentricity Squared...........[0.006674372221277079]
  Second Eccentricity..................[0.08197084108698455]
  Second Eccentricity Squared..........[0.006719218788507778]"""

_CLARKE_1866 = """
Ellipsoid [CLARKE 1866]
  Semi-Major Axis (Equatorial Radius)..[6378206.4]
  Semi-Minor Axis (Polar Radius).......[6356583.8]
  Flattening...........................[0.0033900753039287908]
  Flattening Inverse...................[294.9786982138982]
  First Eccentricity...................[0.0822718542230039]
  First Eccentricity Squared...........[0.006768657997291184]
  Second Eccentricity..................[0.0825517107388772]
  Second Eccentricity Squared..........[0.006814784945915172]"""


_CLARKE_1880 = """
Ellipsoid [CLARKE 1880]
  Semi-Major Axis (Equatorial Radius)..[6378249.145]
  Semi-Minor Axis (Polar Radius).......[6356514.87]
  Flattening...........................[0.003407561308111843]
  Flattening Inverse...................[293.4650060791153]
  First Eccentricity...................[0.08248339919132311]
  First Eccentricity Squared...........[0.0068035111421552025]
  Second Eccentricity..................[0.08276542745958347]
  Second Eccentricity Squared..........[0.006850115982567657]"""


_EVEREST_1830_INDIA = """
Ellipsoid [EVEREST 1830 INDIA]
  Semi-Major Axis (Equatorial Radius)..[6377276.345]
  Semi-Minor Axis (Polar Radius).......[6356075.413]
  Flattening...........................[0.0033244493186534527]
  Flattening Inverse...................[300.8016980102568]
  First Eccentricity...................[0.08147298125166744]
  First Eccentricity Squared...........[0.0066378466740345775]
  Second Eccentricity..................[0.08174473748851525]
  Second Eccentricity Squared..........[0.0066822021070661935]"""

_EVEREST_1830_MALAYSIA = """
Ellipsoid [EVEREST 1830 MALAYSIA]
  Semi-Major Axis (Equatorial Radius)..[6377298.556]
  Semi-Minor Axis (Polar Radius).......[6356097.55]
  Flattening...........................[0.003324449343845343]
  Flattening Inverse...................[300.80169573085277]
  First Eccentricity...................[0.0814729815598449]
  First Eccentricity Squared...........[0.006637846724250804]
  Second Eccentricity..................[0.08174473779978614]
  Second Eccentricity Squared..........[0.0066822021579557725]"""


_EVEREST_1956_INDIA = """
Ellipsoid [EVEREST 1956 INDIA]
  Semi-Major Axis (Equatorial Radius)..[6377301.243]
  Semi-Minor Axis (Polar Radius).......[6356100.228]
  Flattening...........................[0.0033244493543833783]
  Flattening Inverse...................[300.8016947773539]
  First Eccentricity...................[0.08147298168875934]
  First Eccentricity Squared...........[0.006637846745256864]
  Second Eccentricity..................[0.08174473792999415]
  Second Eccentricity Squared..........[0.0066822021792435045]"""


_EVEREST_1964_MALAYSIA_SINGAPORE = """
Ellipsoid [EVEREST 1964 MALAYSIA & SINGAPORE]
  Semi-Major Axis (Equatorial Radius)..[6377304.063]
  Semi-Minor Axis (Polar Radius).......[6356103.039]
  Flattening...........................[0.003324449295589469]
  Flattening Inverse...................[300.8017000971244]
  First Eccentricity...................[0.08147298096952207]
  First Eccentricity Squared...........[0.006637846628060143]
  Second Eccentricity..................[0.08174473720353695]
  Second Eccentricity Squared..........[0.006682202060475285]"""


_EVEREST_1969_MALAYSIA = """
Ellipsoid [EVEREST 1969 MALAYSIA]
  Semi-Major Axis (Equatorial Radius)..[6377295.664]
  Semi-Minor Axis (Polar Radius).......[6356094.668]
  Flattening...........................[0.0033244492833663726]
  Flattening Inverse...................[300.8017012030905]
  First Eccentricity...................[0.08147298081999425]
  First Eccentricity Squared...........[0.006637846603695205]
  Second Eccentricity..................[0.08174473705250837]
  Second Eccentricity Squared..........[0.006682202035783637]"""


_EVEREST_PAKISTAN = """
Ellipsoid [EVEREST PAKISTAN]
  Semi-Major Axis (Equatorial Radius)..[6377309.613]
  Semi-Minor Axis (Polar Radius).......[6356109.571]
  Flattening...........................[0.003324292418982392]
  Flattening Inverse...................[300.81589522323446]
  First Eccentricity...................[0.08147106184331902]
  First Eccentricity Squared...........[0.006637533917877878]
  Second Eccentricity..................[0.08174279880970892]
  Second Eccentricity Squared..........[0.006681885157244453]"""


_FISHER_1960 = """
Ellipsoid [FISHER 1960]
  Semi-Major Axis (Equatorial Radius)..[6378155.0]
  Semi-Minor Axis (Polar Radius).......[6356773.32]
  Flattening...........................[0.003352329944944847]
  Flattening Inverse...................[298.2999932652668]
  First Eccentricity...................[0.08181333493893231]
  First Eccentricity Squared...........[0.006693421773829873]
  Second Eccentricity..................[0.08208852275188877]
  Second Eccentricity Squared..........[0.006738525567587472]"""


_GRS_80 = """
Ellipsoid [GRS 80]
  Semi-Major Axis (Equatorial Radius)..[6378137.0]
  Semi-Minor Axis (Polar Radius).......[6356752.3141]
  Flattening...........................[0.0033528106875095227]
  Flattening Inverse...................[298.2572215381486]
  First Eccentricity...................[0.08181919111988833]
  First Eccentricity Squared...........[0.006694380035512838]
  Second Eccentricity..................[0.08209443822977103]
  Second Eccentricity Squared..........[0.0067394967882615795]"""


_HELMERT_1906 = """
Ellipsoid [HELMERT 1906]
  Semi-Major Axis (Equatorial Radius)..[6378200.0]
  Semi-Minor Axis (Polar Radius).......[6356818.17]
  Flattening...........................[0.0033523298109184524]
  Flattening Inverse...................[298.3000051913226]
  First Eccentricity...................[0.08181333330622664]
  First Eccentricity Squared...........[0.006693421506675721]
  Second Eccentricity..................[0.08208852110265381]
  Second Eccentricity Squared..........[0.006738525296820739]"""


_HOUGH_1906 = """
Ellipsoid [HOUGH 1906]
  Semi-Major Axis (Equatorial Radius)..[6378270.0]
  Semi-Minor Axis (Polar Radius).......[6356794.343]
  Flattening...........................[0.0033670034351006867]
  Flattening Inverse...................[296.99999399320365]
  First Eccentricity...................[0.0819918908067705]
  First Eccentricity Squared...........[0.0067226701580694265]
  Second Eccentricity..................[0.08226889044349586]
  Second Eccentricity Squared..........[0.006768170334803944]"""


_INDONESIAN_1974 = """
Ellipsoid [INDONESIAN 1974]
  Semi-Major Axis (Equatorial Radius)..[6378160.0]
  Semi-Minor Axis (Polar Radius).......[6356774.504]
  Flattening...........................[0.0033529256086395256]
  Flattening Inverse...................[298.2469988070381]
  First Eccentricity...................[0.08182059097282252]
  First Eccentricity Squared...........[0.0066946091071419115]
  Second Eccentricity..................[0.08209585225822184]
  Second Eccentricity Squared..........[0.006739728958003832]"""


_INTERNATIONAL_1924 = """
Ellipsoid [INTERNATIONAL 1924]
  Semi-Major Axis (Equatorial Radius)..[6378388.0]
  Semi-Minor Axis (Polar Radius).......[6356911.946]
  Flattening...........................[0.003367003387062615]
  Flattening Inverse...................[296.9999982305938]
  First Eccentricity...................[0.0819918902228546]
  First Eccentricity Squared...........[0.006722670062316669]
  Second Eccentricity..................[0.08226888985364195]
  Second Eccentricity Squared..........[0.006768170237750658]"""


_KRASSOVSKY_1940 = """
Ellipsoid [KRASSOVSKY 1940]
  Semi-Major Axis (Equatorial Radius)..[6378245.0]
  Semi-Minor Axis (Polar Radius).......[6356863.019]
  Flattening...........................[0.0033523298336767685]
  Flattening Inverse...................[298.30000316622187]
  First Eccentricity...................[0.0818133335834678]
  First Eccentricity Squared...........[0.0066934215520398155]
  Second Eccentricity..................[0.08208852138270177]
  Second Eccentricity Squared..........[0.0067385253427982685]"""

_MODIFIED_AIRY = """
Ellipsoid [MODIFIED AIRY]
  Semi-Major Axis (Equatorial Radius)..[6377340.189]
  Semi-Minor Axis (Polar Radius).......[6356034.448]
  Flattening...........................[0.0033408506318589907]
  Flattening Inverse...................[299.32496546352854]
  First Eccentricity...................[0.08167337375652854]
  First Eccentricity Squared...........[0.006670539980773649]
  Second Eccentricity..................[0.0819471469341423]
  Second Eccentricity Squared..........[0.006715334890645988]
"""

_SOUTH_AMERICAN_1969 = """
Ellipsoid [SOUTH AMERICAN 1969]
  Semi-Major Axis (Equatorial Radius)..[6378160.0]
  Semi-Minor Axis (Polar Radius).......[6356774.719]
  Flattening...........................[0.003352891899858333]
  Flattening Inverse...................[298.249997276158]
  First Eccentricity...................[0.08182018036905428]
  First Eccentricity Squared...........[0.006694541915624534]
  Second Eccentricity..................[0.08209543749645103]
  Second Eccentricity Squared..........[0.006739660857733726]"""


_WGS_72 = """
Ellipsoid [WGS 72]
  Semi-Major Axis (Equatorial Radius)..[6378135.0]
  Semi-Minor Axis (Polar Radius).......[6356750.52]
  Flattening...........................[0.00335277945669078]
  Flattening Inverse...................[298.2599997755319]
  First Eccentricity...................[0.08181881069348491]
  First Eccentricity Squared...........[0.00669431778329633]
  Second Eccentricity..................[0.08209405395108844]
  Second Eccentricity Squared..........[0.006739433694124252]"""


_WGS_84 = """
Ellipsoid [WGS 84]
  Semi-Major Axis (Equatorial Radius)..[6378137.0]
  Semi-Minor Axis (Polar Radius).......[6356752.3142]
  Flattening...........................[0.0033528106718309896]
  Flattening Inverse...................[298.2572229328697]
  First Eccentricity...................[0.08181919092890624]
  First Eccentricity Squared...........[0.006694380004260827]
  Second Eccentricity..................[0.08209443803685366]
  Second Eccentricity Squared..........[0.006739496756586903]"""


## The Ellipsoids of Yor:
ELLIPSOIDS = (_AIRY_1830, _AUSTRALIAN_1965, _BESSEL_1841,
              _BESSEL_1841_NAMIBIA, _CLARKE_1866, _CLARKE_1880,
              _EVEREST_1830_INDIA, _EVEREST_1830_MALAYSIA,
              _EVEREST_1956_INDIA,
              _EVEREST_1964_MALAYSIA_SINGAPORE, _EVEREST_1969_MALAYSIA,
              _EVEREST_PAKISTAN, _FISHER_1960, _GRS_80, _HELMERT_1906,
              _HOUGH_1906, _INDONESIAN_1974, _INTERNATIONAL_1924,
              _INTERNATIONAL_1924, _KRASSOVSKY_1940, _MODIFIED_AIRY,
              _SOUTH_AMERICAN_1969, _WGS_72, _WGS_84)


## Convert string def to a class instance
def _string2ellipse(s, shownames=False):
    """Take a string and parse it into an Ellipsoid"""
    from .ellipsoid import Ellipsoid
    def get_lines(field):
        """split the field into items"""
        for item in s.split("\n"):
            if field in item:  # filter
                return item
            continue
        return ""

    def get_meat(line):
        """get the meat from between the buns"""
        slc = slice(1+line.index("["), line.index("]"))
        return line[slc]

    def get_value(line):
        """evaluate the meat"""
        return eval(get_meat(line))

    model = get_meat(get_lines("Ellipsoid")).replace(" ", "")
    if shownames:  # procedure
        return model
    a = get_value(get_lines("Equatorial"))
    finv = get_value(get_lines("Flattening Inverse."))

    name = s.split('\n')[1].lstrip("Ellipsoid ").lstrip('[').rstrip(']')

    return name, Ellipsoid(a, finv, model=model)


## The ALMANAC of ellipsoids
ALMANAC = dict(map(_string2ellipse, ELLIPSOIDS))


WGS84 = ALMANAC['WGS 84']
## Add Sphereical Earth
ALMANAC["SPHERE"] = WGS84.from_axes(a=pow(WGS84.a*WGS84.b, 0.5),
                                    model="SPHERE")

