import requests
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from collections import OrderedDict
import numpy as np
candidates = Table.read("all-mgk-candidates.csv")

# Only propose to observe V < 14
keep = (~candidates["Vmag"].mask) \
     * (candidates["Vmag"] < 14)

candidates = candidates[keep]
candidates.sort("ra")

def exposure_time(vmag):

    # see exposure-calculations.numbers
    exps = dict([
        (9, 0),
        (10, 200),
        (11, 500),
        (12, 1200),
        (13, 2400),
        (14, 3600)
    ])

    for mag, exp in exps.items():
        if mag > np.round(vmag, 1):
            print(vmag, mag, exp)
            return exp

    return exp

#\target{}{}{}{}
remove_sep = lambda s: s.replace(":", "").replace(".", "")

output = []
total_exp = 0

for candidate in candidates:
    
    position = SkyCoord(candidate["ra"], candidate["dec"], unit=(u.deg, u.deg))

    ra, dec = position.to_string("hmsdms", sep=":").split()
    identifier = "J{}{}".format(remove_sep(ra[:11]), remove_sep(dec[:11]))

    vmag = candidate["Vmag"]
    t_exp = exposure_time(vmag)

    try:
        comment = "V={vmag:.1f} $t_{{\\rm exp}} = {exp:.0f}$ s".format(
            vmag=vmag, exp=t_exp)

    except TypeError:
        comment = "V=\\todo{{XX.X}} $t_{{\\rm exp}} = XX$ s"

    line = "\\target{{{identifier}}}{{{ra}}}{{{dec}}}{{{comment}}}".format(
        identifier=identifier, ra=ra[:11], dec=dec[:11], comment=comment)

    output.append(line)

    total_exp += t_exp

with open("candidates.tex", "w") as fp:
    fp.write("\n".join(output))


candidates.write("candidates.csv")

raise a


# Check SMOKA for each candidate.
smoka_url = "http://smoka.nao.ac.jp/fssearch"

N = len(candidates)

for i, candidate in enumerate(candidates):

    print(i, N)

    #ra = "21h03m52.11s"
    #dec = "-29d42m50.2s"
    ra = "{:.5f}".format(candidate["ra"])
    dec = "{:.5f}".format(candidate["dec"])

    data  = [
        ("object", ""),
        ("resolver", "SIMBAD"),
        ("coordsys", "Equatorial"),
        ("equinox", "J2000"),
        ("fieldofview", "auto"),
        ("RadOrRec", "radius"),
        ("longitudeC", ra),
        ("latitudeC", dec),
        ("radius", "0.08"),
        ("longitudeF", ""),
        ("latitudeF", ""),
        ("longitudeT", ""),
        ("latitudeT", ""),
        ("date_obs", ""),
        ("exptime", ""),
        ("observer", ""),
        ("prop_id", ""),
        ("frameid", ""),
        ("exp_id", ""),
        ("dataset", ""),
        ("asciitable", "Ascii"),
        ("frameorshot", "Frame"),
        ("action", "Search"),
        ("instruments", "HDS"),
        ("multiselect_0", "HDS"),
        ("obs_mod", "IMAG"),
        ("obs_mod", "SPEC"),
        ("obs_mod", "IPOL"),
        ("multiselect_1", "IMAG"),
        ("multiselect_1", "SPEC"),
        ("multiselect_1", "IPOL"),
        ("data_typ", "OBJECT"),
        ("multiselect_2", "OBJECT"),
        ("obs_cat", "Science+Observation"),
        ("multiselect_3", "Science+Observation"),
        ("bandwidth_type", "FILTER"),
        ("band", ""),
        ("dispcol", "FRAMEID"),
        ("dispcol", "DATE_OBS"),
        ("dispcol", "FITS_SIZE"),
        ("dispcol", "OBS_MODE"),
        ("dispcol", "DATA_TYPE"),
        ("dispcol", "OBJECT"),
        ("dispcol", "WVLEN"),
        ("dispcol", "DISPERSER"),
        ("dispcol", "ECLLONG"),
        ("dispcol", "ECLLAT"),
        ("dispcol", "UT_START"),
        ("dispcol", "EXPTIME"),
        ("orderby", "FRAMEID"),
        ("diff", "100"),
        ("output_equinox", "J2000"),
        ("from", "0"),
    ]

    check = requests.post(smoka_url, data=data)

    assert check.ok

    if "No matching frame were found" in check.text:
        print("No data")

    else:
        # You gotta justify this in your proposal. How much data exists?
        raise a
