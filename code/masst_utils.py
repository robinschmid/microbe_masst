import requests
import logging
from datetime import timedelta
import requests_cache
from enum import Enum, auto
import json

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# requests_cache.install_cache("fastmasst_cache", expire_after=timedelta(days=2))

URL = "https://fastlibrarysearch.ucsd.edu/search"

# define the max number of datapoints as an estimate
# limit for URL is 2048
DATAPOINT_CHARACTERS = (10 + 5 + 10)
MAX_N_DATAPOINTS = int((2048 - len("https://fastlibrarysearch.ucsd.edu/search?library=gnpsdata_index&analog=Yes"
                              "&delta_mass_below=130" \
           "&delta_mass_above=200&pm_tolerance=0.05&fragment_tolerance=0.02&cosine_threshold=0.7&query_spectrum=%7B" \
           "%22n_peaks%22%3A+0%2C+%22peaks%22%3A+%5B%5D%2C+%22precursor_mz%22%3A+325.2349%2C+%22precursor_charge%22" \
           "%3A+1%7D")) / DATAPOINT_CHARACTERS)

class DataBase(Enum):
    gnpsdata_index = auto()  # all gnps data
    gnpslibrary = auto()  # gnps library
    massivedata_index = auto()
    massivekb_index = auto()

# based on
# https://github.com/mwang87/GNPS_LCMSDashboard/blob/a9971fa557c735c8e0ccd7681653eebd415a8636/app.py#L1632
# usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00000001556"
def fast_masst(
    usi_or_lib_id,
    precursor_mz_tol=0.05,
    mz_tol=0.02,
    min_cos=0.7,
    analog=False,
    database=DataBase.gnpsdata_index,
):
    if str(usi_or_lib_id).startswith("CCMS"):
        # handle library ID
        usi_or_lib_id = "mzspec:GNPS:GNPS-LIBRARY:accession:{}".format(usi_or_lib_id)

    try:
        params = {
            "usi": usi_or_lib_id,
            "library": str(database.name),
            "analog": "Yes" if analog else "No",
            "pm_tolerance": precursor_mz_tol,
            "fragment_tolerance": mz_tol,
            "cosine_threshold": min_cos,
        }

        return _fast_masst(params)
    # except requests.exceptions.Timeout:
    except Exception:
        logging.exception("Failed fastMASST {}".format(usi_or_lib_id))


def _fast_masst(params):
    search_api_response = requests.get(URL, params=params, timeout=150)
    search_api_response.raise_for_status()
    search_api_response_json = search_api_response.json()
    datasets = search_api_response_json["grouped_by_dataset"]
    dataset_info_dict = dict([(e["Dataset"], e["title"]) for e in datasets])
    matches = search_api_response_json["results"]
    for match in matches:
        match["dataset_title"] = dataset_info_dict.get(match["Dataset"], None)
    return matches



def fast_masst_spectrum(
    mzs,
    intensities,
    precursor_mz,
    precursor_charge=1,
    precursor_mz_tol=0.05,
    mz_tol=0.05,
    min_cos=0.7,
    analog=True,
    database=DataBase.gnpsdata_index,
):
    try:
        max_intensity = max(intensities)
        dps = [[round(mz, 4), round(intensity/max_intensity*100.0, 1)] for mz, intensity in zip(mzs, intensities)]

        total_dps = len(dps)
        if total_dps > MAX_N_DATAPOINTS:
            sorted_dps = sorted(dps, key=lambda dp: -dp[1])[0:MAX_N_DATAPOINTS]
            dps = sorted(sorted_dps, key=lambda dp: dp[0])
            logger.info("Number of signals in spectrum was reduced to meet length limit {}-->{}".format(total_dps,
                                                                                                        len(dps)))


        spec_json = json.dumps(
            {
                "n_peaks": len(dps),
                "peaks": dps,
                "precursor_mz": precursor_mz,
                "precursor_charge": abs(precursor_charge),
            }
        )

        params = {
            "library": database.name,
            "analog": "Yes" if analog else "No",
            "delta_mass_below": 130,
            "delta_mass_above": 200,
            "pm_tolerance": precursor_mz_tol,
            "fragment_tolerance": mz_tol,
            "cosine_threshold": min_cos,
            "query_spectrum": spec_json,
        }
        return _fast_masst(params)
    except Exception:
        logging.exception("Failed fastMASST on spectrum.")


# example
# https://fastlibrarysearch.ucsd.edu/fastsearch/?usi1=mzspec%3AGNPS%3AGNPS-LIBRARY%3Aaccession%3ACCMSLIB00000579622&precursor_mz=183.078&charge=1&library_select=gnpsdata_index&analog_select=No&delta_mass_below=130&delta_mass_above=200&pm_tolerance=0.05&fragment_tolerance=0.05&cosine_threshold=0.7&use_peaks=3#%7B%22peaks%22%3A%20%2280.9734%5Ct955969.8%5Cn81.9816%5Ct542119.2%5Cn98.9841%5Ct483893630.0%5Cn116.9947%5Ct1605324.2%5Cn127.0155%5Ct182958080.0%5Cn131.0102%5Ct878951.4%5Cn155.0467%5Ct73527150.0%5Cn183.0781%5Ct16294011.0%22%7D
if __name__ == "__main__":
    usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005883671"
    matches = fast_masst(usi)

    print(len(matches))
