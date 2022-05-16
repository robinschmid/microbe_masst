import requests
import logging
from datetime import timedelta
import requests_cache
from enum import Enum, auto

requests_cache.install_cache('fastmasst_cache', expire_after=timedelta(days=2))

class DataBase(Enum):
     gnpsdata_index = auto()
     gnpslibrary = auto()
     massivedata_index = auto()
     massivekb_index = auto()

# based on
# https://github.com/mwang87/GNPS_LCMSDashboard/blob/a9971fa557c735c8e0ccd7681653eebd415a8636/app.py#L1632
# usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00000001556"


# http://fastlibrarysearch.ucsd.edu/fastsearch/?usi1=mzspec%3AGNPS%3AGNPS-LIBRARY%3Aaccession%3ACCMSLIB00000001556
# &precursor_mz=None&charge=None&library_select=gnpsdata_index&analog_select=No&delta_mass_below=130&delta_mass_above=200&pm_tolerance=0.05&fragment_tolerance=0.05&cosine_threshold=0.7&use_peaks=0#%7B%22peaks%22%3A%20null%7D
def fast_masst(usi_or_lib_id, precursor_mz_tol=0.05, mz_tol=0.02, min_cos=0.7):
    if str(usi_or_lib_id).startswith("CCMS"):
        # handle library ID
        usi_or_lib_id = "mzspec:GNPS:GNPS-LIBRARY:accession:{}".format(usi_or_lib_id)

    try:
        search_api_url = "https://fastlibrarysearch.ucsd.edu/search?usi={}&library=gnpsdata_index&analog=No&pm_tolerance={}&fragment_tolerance={}&cosine_threshold={}"\
            .format(usi_or_lib_id, precursor_mz_tol, mz_tol, min_cos)
        # search_api_url = "https://fastlibrarysearch.ucsd.edu/search?usi={}&library=massivekb_index&analog=No".format(
        #     usi)
        search_api_response = requests.get(search_api_url, timeout=50)
        search_api_response_json = search_api_response.json()
        matches = search_api_response_json["results"]
        return matches
    # except requests.exceptions.Timeout:
    except Exception:
        logging.exception("Failed fastMASST.")


if __name__ == '__main__':
    usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005883671"
    matches = fast_masst(usi)
    print("done")