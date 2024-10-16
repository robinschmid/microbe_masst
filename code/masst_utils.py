import requests
import logging
from datetime import timedelta
import requests_cache
from enum import Enum, auto
import json
import pandas as pd
from dataclasses import dataclass
import usi_utils

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# requests_cache.install_cache("fastmasst_cache", expire_after=timedelta(days=2))


@dataclass
class SpecialMasst:
    prefix: str
    root: str
    tree_file: str
    metadata_file: str
    tree_node_key: str
    metadata_key: str


MICROBE_MASST = SpecialMasst(
    prefix="microbe",
    root="microbes",
    tree_file="../data/microbe_masst_tree.json",
    metadata_file="../data/microbe_masst_table.csv",
    tree_node_key="NCBI",
    metadata_key="Taxa_NCBI",
)
FOOD_MASST = SpecialMasst(
    prefix="food",
    root="food",
    tree_file="../data/food_masst_tree.json",
    metadata_file="../data/food_masst_table.csv",
    tree_node_key="name",
    metadata_key="node_id",
)
PLANT_MASST = SpecialMasst(
    prefix="plant",
    root="plants",
    tree_file="../data/plant_masst_tree.json",
    metadata_file="../data/plant_masst_table.csv",
    tree_node_key="NCBI",
    metadata_key="Taxa_NCBI",
)
GLOBAL_MASST = SpecialMasst(
    prefix="global",
    root="global",
    tree_file="../data/global_masst_tree.json",
    metadata_file="../data/global_masst_table.csv",
    tree_node_key="ID",
    metadata_key="ID",
)
PERSONALCAREPRODUCT_MASST = SpecialMasst(
    prefix="personalCareProduct",
    root="personal care product",
    tree_file="../data/personalCareProduct_masst_tree.json",
    metadata_file="../data/personalCareProduct_masst_table.csv",
    tree_node_key="name",
    metadata_key="node_id",
)
MICROBIOME_MASST = SpecialMasst(
    prefix="microbiome",
    root="microbiome",
    tree_file="../data/microbiome_masst_tree.json",
    metadata_file="../data/microbiome_masst_table.tsv",
    tree_node_key="name",
    metadata_key="node_id",
)

URL = "https://fasst.gnps2.org/search"
SPECIAL_MASSTS = [FOOD_MASST, MICROBE_MASST, PLANT_MASST, GLOBAL_MASST, PERSONALCAREPRODUCT_MASST, MICROBIOME_MASST]


class DataBase(Enum):
    metabolomicspanrepo_index_latest = auto()  # all gnps data
    gnpsdata_index = auto()  # all gnps data
    gnpsdata_index_11_25_23 = auto()  # all gnps data
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
    analog_mass_below=130,
    analog_mass_above=200,
    database=DataBase.metabolomicspanrepo_index_latest,
):
    if str(usi_or_lib_id).startswith("CCMS"):
        # handle library ID
        usi_or_lib_id = "mzspec:GNPS:GNPS-LIBRARY:accession:{}".format(usi_or_lib_id)

    try:
        # trying to get database name, check if string or enum
        if isinstance(database, DataBase):
            database = database.name

        params = {
            "usi": usi_or_lib_id,
            "library": str(database),
            "analog": "Yes" if analog else "No",
            "delta_mass_below": analog_mass_below,
            "delta_mass_above": analog_mass_above,
            "pm_tolerance": precursor_mz_tol,
            "fragment_tolerance": mz_tol,
            "cosine_threshold": min_cos,
        }

        return _fast_masst(params)
    # except requests.exceptions.Timeout:
    except Exception as e:
        # logging.exception("Failed fastMASST {}".format(usi_or_lib_id))
        raise e


def fast_masst_spectrum(
    mzs,
    intensities,
    precursor_mz,
    precursor_charge=1,
    precursor_mz_tol=0.05,
    mz_tol=0.05,
    min_cos=0.7,
    analog=False,
    analog_mass_below=130,
    analog_mass_above=200,
    database=DataBase.metabolomicspanrepo_index_latest,
    min_signals=3,
):
    """

    :param mzs:
    :param intensities:
    :param precursor_mz:
    :param precursor_charge:
    :param precursor_mz_tol:
    :param mz_tol:
    :param min_cos:
    :param analog:
    :param analog_mass_below:
    :param analog_mass_above:
    :param database:
    :return: (MASST results as json, filtered data points as array of array [[x,y],[...]]
    """
    # relative intensity and precision
    dps = [[round(mz, 5), intensity] for mz, intensity in zip(mzs, intensities)]
    spec_dict = {
        "n_peaks": len(dps),
        "peaks": dps,
        "precursor_mz": precursor_mz,
        "precursor_charge": abs(precursor_charge),
    }
    return fast_masst_spectrum_dict(
        spec_dict,
        precursor_mz_tol,
        mz_tol,
        min_cos,
        analog,
        analog_mass_below,
        analog_mass_above,
        database,
        min_signals,
    )


def fast_masst_spectrum_dict(
    spec_dict: dict,
    precursor_mz_tol=0.05,
    mz_tol=0.05,
    min_cos=0.7,
    analog=False,
    analog_mass_below=130,
    analog_mass_above=200,
    database=DataBase.metabolomicspanrepo_index_latest,
    min_signals=3,
):
    """

    :param spec_dict: dictionary with GNPS json like spectrum. Example:
    {"n_peaks":8,"peaks":[[80.97339630126953,955969.8125],[81.98159790039062,542119.1875],[98.98410034179688,
    483893632.0],[116.99469757080078,1605324.25],[127.0155029296875,182958080.0],[131.01019287109375,878951.375],
    [155.0467071533203,73527152.0],[183.07809448242188,16294011.0]],"precursor_charge":1,"precursor_mz":183.078,
    "splash":"splash10-0002-9500000000-fe91737b5df956c7e69e"}
    :param precursor_mz_tol:
    :param mz_tol:
    :param min_cos:
    :param analog:
    :param analog_mass_below:
    :param analog_mass_above:
    :param database:
    :return: (MASST results as json, filtered data points as array of array [[x,y],[...]]
    """
    try:
        max_intensity = max([v[1] for v in spec_dict["peaks"]])
        dps = [
            [round(dp[0], 5), round(dp[1] / max_intensity * 100.0, 1)]
            for dp in spec_dict["peaks"]
        ]
        dps = [dp for dp in dps if dp[1] >= 0.1]

        spec_dict["peaks"] = dps
        spec_dict["n_peaks"] = len(dps)
        if spec_dict["n_peaks"] < min_signals:
            return None, dps

        spec_json = json.dumps(spec_dict)

        # trying to get database name, check if string or enum
        if isinstance(database, DataBase):
            database = database.name

        params = {
            "library": str(database),
            "analog": "Yes" if analog else "No",
            "delta_mass_below": analog_mass_below,
            "delta_mass_above": analog_mass_above,
            "pm_tolerance": precursor_mz_tol,
            "fragment_tolerance": mz_tol,
            "cosine_threshold": min_cos,
            "query_spectrum": spec_json,
        }
        return _fast_masst(params), dps
    except Exception as e:
        # logging.exception("Failed fastMASST on spectrum.")
        raise e


def _fast_masst(params):
    """

    :param params: dict of the query input and parameters
    :return: dict with the masst results. [results] contains the individual matches, [grouped_by_dataset] contains
    all datasets and their titles
    """
    search_api_response = requests.post(URL, data=params, timeout=300)
    logging.debug("fastMASST response={}".format(search_api_response.status_code))
    search_api_response.raise_for_status()
    search_api_response_json = search_api_response.json()
    return search_api_response_json


def filter_matches(df, precursor_mz_tol, min_matched_signals, analog):
    # DO NOT FILTER BY MZ FOR ANALOG
    if analog:
        return df.loc[df["Matching Peaks"] >= min_matched_signals]
    else:
        return df.loc[
            (
                df["Delta Mass"].between(
                    -precursor_mz_tol, precursor_mz_tol, inclusive="both"
                )
            )
            & (df["Matching Peaks"] >= min_matched_signals)
        ]


def extract_matches_from_masst_results(
    results_dict,
    precursor_mz_tol,
    min_matched_signals,
    analog,
    limit_to_best_match_in_file: bool = False,
    add_dataset_titles=False,
) -> pd.DataFrame:
    """
    :param results_dict: masst results
    :param add_dataset_titles: add dataset titles to each row
    :return: DataFrame of the individual matches
    """
    masst_df = pd.DataFrame(results_dict["results"])
    try:
        masst_df.drop(
            columns=[
                "Unit Delta Mass",
                "Query Scan",
                "Query Filename",
                "Index UnitPM",
                "Index IdxInUnitPM",
                "Filtered Input Spectrum Path",
            ],
            inplace=True,
            axis=1,
        )
    except Exception as e:
        # fastMASST response is sometimes empty
        return masst_df

    masst_df = filter_matches(masst_df, precursor_mz_tol, min_matched_signals, analog)
    if add_dataset_titles:
        datasets = results_dict["grouped_by_dataset"]
        dataset_info_dict = dict([(e["Dataset"], e["title"]) for e in datasets])
        # might not be in the df
        masst_df.drop(columns=["mzs", "intensities"], inplace=True, axis=1)

        for match in masst_df:
            match["dataset_title"] = dataset_info_dict.get(match["Dataset"], None)

    if limit_to_best_match_in_file:
        # create a usi column that only points to the dataset:file (not scan)
        masst_df["file_usi"] = [
            usi_utils.ensure_simple_file_usi(usi) for usi in masst_df["USI"]
        ]
        masst_df = masst_df.sort_values(
            by=["Cosine", "Matching Peaks"], ascending=[False, False]
        ).drop_duplicates("file_usi")
    return masst_df


def extract_datasets_from_masst_results(
    results_dict, matches_df: pd.DataFrame
) -> pd.DataFrame:
    datasets_df = pd.DataFrame(results_dict["grouped_by_dataset"])
    # recalc frequency with filtered MASST results
    new_dataset_df = matches_df.groupby("Dataset").size().reset_index(name="Frequency")
    # transfer dataset title
    new_dataset_df.merge(datasets_df, on="Dataset", how="left")
    return new_dataset_df


# example
# https://fastlibrarysearch.ucsd.edu/fastsearch/?usi1=mzspec%3AGNPS%3AGNPS-LIBRARY%3Aaccession%3ACCMSLIB00000579622
# &precursor_mz=183.078&charge=1&library_select=metabolomicspanrepo_index_latest&analog_select=No&delta_mass_below=130&delta_mass_above
# =200&pm_tolerance=0.05&fragment_tolerance=0.05&cosine_threshold=0.7&use_peaks=3#%7B%22peaks%22%3A%20%2280.9734
# %5Ct955969.8%5Cn81.9816%5Ct542119.2%5Cn98.9841%5Ct483893630.0%5Cn116.9947%5Ct1605324.2%5Cn127.0155%5Ct182958080.0
# %5Cn131.0102%5Ct878951.4%5Cn155.0467%5Ct73527150.0%5Cn183.0781%5Ct16294011.0%22%7D
if __name__ == "__main__":
    usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005883948"
    matches = fast_masst(usi)

    print(len(matches))
