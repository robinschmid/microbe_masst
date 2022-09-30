import sys
import logging
import pandas as pd
from tqdm import tqdm
import re
import argparse

from masst_tree import create_enriched_masst_tree
import masst_utils as masst
import usi_utils

MATCH_COLUMNS = ["USI", "Cosine", "Matching Peaks", "Status"]

LIB_COLUMNS = ["USI", "GNPSLibraryAccession", "Cosine", "Matching Peaks", "CompoundName", "Adduct", "Charge"]

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# activate pandas tqdm progress_apply
tqdm.pandas()


def process_matches(file_name, compound_name, matches, library_matches, precursor_mz_tol, min_matched_signals,
                    input_label, params_label, usi=None):
    common_file = "../{}_{}".format(file_name, compound_name.replace(" ", "_"))

    # extract results
    matches_df = masst.extract_matches_from_masst_results(matches,
                                                          precursor_mz_tol,
                                                          min_matched_signals,
                                                          False
                                                          )

    # create a usi column that only points to the dataset:file (not scan)
    matches_df["file_usi"] = [usi_utils.ensure_simple_file_usi(usi) for usi in matches_df["USI"]]
    matches_df = matches_df.sort_values(by=["Cosine", "Matching Peaks"], ascending=[False, False]).drop_duplicates(
        "file_usi")

    matches_df[MATCH_COLUMNS].to_csv("{}_{}.tsv".format(common_file, "matches"), index=False, sep="\t")

    lib_matches_df = masst.extract_matches_from_masst_results(library_matches,
                                                              precursor_mz_tol,
                                                              min_matched_signals,
                                                              False
                                                              )
    lib_matches_df[LIB_COLUMNS].to_csv("{}_{}.tsv".format(common_file, "library"), index=False, sep="\t")

    if "grouped_by_dataset" not in matches:
        logger.debug("Missing datasets")
    # extract matches
    datasets_df = masst.extract_datasets_from_masst_results(matches, matches_df)
    datasets_df.to_csv("{}_{}.tsv".format(common_file, "datasets"), index=False, sep="\t")

    # add library matches to table
    lib_match_json = lib_matches_df.to_json(orient="records")

    # microbeMASST
    logger.debug("Exporting microbeMASST %s", compound_name)
    create_enriched_masst_tree(matches_df,
                               masst.MICROBE_MASST,
                               common_file=common_file,
                               lib_match_json=lib_match_json,
                               input_str=input_label,
                               parameter_str=params_label,
                               usi=usi,
                               format_out_json=True,
                               compress_out_html=False
                               )

    # foodMASST
    logger.debug("Exporting foodMASST %s", compound_name)
    create_enriched_masst_tree(matches_df,
                               masst.FOOD_MASST,
                               common_file=common_file,
                               lib_match_json=lib_match_json,
                               input_str=input_label,
                               parameter_str=params_label,
                               usi=usi,
                               format_out_json=True,
                               compress_out_html=False
                               )

    return matches_df


def query_usi_or_id(file_name, usi_or_lib_id, compound_name,
                    precursor_mz_tol=0.05,
                    mz_tol=0.02,
                    min_cos=0.7,
                    min_matched_signals=3,
                    analog=False,
                    analog_mass_below=150,
                    analog_mass_above=200):
    # might raise exception for service
    try:
        logger.debug("Query fastMASST id:%s  of %s", usi_or_lib_id, compound_name)
        matches = masst.fast_masst(usi_or_lib_id, precursor_mz_tol=precursor_mz_tol, mz_tol=mz_tol, min_cos=min_cos,
                                   analog=analog, analog_mass_below=analog_mass_below,
                                   analog_mass_above=analog_mass_above, database=masst.DataBase.gnpsdata_index)

        if not matches or "results" not in matches or len(matches["results"]) == 0:
            logger.debug("Empty fastMASST response for compound %s with id %s", compound_name, usi_or_lib_id)
            return None

        library_matches = masst.fast_masst(usi_or_lib_id, precursor_mz_tol=precursor_mz_tol, mz_tol=mz_tol,
                                           min_cos=min_cos,
                                           analog=False, database=masst.DataBase.gnpslibrary)

        params_label = create_params_label(analog, analog_mass_above, analog_mass_below, min_cos, min_matched_signals,
                                           mz_tol, precursor_mz_tol)
        input_label = "ID: {};  Descriptor: {}".format(usi_or_lib_id, compound_name)
        return process_matches(file_name, compound_name, matches, library_matches, precursor_mz_tol,
                               min_matched_signals, input_label, params_label, usi_utils.ensure_usi(usi_or_lib_id))
    except Exception as e:
        logger.exception(e)
        return None


def query_spectrum(file_name, compound_name, precursor_mz, precursor_charge, mzs, intensities,
                   precursor_mz_tol=0.05,
                   mz_tol=0.02,
                   min_cos=0.7,
                   min_matched_signals=3,
                   analog=False,
                   analog_mass_below=150,
                   analog_mass_above=200):
    # might raise exception for service
    try:
        matches, filtered_dps = masst.fast_masst_spectrum(
            mzs=mzs,
            intensities=intensities,
            precursor_mz=precursor_mz,
            precursor_charge=precursor_charge,
            precursor_mz_tol=precursor_mz_tol,
            mz_tol=mz_tol,
            min_cos=min_cos,
            analog=analog,
            analog_mass_below=analog_mass_below,
            analog_mass_above=analog_mass_above,
            database=masst.DataBase.gnpsdata_index
        )
        if not matches or "results" not in matches or len(matches["results"]) == 0:
            logger.debug("Empty fastMASST response for spectrum %s", compound_name)
            return None

        library_matches, _ = masst.fast_masst_spectrum(
            mzs=mzs,
            intensities=intensities,
            precursor_mz=precursor_mz,
            precursor_charge=precursor_charge,
            precursor_mz_tol=precursor_mz_tol,
            mz_tol=mz_tol,
            min_cos=min_cos,
            analog=False,
            database=masst.DataBase.gnpslibrary
        )

        params_label = create_params_label(analog, analog_mass_above, analog_mass_below, min_cos, min_matched_signals,
                                           mz_tol, precursor_mz_tol)
        input_label = "Descriptor: {};  Precursor m/z: {};  Data points:{}".format(compound_name, round(precursor_mz,
                                                                                                        5),
                                                                                   len(filtered_dps))
        return process_matches(file_name, compound_name, matches, library_matches, precursor_mz_tol,
                               min_matched_signals, input_label, params_label)
    except Exception as e:
        return None


def path_safe(file):
    return re.sub('[^-a-zA-Z0-9_.() ]+', '_', file)


def save_masst_results(jobs_df, out_filename_no_ext):
    finished_jobs_tsv = "../output/finished.tsv"
    masst_results_tsv = "../{}.tsv".format(out_filename_no_ext)
    # explode rows and export
    # export_masst_results(jobs_df, masst_results_tsv)
    # add tree links
    jobs_df["Result"] = jobs_df.progress_apply(
        lambda row: "Success" if row["fastMASST"] is not None else None,
        axis=1)
    # export the list of links
    logger.info("Exporting finished_jobs_tsv %s", finished_jobs_tsv)
    jobs_df.drop(columns=["fastMASST"], axis=1, inplace=True)
    jobs_df.to_csv(finished_jobs_tsv, sep="\t", index=False)


def export_masst_results(jobs_df: pd.DataFrame, masst_results_tsv: str):
    # export masst results
    logger.info("Exporting masst results summary to %s", masst_results_tsv)
    # jobs_df = explode_masst_columns(jobs_df)
    # drop data points
    try:
        jobs_df.drop(columns=["Query Scan", "Query Filename", "Index UnitPM", "Index IdxInUnitPM",
                              "Filtered Input Spectrum Path"], inplace=True, axis=1)
        # might not be in the df
        jobs_df.drop(columns=["mzs", "intensities"], inplace=True, axis=1)
    except Exception as e:
        logger.exception("Error removing columns from export table", e)
    jobs_df.to_csv(masst_results_tsv, sep="\t", index=False)


def create_params_label(analog, analog_mass_above, analog_mass_below, min_cos, min_matched_signals, mz_tol,
                        precursor_mz_tol):
    params_label = "min matched signals: {};  min cosine: {};  precursor m/z tolerance: {};  m/z " \
                   "tolerance: {};  analog: {}".format(min_matched_signals, min_cos, precursor_mz_tol, mz_tol,
                                                       analog)
    if analog:
        params_label += ";  analogs below m/z:{};  analogs above m/z:{}".format(analog_mass_below,
                                                                                analog_mass_above)
    return params_label


if __name__ == '__main__':
    # parsing the arguments (all optional)
    parser = argparse.ArgumentParser(description='Run fast microbeMASST in batch')
    parser.add_argument('--usi_or_lib_id', type=str,
                        help='usi or GNPS library ID',
                        default="CCMSLIB00006582001")
    parser.add_argument('--out_file', type=str, help='output html and other files, name without extension',
                        default="output/aaaafastMASST_")

    # only for USI or lib ID file
    parser.add_argument('--compound_name', type=str, help='compound name',
                        default="phe_ca")

    # MASST params
    parser.add_argument('--precursor_mz_tol', type=float,
                        help='precursor mz tolerance',
                        default="0.05")
    parser.add_argument('--mz_tol', type=float,
                        help='mz tolerance to match signals',
                        default="0.05")
    parser.add_argument('--min_cos', type=float,
                        help='Minimum cosine score for a match',
                        default="0.7")
    parser.add_argument('--min_matched_signals', type=int,
                        help='Minimum matched signals',
                        default="3")
    parser.add_argument('--analog', type=bool,
                        help='Search for analogs within mass window',
                        default=False)
    parser.add_argument('--analog_mass_below', type=float,
                        help='Maximum mass delta for analogs',
                        default="150")
    parser.add_argument('--analog_mass_above', type=float,
                        help='Maximum mass delta for analogs',
                        default="200")

    args = parser.parse_args()

    try:
        usi_or_lib_id = args.usi_or_lib_id
        compound_name = args.compound_name
        out_file = args.out_file

        # MASST parameters
        precursor_mz_tol = args.precursor_mz_tol
        mz_tol = args.mz_tol
        min_cos = args.min_cos
        min_matched_signals = args.min_matched_signals
        analog = args.analog
        analog_mass_below = args.analog_mass_below
        analog_mass_above = args.analog_mass_above

        query_usi_or_id(
            out_file,
            usi_or_lib_id,
            compound_name,
            precursor_mz_tol=precursor_mz_tol,
            mz_tol=mz_tol,
            min_cos=min_cos,
            min_matched_signals=min_matched_signals,
            analog=analog,
            analog_mass_below=analog_mass_below,
            analog_mass_above=analog_mass_above
        )
    except Exception as e:
        # exit with error
        logger.exception(e)
        sys.exit(1)

    # exit with OK
    sys.exit(0)
