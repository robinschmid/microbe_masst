import sys
import logging
import pandas as pd
from tqdm import tqdm
import re
import argparse
import pyteomics.mgf
import bundle_to_html

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait

import microbe_masst as microbemasst
import masst_utils as masst

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# activate pandas tqdm progress_apply
tqdm.pandas()


def process_matches(file_name, compound_name, matches, library_matches, precursor_mz_tol, min_matched_signals,
                    input_label, params_label):
    common_file = "../{}_{}".format(file_name, compound_name.replace(" ", "_"))

    # extract results
    matches_df = masst.extract_matches_from_masst_results(matches,
                                                          precursor_mz_tol,
                                                          min_matched_signals,
                                                          False
                                                          )
    matches_df.to_csv("{}_{}.tsv".format(common_file, "matches"), index=False, sep="\t")

    lib_matches_df = masst.extract_matches_from_masst_results(library_matches,
                                                              precursor_mz_tol,
                                                              min_matched_signals,
                                                              False
                                                              )
    lib_matches_df.to_csv("{}_{}.tsv".format(common_file, "library"), index=False, sep="\t")

    if "grouped_by_dataset" not in matches:
        logger.debug("Missing datasets")
    # extract matches
    datasets_df = masst.extract_datasets_from_masst_results(matches, matches_df)
    datasets_df.to_csv("{}_{}.tsv".format(common_file, "datasets"), index=False, sep="\t")

    # add library matches to table
    lib_match_json = lib_matches_df.to_json(orient="records")

    # microbeMASST
    out_html = "{}_{}.html".format(common_file, "microbes")
    out_json_tree = "{}_{}.json".format(common_file, "microbes")
    out_counts = "{}_{}.json".format(common_file, "counts_microbes")
    replace_dict = {
        "PLACEHOLDER_JSON_DATA": out_json_tree,
        "LIBRARY_JSON_DATA_PLACEHOLDER": lib_match_json,
        "INPUT_LABEL_PLACEHOLDER": input_label,
        "PARAMS_PLACEHOLDER": params_label
    }

    logger.debug("Exporting microbeMASST %s", compound_name)
    result = microbemasst.run_microbe_masst_for_matches(matches_df,
                                                        in_html="../code/collapsible_tree_v3.html",
                                                        in_ontology="../data/ncbi.json",
                                                        metadata_file="../data/microbe_masst_table.csv",
                                                        out_counts_file=out_counts,
                                                        out_json_tree=out_json_tree,
                                                        format_out_json=True,
                                                        out_html=out_html,
                                                        compress_out_html=True,
                                                        replace_dict=replace_dict,
                                                        node_key="NCBI",
                                                        data_key="ncbi"
                                                        )

    # foodMASST
    out_html = "{}_{}.html".format(common_file, "food")
    out_json_tree = "{}_{}.json".format(common_file, "food")
    out_counts = "{}_{}.json".format(common_file, "counts_food")
    replace_dict["PLACEHOLDER_JSON_DATA"] = out_json_tree

    logger.debug("Exporting foodMASST %s", compound_name)
    result = microbemasst.run_food_masst_for_matches(matches_df,
                                                     in_html="../code/collapsible_tree_v3.html",
                                                     out_counts_file=out_counts,
                                                     out_json_tree=out_json_tree,
                                                     format_out_json=True,
                                                     out_html=out_html,
                                                     compress_out_html=True,
                                                     replace_dict=replace_dict,
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
            return None

        library_matches = masst.fast_masst(usi_or_lib_id, precursor_mz_tol=precursor_mz_tol, mz_tol=mz_tol,
                                           min_cos=min_cos,
                                           analog=False, database=masst.DataBase.gnpslibrary)

        params_label = create_params_label(analog, analog_mass_above, analog_mass_below, min_cos, min_matched_signals,
                                           mz_tol, precursor_mz_tol)
        input_label = "ID: {};  Descriptor: {}".format(usi_or_lib_id, compound_name)
        return process_matches(file_name, compound_name, matches, library_matches, precursor_mz_tol,
                               min_matched_signals, input_label, params_label)
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


def explode_masst_columns(jobs_df):
    # multiple matches for each row split into multiple rows
    masst_df = jobs_df.explode("fastMASST", ignore_index=True)
    # each match is a dict so split into columns
    masst_df = pd.concat([masst_df.drop(['fastMASST'], axis=1), masst_df['fastMASST'].apply(pd.Series)], axis=1)
    return masst_df


def run_on_usi_and_id_list(input_file,
                           out_filename_no_ext,
                           usi_or_lib_id="Output USI",
                           compound_name_header="COMPOUND_NAME",
                           sep=",",
                           precursor_mz_tol=0.05,
                           mz_tol=0.02,
                           min_cos=0.7,
                           min_matched_signals=3,
                           analog=False,
                           analog_mass_below=150,
                           analog_mass_above=200,
                           parallel_queries=100
                           ):
    jobs_df = pd.read_csv(input_file, sep=sep)
    jobs_df.rename(columns={usi_or_lib_id: 'input_id', compound_name_header: 'Compound'}, inplace=True)
    jobs_df = jobs_df[jobs_df["input_id"].notnull()]
    jobs_df = jobs_df.astype({'Compound': 'string'})
    jobs_df["Compound"] = jobs_df["Compound"].apply(path_safe)

    logger.debug("Running fast microbe masst on input")

    with ThreadPoolExecutor(parallel_queries) as executor:
        futures = [executor.submit(query_usi_or_id,
                                   out_filename_no_ext,
                                   id,
                                   name,
                                   precursor_mz_tol=precursor_mz_tol,
                                   mz_tol=mz_tol,
                                   min_cos=min_cos,
                                   min_matched_signals=min_matched_signals,
                                   analog=analog,
                                   analog_mass_below=analog_mass_below,
                                   analog_mass_above=analog_mass_above
                                   ) for id, name in zip(jobs_df["input_id"], jobs_df["Compound"])]

        wait(futures)
        jobs_df["fastMASST"] = [f.result() for f in futures]

    save_masst_results(jobs_df, out_filename_no_ext)


def run_on_mgf(input_file,
               out_filename_no_ext,
               precursor_mz_tol=0.05,
               mz_tol=0.02,
               min_cos=0.7,
               min_matched_signals=3,
               analog=False,
               analog_mass_below=150,
               analog_mass_above=200,
               parallel_queries=100
               ):
    ids, precursor_mzs, precursor_charges = [], [], []
    mzs, intensities = [], []

    with pyteomics.mgf.MGF(input_file) as f_in:
        for spectrum_dict in tqdm(f_in):
            abundances = spectrum_dict["intensity array"]
            if len(abundances) >= min_matched_signals:
                ids.append(str(spectrum_dict["params"]["scans"]))
                precursor_mzs.append(float(spectrum_dict["params"]["pepmass"][0]))
                precursor_charges.append(int(spectrum_dict["params"]["charge"][0]))
                mzs.append(spectrum_dict["m/z array"])
                intensities.append(abundances)

    jobs_df = pd.DataFrame(
        {
            "Compound": ids,
            "precursor_mz": precursor_mzs,
            "precursor_charge": precursor_charges,
            "mzs": mzs,
            "intensities": intensities,
        }
    )

    logger.info("Running fast microbe masst on input n={} spectra".format(len(jobs_df)))

    with ThreadPoolExecutor(parallel_queries) as executor:
        futures = [executor.submit(query_spectrum,
                                   out_filename_no_ext,
                                   name,
                                   prec_mz,
                                   prec_charge,
                                   mz_array,
                                   intensity_array,
                                   precursor_mz_tol=precursor_mz_tol,
                                   mz_tol=mz_tol,
                                   min_cos=min_cos,
                                   min_matched_signals=min_matched_signals,
                                   analog=analog,
                                   analog_mass_below=analog_mass_below,
                                   analog_mass_above=analog_mass_above
                                   ) for name, prec_mz, prec_charge, mz_array, intensity_array in
                   zip(ids, precursor_mzs, precursor_charges, mzs, intensities)]

        wait(futures)
        jobs_df["fastMASST"] = [f.result() for f in futures]

    save_masst_results(jobs_df, out_filename_no_ext)


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
    parser.add_argument('--in_file', type=str,
                        help='input file either mgf with spectra or table that contains the two columns specified by '
                             'usi_or_lib_id and compound_header',
                        # default="../casmi_pos_sirius/bifido.mgf")
                        # default="../casmi_pos_sirius/small.mgf")
                        default="../examples/example_links.tsv")
    parser.add_argument('--out_file', type=str, help='output html and other files, name without extension',
                        default="output/fastMASST")

    # only for USI or lib ID file
    parser.add_argument('--usi_or_lib_id', type=str, help='specify the usi or GNPS library id to search',
                        default="USI")
    parser.add_argument('--compound_header', type=str,
                        help='defines the header of the compound names, make sure to be path safe',
                        default="Compound")
    parser.add_argument('--separator', type=str,
                        help='separator for input file, e.g., \\t for tab',
                        default="\t")

    # MASST params
    parser.add_argument('--precursor_mz_tol', type=float,
                        help='precursor mz tolerance',
                        default="0.05")
    parser.add_argument('--mz_tol', type=float,
                        help='mz tolerance to match signals',
                        default="0.02")
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
    parser.add_argument('--parallel_queries', type=int,
                        help='the number of async queries. fastMASST step is IO bound so higher number than CPU '
                             'speeds up the process',
                        default="50")

    args = parser.parse_args()

    try:
        in_file = args.in_file
        out_file = args.out_file
        lib_id = args.usi_or_lib_id
        compound_header = args.compound_header
        sep = args.separator

        # MASST parameters
        precursor_mz_tol = args.precursor_mz_tol
        mz_tol = args.mz_tol
        min_cos = args.min_cos
        min_matched_signals = args.min_matched_signals
        analog = args.analog
        analog_mass_below = args.analog_mass_below
        analog_mass_above = args.analog_mass_above
        parallel_queries = args.parallel_queries

        if str(in_file).endswith(".mgf"):
            run_on_mgf(input_file=in_file,
                       out_filename_no_ext=out_file,
                       precursor_mz_tol=precursor_mz_tol,
                       mz_tol=mz_tol,
                       min_cos=min_cos,
                       min_matched_signals=min_matched_signals,
                       analog=analog,
                       analog_mass_below=analog_mass_below,
                       analog_mass_above=analog_mass_above,
                       parallel_queries=parallel_queries
                       )
        else:
            run_on_usi_and_id_list(input_file=in_file,
                                   out_filename_no_ext=out_file,
                                   usi_or_lib_id=lib_id,
                                   compound_name_header=compound_header,
                                   sep=sep,
                                   precursor_mz_tol=precursor_mz_tol,
                                   mz_tol=mz_tol,
                                   min_cos=min_cos,
                                   min_matched_signals=min_matched_signals,
                                   analog=analog,
                                   analog_mass_below=analog_mass_below,
                                   analog_mass_above=analog_mass_above,
                                   parallel_queries=parallel_queries
                                   )
    except Exception as e:
        # exit with error
        logger.exception(e)
        sys.exit(1)

    # exit with OK
    sys.exit(0)
