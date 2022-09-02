import sys
import logging
import pandas as pd
from urllib import parse

from pyparsing import col
from tqdm import tqdm
import re
from datetime import timedelta
import requests_cache
import argparse
import pyteomics.mgf

import microbe_masst as microbemasst

requests_cache.install_cache('fastmasst_cache', expire_after=timedelta(days=2))

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# activate pandas tqdm progress_apply
tqdm.pandas()

example_link = "https://robinschmid.github.io/microbe_masst/{}_{}.html"


def get_html_link(file_name, compound_name):
    return example_link.format(file_name, parse.quote(compound_name))


def run_job(file_name, usi_or_lib_id, compound_name):
    out_html = "../{}_{}.html".format(file_name, compound_name.replace(" ", "_"))
    out_json_tree = "../{}_{}.json".format(file_name, compound_name.replace(" ", "_"))

    result = microbemasst.run_microbe_masst(usi_or_lib_id, precursor_mz_tol=0.05, mz_tol=0.02, min_cos=0.7,
                                            in_html="../code/collapsible_tree_v3.html",
                                            in_ontology="../data/ncbi.json",
                                            metadata_file="../data/microbe_masst_table.csv",
                                            out_counts_file="../output/microbe_masst_counts.tsv",
                                            out_json_tree=out_json_tree,
                                            format_out_json=True,
                                            out_html=out_html, compress_out_html=True, node_key="NCBI", data_key="ncbi"
                                            )
    return result


def run_job_with_spectrum(file_name, compound_name, precursor_mz, precursor_charge, mzs, intensities):
    out_html = "../{}_{}.html".format(file_name, compound_name.replace(" ", "_"))
    out_json_tree = "../{}_{}.json".format(file_name, compound_name.replace(" ", "_"))

    matches = microbemasst.run_microbe_masst_for_spectrum(precursor_mz, precursor_charge, mzs, intensities,
                                                          precursor_mz_tol=0.05, mz_tol=0.02, min_cos=0.7,
                                                          in_html="../code/collapsible_tree_v3.html",
                                                          in_ontology="../data/ncbi.json",
                                                          metadata_file="../data/microbe_masst_table.csv",
                                                          out_counts_file="../output/microbe_masst_counts.tsv",
                                                          out_json_tree=out_json_tree,
                                                          format_out_json=True,
                                                          out_html=out_html, compress_out_html=True, node_key="NCBI",
                                                          data_key="ncbi"
                                                          )

    # food
    out_html = "../{}_{}.html".format(file_name, compound_name.replace(" ", "_")).replace("microbe", "food")
    out_json_tree = "../{}_{}.json".format(file_name, compound_name.replace(" ", "_")).replace("microbe", "food")
    matches = microbemasst.run_food_masst_for_matches(matches,
                                                      in_html="../code/collapsible_tree_v3.html",
                                                      out_json_tree=out_json_tree,
                                                      format_out_json=True,
                                                      out_html=out_html, compress_out_html=True,
                                                      )

    return matches


def path_safe(file):
    return re.sub('[^-a-zA-Z0-9_.() ]+', '_', file)


def explode_masst_columns(jobs_df):
    # multiple matches for each row split into multiple rows
    masst_df = jobs_df.explode("fastMASST", ignore_index=True)
    # each match is a dict so split into columns
    masst_df = pd.concat([masst_df.drop(['fastMASST'], axis=1), masst_df['fastMASST'].apply(pd.Series)], axis=1)
    return masst_df


def run_batch_fastmicrobe_masst(input_file, out_filename_no_ext, usi_or_lib_id="Output USI",
                                compound_name_header="COMPOUND_NAME", sep=","):
    finished_jobs_tsv = "../output/finished.tsv"
    masst_results_tsv = "../{}.tsv".format(out_filename_no_ext)

    jobs_df = pd.read_csv(input_file, sep=sep)
    jobs_df.rename(columns={usi_or_lib_id: 'input_id', compound_name_header: 'Compound'}, inplace=True)
    jobs_df = jobs_df[jobs_df["input_id"].notnull()]
    jobs_df = jobs_df.astype({'Compound': 'string'})
    jobs_df["Compound"] = jobs_df["Compound"].apply(path_safe)

    logger.debug("Running fast microbe masst on input")
    jobs_df["fastMASST"] = jobs_df.progress_apply(
        lambda row: run_job(out_filename_no_ext, row["input_id"], row["Compound"]), axis=1)

    # explode rows and export
    export_masst_results(jobs_df, masst_results_tsv)

    # add tree links
    jobs_df["Tree"] = jobs_df.progress_apply(
        lambda row: get_html_link(out_filename_no_ext, row["Compound"]) if row["fastMASST"] is not None else None,
        axis=1)
    # export the list of links
    logger.info("Exporting finished_jobs_tsv %s", finished_jobs_tsv)
    jobs_df.drop(columns=["fastMASST"], axis=1, inplace=True)
    jobs_df.to_csv(finished_jobs_tsv, sep="\t", index=False)


def run_on_mgf(input_file, out_filename_no_ext):
    finished_jobs_tsv = "../output/finished.tsv"
    masst_results_tsv = "../{}.tsv".format(out_filename_no_ext)

    ids, precursor_mzs, precursor_charges = [], [], []
    mzs, intensities = [], []

    with pyteomics.mgf.MGF(input_file) as f_in:
        for spectrum_dict in tqdm(f_in):
            ids.append(str(spectrum_dict["params"]["scans"]))
            precursor_mzs.append(float(spectrum_dict["params"]["pepmass"][0]))
            precursor_charges.append(int(spectrum_dict["params"]["charge"][0]))
            mzs.append(spectrum_dict["m/z array"])
            intensities.append(spectrum_dict["intensity array"])

    jobs_df = pd.DataFrame(
        {
            "Compound": ids,
            "precursor_mz": precursor_mzs,
            "precursor_charge": precursor_charges,
            "mzs": mzs,
            "intensities": intensities,
        }
    )

    logger.debug("Running fast microbe masst on input")
    jobs_df["fastMASST"] = jobs_df.progress_apply(
        lambda row: run_job_with_spectrum(out_filename_no_ext, row["Compound"], row["precursor_mz"],
                                          row["precursor_charge"], row["mzs"], row["intensities"]), axis=1)

    # explode rows and export
    export_masst_results(jobs_df, masst_results_tsv)

    # add tree links
    jobs_df["Tree"] = jobs_df.progress_apply(
        lambda row: get_html_link(out_filename_no_ext, row["Compound"]) if row["fastMASST"] is not None else None,
        axis=1)
    # export the list of links
    logger.info("Exporting finished_jobs_tsv %s", finished_jobs_tsv)
    jobs_df.drop(columns=["fastMASST"], axis=1, inplace=True)
    jobs_df.to_csv(finished_jobs_tsv, sep="\t", index=False)


def export_masst_results(jobs_df: pd.DataFrame, masst_results_tsv: str):
    # export masst results
    logger.info("Exporting masst results summary to %s", masst_results_tsv)
    masst_df = explode_masst_columns(jobs_df)
    # drop data points
    try:
        masst_df.drop(columns=["Query Scan", "Query Filename", "Index UnitPM", "Index IdxInUnitPM",
                               "Filtered Input Spectrum Path"], inplace=True, axis=1)
        # might not be in the df
        masst_df.drop(columns=["mzs", "intensities"], inplace=True, axis=1)
    except Exception as e:
        logger.exceptions("Error removing columns from export table", e)
    masst_df.to_parquet(masst_results_tsv.replace(".tsv", ".parquet"))


if __name__ == '__main__':
    # parsing the arguments (all optional)
    parser = argparse.ArgumentParser(description='Run fast microbeMASST in batch')
    parser.add_argument('--in_file', type=str,
                        help='input file either mgf with spectra or table that contains the two columns specified by '
                             'usi_or_lib_id and compound_header',
                        # default="../examples/alexey.csv")
    default="../casmi_pos_sirius/small.mgf")
    # default="../examples/example_links.tsv")
    parser.add_argument('--out_file', type=str, help='output html and other files, name without extension',
                        default="output/small_fast_microbeMasst")
    parser.add_argument('--usi_or_lib_id', type=str, help='specify the usi or GNPS library id to search',
                        default="USI")
    parser.add_argument('--compound_header', type=str,
                        help='defines the header of the compound names, make sure to be path safe',
                        default="Compound")
    parser.add_argument('--separator', type=str,
                        help='separator for input file, e.g., \\t for tab',
                        default=",")
    args = parser.parse_args()

    try:
        in_file = args.in_file
        out_file = args.out_file
        lib_id = args.usi_or_lib_id
        compound_header = args.compound_header
        sep = args.separator

        if str(in_file).endswith(".mgf"):
            run_on_mgf(input_file=in_file, out_filename_no_ext=out_file)
        else:
            run_batch_fastmicrobe_masst(input_file=in_file, out_filename_no_ext=out_file,
                                        usi_or_lib_id=lib_id, compound_name_header=compound_header,
                                        sep=sep)
    except Exception as e:
        # exit with error
        logger.exception(e)
        sys.exit(1)

    # exit with OK
    sys.exit(0)
