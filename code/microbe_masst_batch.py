import sys
import logging
import csv
import pandas as pd
from urllib import parse
from tqdm import tqdm
import re

import microbe_masst as microbemasst

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# activate pandas tqdm progress_apply
tqdm.pandas()

example_link = "https://robinschmid.github.io/microbe_masst/{}_{}.html"


def get_html_link(file_name, compound_name):
    return example_link.format(file_name, parse.quote(compound_name))


def run_job(file_name, usi_or_lib_id, compound_name):
    out_html = "../{}_{}.html".format(file_name, compound_name.replace(" ", "_"))

    result = microbemasst.run_microbe_masst(usi_or_lib_id, precursor_mz_tol=0.05, mz_tol=0.02, min_cos=0.7,
                                          in_html="collapsible_tree_v3.html",
                                          in_ontology="../data/ncbi.json",
                                          metadata_file="../data/microbe_masst_table.csv",
                                          out_counts_file="../output/microbe_masst_counts.tsv",
                                          out_json_tree="../output/merged_ncbi_ontology_data.json",
                                          format_out_json=True,
                                          out_html=out_html, compress_out_html=True, node_key="NCBI", data_key="ncbi"
                                          )
    return result


def path_safe(file):
    return re.sub('[^-a-zA-Z0-9_.() ]+', '_', file)


def explode_masst_columns(jobs_df):
    masst_df = jobs_df.explode("fastMASST", ignore_index=True)
    masst_df = masst_df.join(pd.DataFrame([masst for masst in masst_df['fastMASST']]))
    masst_df.drop(columns=["fastMASST"], inplace=True)
    return masst_df

if __name__ == '__main__':
    input_file = "../examples/emily.csv"
    # output
    file_name = "emily/fast_microbeMasst"
    finished_jobs_tsv = "../emily/finished.tsv"
    masst_results_tsv = "../{}.tsv".format(file_name)

    try:
        finsihed_jobs_df = pd.read_csv(finished_jobs_tsv, sep="\t")
    except:
        logger.info("No example list or error during read")
    # TODO check which jobs already finished

    # read list of jobs
    # jobs_df = pd.read_csv("../emily/Combinatorial_reactions_USIs - fatty acid amides.tsv", sep="\t")
    jobs_df = pd.read_csv(input_file, sep=",")
    jobs_df.rename(columns={'USI': 'input_id', 'COMPOUND_NAME': 'Compound'}, inplace=True)
    # jobs_df.rename(columns={'Output USI': 'ID', 'COMPOUND_NAME': 'Compound'}, inplace=True)
    jobs_df["Compound"] = jobs_df["Compound"].apply(path_safe)

    jobs_df["fastMASST"] = jobs_df.progress_apply(lambda row: run_job(file_name, row["input_id"], row["Compound"]), axis=1)
    jobs_df["Tree"] = jobs_df.progress_apply(lambda row:
                                             get_html_link(file_name, row["Compound"]) if row["fastMASST"] is not
                                                                                          None else None,
                                             axis=1)

    # export masst results
    masst_df = explode_masst_columns(jobs_df)
    masst_df.to_csv(masst_results_tsv, sep="\t", index=False)

    # export the list of links
    jobs_df.drop(columns=["fastMASST"], inplace=True)
    jobs_df.to_csv(finished_jobs_tsv, sep="\t", index=False)

    sys.exit(0)
