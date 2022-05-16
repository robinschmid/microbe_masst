import sys
import logging
import csv
import pandas as pd
from urllib import parse
from tqdm import tqdm

import microbe_masst as micromasst

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# activate pandas tqdm progress_apply
tqdm.pandas()

example_link = "https://robinschmid.github.io/GFOPontology/{}_{}.html"


def run_job(file_name, usi_or_lib_id, compound_name):
    out_html = "../{}_{}.html".format(file_name, compound_name.replace(" ", "_"))

    result = micromasst.run_microbe_masst(usi_or_lib_id, precursor_mz_tol=0.05, mz_tol=0.02, min_cos=0.7,
                                          in_html="collapsible_tree_v3.html",
                                          in_ontology="../data/microbe_masst/ncbi.json",
                                          metadata_file="../data/microbe_masst/microbe_masst_table.csv",
                                          out_counts_file="dist/microbe_masst_counts.tsv",
                                          out_json_tree="dist/merged_ncbi_ontology_data.json", format_out_json=True,
                                          out_html=out_html, compress_out_html=True, node_key="NCBI", data_key="ncbi"
                                          )

    if result is not None:
        return example_link.format(file_name, parse.quote(compound_name))
    else:
        return "NO_SUCCESS"


# jobs = {
#         "CCMSLIB00005883671": "GABA",
#         "CCMSLIB00006582001": "phe_CA",
#         "CCMSLIB00006581985": "glu_CA",
#         "CCMSLIB00000006885": "surfactin_C13",
#         "CCMSLIB00000006895": "surfactin_C15",
#         "CCMSLIB00005721043": "salinosporamide_A",
#         "CCMSLIB00006694017": "salinomycin",
#         "CCMSLIB00000006871": "actinomycin_d",
#         "CCMSLIB00003134958": "actinomycin_d_nist14_match",
#         "CCMSLIB00000075066": "arylomycin_a4",
#         "CCMSLIB00005772087": "kanamycin_a",
#         "CCMSLIB00005464333": "PANTOTHENATE",
#         "CCMSLIB00005723628": "acyl_ferrioxamine_2",
#         "CCMSLIB00000072229": "PQS",
#         "CCMSLIB00000072137": "NHQ",
#         "CCMSLIB00000072038": "putative quinolone",
#         "CCMSLIB00003136275": "lovastatin",
#         "CCMSLIB00005435739": "lovastatin_2",
#         "CCMSLIB00004679270": "peptaibol",
#         "CCMSLIB00003134635": "azithromycin",
#         "CCMSLIB00000579271": "Surugamide A",
#         "CCMSLIB00000001621": "Desferrioxamine E",
#         "CCMSLIB00000072100": "Desferrioxamine B",
#         "CCMSLIB00005716848": "Acyl_Desferrioxamine_C13_Promicroferrioxamine",
#         "CCMSLIB00003739952": "Penicillin G",
#         "CCMSLIB00005435755": "Ferrichrome",
#         "CCMSLIB00000070253": "Barbamide",
#         "CCMSLIB00000001562": "Jamaicamide A",
#         "CCMSLIB00000075016": "Napsamycin",
#         "CCMSLIB00001059079": "triacylfusarin",
#         "CCMSLIB00000074975": "Fusarin A",
#         "CCMSLIB00004721498": "Polymyxin B",
#         "CCMSLIB00005732667": "Microcystin LR",
#         "CCMSLIB00000005120": "Mycophenolic",
#         "CCMSLIB00005755738": "Mycotoxin F2",
#         "CCMSLIB00005723573": "Beauvericin",
#         "CCMSLIB00005727552": "Beauvericin_2"
#     }

if __name__ == '__main__':
    # jobs = {k: v.replace(" ", "_") for (k, v) in jobs.items()}
    file_name ="emily/fast_microbeMasst"
    finished_jobs_tsv = "../emily/example_links.tsv"
    try:
        finsihed_jobs_df = pd.read_csv(finished_jobs_tsv, sep="\t")
    except:
        logger.info("No example list or error during read")
    # TODO check which jobs already finished

    # read list of jobs
    # jobs_df = pd.read_csv("../emily/Combinatorial_reactions_USIs - fatty acid amides.tsv", sep="\t")
    jobs_df = pd.read_csv("../emily/Combinatorial_reactions_USIs - sulfated compounds.tsv", sep="\t")
    jobs_df.rename(columns={'Output USI': 'ID', 'COMPOUND_NAME': 'Compound'}, inplace=True)

    jobs_df["Tree"] = jobs_df.progress_apply(lambda row: run_job(file_name, row["ID"], row["Compound"]), axis=1)

    jobs_df.to_csv(finished_jobs_tsv, sep="\t")

    sys.exit(0)
