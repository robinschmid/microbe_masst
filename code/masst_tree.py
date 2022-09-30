import logging
from pathlib import Path
import pandas as pd

from masst_utils import SpecialMasst
import bundle_to_html
import json_ontology_extender

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def create_enriched_masst_tree(matches_df,
                               special_masst: SpecialMasst,
                               common_file,
                               lib_match_json,
                               input_str,
                               parameter_str,
                               usi: str = None,
                               in_html="../code/collapsible_tree_v3.html",
                               format_out_json=True,
                               compress_out_html=True,
                               ):
    if (matches_df is None) or (len(matches_df) <= 0):
        return False

    try:
        out_html = "{}_{}.html".format(common_file, special_masst.prefix)
        out_json_tree = "{}_{}.json".format(common_file, special_masst.prefix)
        out_counts_file = "{}_counts_{}.tsv".format(common_file, special_masst.prefix)
        replace_dict = {
            "PLACEHOLDER_JSON_DATA": out_json_tree,
            "LIBRARY_JSON_DATA_PLACEHOLDER": lib_match_json,
            "INPUT_LABEL_PLACEHOLDER": input_str,
            "USI_LABEL_PLACEHOLDER": usi if usi else "",
            "PARAMS_PLACEHOLDER": parameter_str
        }

        prepare_paths(out_counts_file, out_html, out_json_tree)

        # exports the counts file for all matches
        results_df = export_metadata_matches(special_masst, matches_df, out_counts_file)
        results_df = group_matches(special_masst, results_df)
        # adds them to the json ontology
        json_ontology_extender.add_data_to_ontology_file(special_masst=special_masst,
                                                         output=out_json_tree,
                                                         meta_matched_df=results_df,
                                                         format_out_json=format_out_json
                                                         )
        # bundles the final html
        return bundle_to_html.build_dist_html(in_html, out_html, replace_dict, compress_out_html)
    except Exception as e:
        # exit with error
        logger.exception(e)
    # default return None
    return None


def export_metadata_matches(special_masst: SpecialMasst, matches_df: pd.DataFrame, out_tsv_file) -> pd.DataFrame:
    metadata_file = special_masst.metadata_file
    if str(metadata_file).endswith(".tsv"):
        metadata_df = pd.read_csv(metadata_file, sep="\t")
    else:
        metadata_df = pd.read_csv(metadata_file)

    # join on the file usi
    results_df = pd.concat([matches_df.set_index('file_usi'), metadata_df.set_index('file_usi')], axis=1,
                           join='inner').reset_index()

    # export file with ncbi, matched_size,
    results_df.to_csv(out_tsv_file, index=False, sep="\t")
    return results_df


def group_matches(special_masst: SpecialMasst, results_df) -> pd.DataFrame:
    grouped = results_df.groupby(special_masst.metadata_key)
    results_df = grouped.agg(matched_size=(special_masst.metadata_key, "size"))
    results_df["matches_json"] = grouped[['USI', 'Cosine', 'Matching Peaks']].apply(lambda x: x.to_json(
        orient='records'))
    return results_df.reset_index()


def prepare_paths(out_counts_file, out_html, out_json_tree):
    try:
        # ensure output paths
        Path(out_html).parent.mkdir(parents=True, exist_ok=True)
        Path(out_counts_file).parent.mkdir(parents=True, exist_ok=True)
        Path(out_json_tree).parent.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.exception(e)
