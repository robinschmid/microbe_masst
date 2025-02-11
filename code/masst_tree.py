from pathlib import Path
import pandas as pd

from masst_utils import SpecialMasst
from masst_utils import SPECIAL_MASSTS
from utils import prepare_paths
import bundle_to_html
import json_ontology_extender
import logging
import json

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def create_enriched_masst_tree(
    matches_df,
    special_masst: SpecialMasst,
    common_file,
    lib_match_json,
    input_str,
    parameter_str,
    usi: str = None,
    in_html="../code/collapsible_tree_v3.html",
    format_out_json=False,
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
            "PARAMS_PLACEHOLDER": parameter_str,
        }

        prepare_paths(files=[out_counts_file, out_html, out_json_tree])

        # exports the counts file for all matches
        results_df = export_metadata_matches(special_masst, matches_df, out_counts_file)
        if len(results_df) <= 0:
            return None

        results_df = group_matches(special_masst, results_df)
        # adds them to the json ontology
        json_ontology_extender.add_data_to_ontology_file(
            special_masst=special_masst,
            output=out_json_tree,
            meta_matched_df=results_df,
            format_out_json=format_out_json,
        )
        # bundles the final html
        return bundle_to_html.build_dist_html(
            in_html, out_html, replace_dict, compress_out_html
        )
    except Exception as e:
        # exit with error
        logger.exception(e)
    # default return None
    return None


def create_combined_masst_tree(
    matches_df,
    common_file,
    lib_match_json,
    input_str,
    parameter_str,
    usi: str = None,
    in_html="../code/collapsible_tree_v3.html",
    format_out_json=False,
    compress_out_html=True,
):
    if (matches_df is None) or (len(matches_df) <= 0):
        return False

    tree_roots = []
    for special_masst in SPECIAL_MASSTS:
        try:
            out_json_tree = "{}_{}.json".format(common_file, special_masst.prefix)
            with open(out_json_tree) as json_file:
                # load all trees, add masst_type to identify, rename root
                treeRoot = json.load(json_file)
                json_ontology_extender.set_field_in_all_nodes(
                    treeRoot, "masst_type", special_masst.root
                )
                treeRoot["name"] = special_masst.root
                tree_roots.append(treeRoot)
        except:
            pass

    # skip if no or only one tree was detected
    if len(tree_roots) <= 1:
        return False

    combined_root = {"name": "root", "children": tree_roots}

    # add values to root
    json_ontology_extender.calc_root_stats(combined_root)
    json_ontology_extender.add_pie_data_to_node_and_children(combined_root, False)

    try:
        combined_prefix = "combined"
        out_json_tree = "{}_{}.json".format(common_file, combined_prefix)
        prepare_paths(files=[out_json_tree])

        with open(out_json_tree, "w") as file:
            if format_out_json:
                out_tree = json.dumps(
                    combined_root, indent=2, cls=json_ontology_extender.NpEncoder
                )
            else:
                out_tree = json.dumps(
                    combined_root, cls=json_ontology_extender.NpEncoder
                )
            print(out_tree, file=file)

        out_html = "{}_{}.html".format(common_file, combined_prefix)
        replace_dict = {
            "PLACEHOLDER_JSON_DATA": out_json_tree,
            "LIBRARY_JSON_DATA_PLACEHOLDER": lib_match_json,
            "INPUT_LABEL_PLACEHOLDER": input_str,
            "USI_LABEL_PLACEHOLDER": usi if usi else "",
            "PARAMS_PLACEHOLDER": parameter_str,
        }

        # bundles the final html
        return bundle_to_html.build_dist_html(
            in_html, out_html, replace_dict, compress_out_html
        )
    except Exception as e:
        # exit with error
        logger.exception(e)
    # default return None
    return None


def export_metadata_matches(
    special_masst: SpecialMasst, matches_df: pd.DataFrame, out_tsv_file
) -> pd.DataFrame:
    metadata_file = special_masst.metadata_file
    if str(metadata_file).endswith(".tsv"):
        metadata_df = pd.read_csv(metadata_file, sep="\t")
    else:
        metadata_df = pd.read_csv(metadata_file)

    # join on the file usi
    # results_df = matches_df.merge(metadata_df, on="file_usi", how="inner")
    results_df = pd.merge(matches_df, metadata_df, on="file_usi", how="inner")

    # export file with ncbi, matched_size,
    if len(results_df) > 0:
        results_df.to_csv(out_tsv_file, index=False, sep="\t")
    return results_df


def group_matches(special_masst: SpecialMasst, results_df) -> pd.DataFrame:
    grouped = results_df.groupby(special_masst.metadata_key)
    results_df = grouped.agg(matched_size=(special_masst.metadata_key, "size"))
    results_df["matches_json"] = grouped[["USI", "Cosine", "Matching Peaks"]].apply(
        lambda x: x.to_json(orient="records")
    )
    return results_df.reset_index()
