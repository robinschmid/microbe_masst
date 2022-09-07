import sys
import argparse
import logging
from pathlib import Path
import masst_utils as masst
import build_microbe_masst_tree as mmtree
import build_food_masst_tree as foodtree

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def run_microbe_masst(usi_or_lib_id, precursor_mz_tol=0.05, mz_tol=0.02, min_cos=0.7,
                      in_html="../code/collapsible_tree_v3.html", in_ontology="../data/ncbi.json",
                      metadata_file="../data/microbe_masst_table.csv",
                      out_counts_file="../output/microbe_masst_counts.tsv",
                      out_json_tree="../output/merged_ncbi_ontology_data.json", format_out_json=True,
                      out_html="../output/oneindex.html", compress_out_html=True, node_key="NCBI", data_key="ncbi"
                      ):
    prepare_paths(out_counts_file, out_html, out_json_tree)

    try:
        matches = masst.fast_masst(usi_or_lib_id, precursor_mz_tol, mz_tol, min_cos)
        return run_microbe_masst_for_matches(matches, in_html, in_ontology, metadata_file, out_counts_file,
                                             out_json_tree, format_out_json, out_html, compress_out_html, node_key,
                                             data_key)
    except Exception as e:
        # exit with error
        logger.exception(e)
    # default return None
    return None


def run_microbe_masst_for_spectrum(precursor_mz, precursor_charge, mzs, intensities, precursor_mz_tol=0.05, mz_tol=0.02,
                                   min_cos=0.7,
                                   in_html="../code/collapsible_tree_v3.html", in_ontology="../data/ncbi.json",
                                   metadata_file="../data/microbe_masst_table.csv",
                                   out_counts_file="../output/microbe_masst_counts.tsv",
                                   out_json_tree="../output/merged_ncbi_ontology_data.json", format_out_json=True,
                                   out_html="../output/oneindex.html", compress_out_html=True, node_key="NCBI",
                                   data_key="ncbi"
                                   ):
    prepare_paths(out_counts_file, out_html, out_json_tree)

    try:
        matches, _ = masst.fast_masst_spectrum(mzs, intensities, precursor_mz, precursor_charge, precursor_mz_tol,
                                               mz_tol,
                                               min_cos)
        return run_microbe_masst_for_matches(matches, in_html, in_ontology, metadata_file, out_counts_file,
                                             out_json_tree, format_out_json, out_html, compress_out_html, node_key,
                                             data_key)
    except Exception as e:
        # exit with error
        logger.exception(e)
    # default return None
    return None


def run_microbe_masst_for_matches(masst_matches,
                                  in_html="../code/collapsible_tree_v3.html", in_ontology="../data/ncbi.json",
                                  metadata_file="../data/microbe_masst_table.csv",
                                  out_counts_file="../output/microbe_masst_counts.tsv",
                                  out_json_tree="../output/merged_ncbi_ontology_data.json", format_out_json=True,
                                  out_html="../output/oneindex.html", compress_out_html=True, node_key="NCBI",
                                  data_key="ncbi",
                                  replace_dict=None
                                  ):
    prepare_paths(out_counts_file, out_html, out_json_tree)

    try:
        if (masst_matches is not None) and (len(masst_matches) > 0):
            mmtree.create_tree_html(in_html, in_ontology, metadata_file, None, masst_matches["USI"], out_counts_file,
                                    out_json_tree, format_out_json, out_html, compress_out_html, node_key, data_key,
                                    replace_dict=replace_dict)
            return masst_matches
    except Exception as e:
        # exit with error
        logger.exception(e)
    # default return None
    return None


def run_food_masst_for_matches(masst_matches,
                               in_html="../code/collapsible_tree_v3.html", in_ontology="../data/GFOP.json",
                               metadata_file="../data/foodmasst_filtered.tsv",
                               out_counts_file="../output/food_masst_counts.tsv",
                               out_json_tree="../output/merged_gfop_ontology_data.json", format_out_json=True,
                               out_html="../output/oneindex.html", compress_out_html=True, node_key="name",
                               data_key="ontology_term",
                               replace_dict=None
                               ):
    prepare_paths(out_counts_file, out_html, out_json_tree)

    try:
        if (masst_matches is not None) and (len(masst_matches) > 0):
            foodtree.create_tree_html(in_html, in_ontology, metadata_file, None, masst_matches["USI"], out_counts_file,
                                      out_json_tree, format_out_json, out_html, compress_out_html, node_key, data_key,
                                      replace_dict=replace_dict)
            return masst_matches
    except Exception as e:
        # exit with error
        logger.exception(e)
    # default return None
    return None


def prepare_paths(out_counts_file, out_html, out_json_tree):
    try:
        # ensure output paths
        Path(out_html).parent.mkdir(parents=True, exist_ok=True)
        Path(out_counts_file).parent.mkdir(parents=True, exist_ok=True)
        Path(out_json_tree).parent.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.exception(e)


if __name__ == '__main__':
    # parsing the arguments (all optional)
    parser = argparse.ArgumentParser(description='Create tree data by merging extra data into an ontology. Then '
                                                 'create a distributable html file that internalizes all scripts, '
                                                 'data, etc. ')
    parser.add_argument('--usi_or_lib_id', type=str,
                        help='universal spectrum identifier or GNPS library ID to search by fastMASST',
                        default="mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005883671")
    # default="mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00000001556")
    parser.add_argument('--in_html', type=str, help='The input html file',
                        default="../code/collapsible_tree_v3.html")
    parser.add_argument('--ontology', type=str, help='the json ontology file with children',
                        default="../data/ncbi.json")
    parser.add_argument('--metadata_file', type=str, help='microbe masst metadata',
                        default="../data/microbe_masst_table.csv")
    parser.add_argument('--out_counts_file', type=str, help='the intermediate counts (matches) file. automatic: use '
                                                            'the masst_file name with suffix: _counts',
                        default="../output/microbe_masst_counts.tsv")
    parser.add_argument('--out_html', type=str, help='output html file', default="../output/microbeMasst.html")
    parser.add_argument('--compress', type=bool, help='Compress output file (needs minify_html)',
                        default=True)
    parser.add_argument('--out_tree', type=str, help='output file', default="../output/merged_ncbi_ontology_data.json")
    parser.add_argument('--format', type=bool, help='Format the json output False or True',
                        default=True)
    parser.add_argument('--node_key', type=str, help='the field in the ontology to be compare to the field in the '
                                                     'data file', default="NCBI")
    parser.add_argument('--data_key', type=str,
                        help='the field in the data file to be compared to the field in the ontology',
                        default="ncbi")
    args = parser.parse_args()

    # is a url - try to download file
    # something like https://raw.githubusercontent.com/robinschmid/GFOPontology/master/data/GFOP.owl
    # important use raw file on github!
    try:
        run_microbe_masst(args.usi_or_lib_id, 0.05, 0.02, 0.7,
                          # tree generation
                          args.in_html, args.ontology, args.metadata_file, args.out_counts_file,
                          args.out_tree, args.format, args.out_html, args.compress, args.node_key, args.data_key)
    except Exception as e:
        # exit with error
        logger.exception(e)
        sys.exit(1)

    # exit with OK
    sys.exit(0)
