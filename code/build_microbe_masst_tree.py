import sys
import argparse
import logging
from pathlib import Path

import bundle_to_html
import json_ontology_extender
import microbe_masst_results

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def create_tree_html(in_html="collapsible_tree_v3.html", in_ontology="../data/microbe_masst/ncbi.json",
                     metadata_file="../data/microbe_masst/microbe_masst_table.csv",
                     masst_file="../examples/phelylglycocholic_acid.tsv", matching_usi_list=None,
                     out_counts_file="dist/microbe_masst_counts.tsv",
                     out_json_tree="dist/merged_ncbi_ontology_data.json", format_out_json=True,
                     out_html="dist/oneindex.html", compress_out_html=True, node_key="NCBI", data_key="ncbi"):
    """
    Merges extra data into an ontology and creates a single distributable html file. Compression reduces the size of
    the html file.

    :param masst_file: imports masst results from file (classical masst). Is not used if a matching_usi_list is
    provided, e.g., from fastMASST
    :param matching_usi_list: matching universal spectrum identifier list that refer to masst matches in datasets.
    e.g., from fastMASST
    :param in_html: the base html file with different dependencies
    :param in_ontology: the input ontology
    :param out_json_tree: the merged tree data is exported to a json file
    :param out_html: the final distributable html file, merged with all dependencies
    :param compress_out_html: apply compression (reduces readability)
    """
    if out_counts_file == "auto" or out_counts_file == "automatic":
        out_counts_file = "dist/{}_counts.tsv".format(Path(masst_file).stem)

    if matching_usi_list is not None:
        microbe_masst_results.create_counts_file_from_usi(metadata_file, matching_usi_list, out_counts_file)
    else:
        microbe_masst_results.create_counts_file(metadata_file, masst_file, out_counts_file)

    json_ontology_extender.add_data_to_ontology_file(out_json_tree, in_ontology, out_counts_file, node_key, data_key,
                                                     format_out_json)
    return bundle_to_html.build_dist_html(in_html, out_html, out_json_tree, compress_out_html)


if __name__ == '__main__':
    # parsing the arguments (all optional)
    parser = argparse.ArgumentParser(description='Create tree data by merging extra data into an ontology. Then '
                                                 'create a distributable html file that internalizes all scripts, '
                                                 'data, etc. ')
    parser.add_argument('--in_html', type=str, help='The input html file',
                        default="collapsible_tree_v3.html")
    parser.add_argument('--ontology', type=str, help='the json ontology file with children',
                        default="../data/microbe_masst/ncbi.json")
    parser.add_argument('--metadata_file', type=str, help='microbe masst metadata',
                        default="../data/microbe_masst/microbe_masst_table.csv")
    parser.add_argument('--masst_file', type=str, help='a tab separated file with additional data that is added to '
                                                       'metadata file',
                        default="../examples/coprogen_b.tsv")
    parser.add_argument('--out_counts_file', type=str, help='the intermediate counts (matches) file. automatic: use '
                                                            'the masst_file name with suffix: _counts',
                        default="dist/microbe_masst_counts.tsv")
    parser.add_argument('--out_html', type=str, help='output html file', default="dist/microbeMasst_coprogen_b.html")
    parser.add_argument('--compress', type=bool, help='Compress output file (needs minify_html)',
                        default=True)
    parser.add_argument('--out_tree', type=str, help='output file', default="dist/merged_ncbi_ontology_data.json")
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
        create_tree_html(args.in_html, args.ontology, args.metadata_file, args.masst_file, args.out_counts_file,
                         args.out_tree, args.format, args.out_html, args.compress, args.node_key, args.data_key)
    except Exception as e:
        # exit with error
        logger.exception(e)
        sys.exit(1)

    # exit with OK
    sys.exit(0)
