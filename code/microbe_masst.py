import sys
import argparse
import logging
from pathlib import Path

import masst_utils as masst
import build_microbe_masst_tree as mmtree

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def run_microbe_masst(usi_or_lib_id, precursor_mz_tol=0.05, mz_tol=0.02, min_cos=0.7,
                      in_html="collapsible_tree_v3.html", in_ontology="../data/microbe_masst/ncbi.json",
                      metadata_file="../data/microbe_masst/microbe_masst_table.csv",
                      out_counts_file="dist/microbe_masst_counts.tsv",
                      out_json_tree="dist/merged_ncbi_ontology_data.json", format_out_json=True,
                      out_html="dist/oneindex.html", compress_out_html=True, node_key="NCBI", data_key="ncbi"
                      ):
    try:
        matches = masst.fast_masst(usi_or_lib_id, precursor_mz_tol, mz_tol, min_cos)
        if (matches is not None) and (len(matches) > 0):
            match_usi_list = [match["USI"] for match in matches]
            mmtree.create_tree_html(in_html, in_ontology, metadata_file, None, match_usi_list, out_counts_file,
                                    out_json_tree, format_out_json, out_html, compress_out_html, node_key, data_key)
            return matches
    except Exception as e:
        # exit with error
        logger.exception(e)
    # default return None
    return None


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
                        default="collapsible_tree_v3.html")
    parser.add_argument('--ontology', type=str, help='the json ontology file with children',
                        default="../data/microbe_masst/ncbi.json")
    parser.add_argument('--metadata_file', type=str, help='microbe masst metadata',
                        default="../data/microbe_masst/microbe_masst_table.csv")
    parser.add_argument('--out_counts_file', type=str, help='the intermediate counts (matches) file. automatic: use '
                                                            'the masst_file name with suffix: _counts',
                        default="dist/microbe_masst_counts.tsv")
    parser.add_argument('--out_html', type=str, help='output html file', default="dist/microbeMasst.html")
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
