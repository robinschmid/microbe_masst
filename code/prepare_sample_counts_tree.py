import sys
import argparse
import logging
import pandas as pd
import json
import numpy as np

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)


def update_group_size(node, metadata_df, node_key="NCBI", data_key="Taxa_NCBI"):
    group_size = (metadata_df[data_key].values == node[node_key]).sum()
    # apply to all children
    if "children" in node:
        for child in node["children"]:
            group_size += update_group_size(child, metadata_df, node_key, data_key)
    node["group_size"] = group_size
    return group_size


def get_all_ids(ncbi_ids, node, node_key="NCBI"):
    ncbi_ids.append(node[node_key])
    # apply to all children
    if "children" in node:
        for child in node["children"]:
            get_all_ids(ncbi_ids, child, node_key)
    return ncbi_ids


def update_metadata_on_tree(
    in_ontology="../data/microbe_masst_tree.json",
    metadata_file="../data/microbe_masst_table.csv",
    node_key="NCBI",
    data_key="Taxa_NCBI",
):
    try:
        with open(in_ontology) as json_file:
            treeRoot = json.load(json_file)
            df = pd.read_csv(metadata_file, sep=",")
            # ensure that the grouping columns are strings as we usually match string ids
            df[data_key] = df[data_key].astype(str)

            update_group_size(treeRoot, df, node_key, data_key)
            ids = get_all_ids(list(), treeRoot, node_key)
            not_in_tree_df = df[~df[data_key].isin(ids)]
            not_in_tree_df.to_csv("../data/not_in_tree.csv", index=False)

        with open(in_ontology, "w") as file:
            out_tree = json.dumps(treeRoot, indent=2, cls=NpEncoder)
            print(out_tree, file=file)
        return True
    except Exception as e:
        logger.exception(e)
        return False


if __name__ == "__main__":
    # parsing the arguments (all optional)
    parser = argparse.ArgumentParser(
        description="Update tree ontology with new metadata file"
    )

    # microbemasst
    #parser.add_argument("--ontology", type=str, help="the json ontology file with children",
    #    default="../data/microbe_masst_tree.json")
    #parser.add_argument("--metadata_file", type=str, help="microbe masst metadata",
    #    default="../data/microbe_masst_table.csv")
    #parser.add_argument( "--node_key", type=str, help="the field in the ontology to be compare to the field in the " "data file",
    #    default="NCBI")
    #parser.add_argument("--data_key", type=str, help="the field in the data file to be compared to the field in the ontology",
    #    default="Taxa_NCBI")

    # plantmasst
    #parser.add_argument('--ontology', type=str, help='the json ontology file with children',
    #                    default="../data/plant_masst_tree.json")
    #parser.add_argument('--metadata_file', type=str, help='masst file metadata',
    #                    default="../data/plant_masst_table.csv")
    #parser.add_argument('--node_key', type=str, help='the field in the ontology to be compare to the field in the '
    #                                                 'data file', default="NCBI")
    #parser.add_argument('--data_key', type=str,
    #                    help='the field in the data file to be compared to the field in the ontology',
    #                   default="Taxa_NCBI")

    # foodmasst
    # parser.add_argument("--ontology", type=str, help="the json ontology file with children",
    #     default="../data/food_masst_tree.json")
    # parser.add_argument("--metadata_file", type=str, help="microbe masst metadata",
    #     default="../data/food_masst_table.csv")
    # parser.add_argument("--node_key", type=str, help="the field in the ontology to be compare to the field in the " "data file",
    #     default="name")
    # parser.add_argument("--data_key", type=str, help="the field in the data file to be compared to the field in the ontology",
    #     default="node_id")

    # tissuemasst
    parser.add_argument('--ontology', type=str, help='the json ontology file with children',
                        default="../data/tissue_masst_tree.json")
    parser.add_argument('--metadata_file', type=str, help='masst file metadata',
                        default="../data/tissue_masst_table.csv")
    parser.add_argument('--node_key', type=str, help='the field in the ontology to be compare to the field in the '
                                                     'data file', default="ID")
    parser.add_argument('--data_key', type=str,
                        help='the field in the data file to be compared to the field in the ontology',
                        default="ID")


    args = parser.parse_args()

    try:
        update_metadata_on_tree(
            args.ontology, args.metadata_file, args.node_key, args.data_key
        )
    except Exception as e:
        # exit with error
        logger.exception(e)
        sys.exit(1)

    # exit with OK
    sys.exit(0)