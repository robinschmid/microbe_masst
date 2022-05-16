import argparse
import sys
import pandas as pd
import json
import numpy as np
import logging

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


def add_data_to_node(node, df, node_field, data_field):
    """
    Merge data into node and apply to all children
    :param node: the current node in a tree structure with ["children"] property
    :param df: the data frame with additional data
    :param node_field: node[field] determines the key to align tree and additional data
    :param data_field: data[field] determines the key to align tree and additional data
    """
    try:
        ncbi = node.get(node_field)
        if ncbi is None:
            logger.warning("node has no id {}".format(node.get("name", "NONAME")))
        else:
            # use string for comparison of IDs
            filtered = df[df[data_field] == str(ncbi)]
            if len(filtered) > 0:
                rowi = filtered.index[0]
                for col, value in df.iteritems():
                    if col != data_field:
                        # print("%s is %s" % (col, value[rowi]))
                        node[col] = value[rowi]
    except Exception as ex:
        logger.exception(ex)
    # apply to all children
    if "children" in node:
        for child in node["children"]:
            add_data_to_node(child, df, node_field, data_field)


def accumulate_field_in_parents(node, field):
    """
    check count and group size fields and accumulate over tree
    :param field: field to propagate
    :param node: the current node in a tree structure with ["children"] property
    """
    node_value = node.get(field, 0)
    # apply to all children
    if "children" in node:
        for child in node["children"]:
            node_value += accumulate_field_in_parents(child, field)
    node[field] = node_value
    return node_value


def field_missing(node, field, report_missing=False, replace_with_field=None):
    """
    check if the node or any children down the tree lacks the field
    :param field: the field to be searched
    :param node: the current node in a tree structure with ["children"] property
    """
    missing = 0
    if node.get(field, None) is None:
        if report_missing:
            logger.error("Missing: {}".format(node.get("name", "NONAME")))
        if replace_with_field is not None:
            node[field] = node.get(replace_with_field, "")
        missing = 1
    if "children" in node:
        for child in node["children"]:
            missing += field_missing(child, field)

    return missing


def add_pie_data_to_node_and_children(node):
    # the pie data needs an array with multiple entries - therefore use fraction and 1-fraction
    node["pie_data"] = [{}, {}];
    node["pie_data"][0]["occurrence_fraction"] = node["occurrence_fraction"];
    node["pie_data"][0]["index"] = 0;
    node["pie_data"][0]["group_size"] = node["group_size"];
    node["pie_data"][0]["matched_size"] = node["matched_size"];
    node["pie_data"][1]["occurrence_fraction"] = 1.0 - node["occurrence_fraction"];
    node["pie_data"][1]["index"] = 1;
    node["pie_data"][1]["group_size"] = node["group_size"];
    node["pie_data"][1]["matched_size"] = node["matched_size"];

    # apply to all children
    if "children" in node:
        for child in node["children"]:
            add_pie_data_to_node_and_children(child)


def add_data_to_ontology_file(output="dist/merged_ontology_data.json", ontology_file="../data/GFOP.json",
                              in_data="../examples/caffeic_acid.tsv", node_key="name", data_key="group_value",
                              format_out_json=True):
    with open(ontology_file) as json_file:
        treeRoot = json.load(json_file)

        # read the additional data
        df = pd.read_csv(in_data, sep='\t')
        # ensure that the grouping columns are strings as we usually match string ids
        df[data_key] = df[data_key].astype(str)

        # print(df)

        # loop over all children
        add_data_to_node(treeRoot, df, node_key, data_key)

        # check if group_size is available otherwise propagate
        if field_missing(treeRoot, "NCBI", report_missing=True, replace_with_field="name") > 0:
            logger.error("NCBI id is missing in a node")
        if field_missing(treeRoot, "group_size") > 0:
            accumulate_field_in_parents(treeRoot, "group_size")
        if field_missing(treeRoot, "matched_size") > 0:
            accumulate_field_in_parents(treeRoot, "matched_size")

        calc_stats(treeRoot)

        # calc gfop specific data for root
        calc_root_stats(treeRoot)
        # add data in format for pie charts
        add_pie_data_to_node_and_children(treeRoot)

        print("Writing to {}".format(output))
        with open(output, "w") as file:
            if format_out_json:
                out_tree = json.dumps(treeRoot, indent=2, cls=NpEncoder)
            else:
                out_tree = json.dumps(treeRoot, cls=NpEncoder)
            print(out_tree, file=file)


def calc_stats(node):
    if "children" in node:
        for child in node["children"]:
            calc_stats(child)
    if node["group_size"] == 0:
        node["occurrence_fraction"] = 0
    else:
        node["occurrence_fraction"] = node["matched_size"] / node["group_size"]


def calc_root_stats(treeRoot):
    treeRoot["group_size"] = 0
    treeRoot["matched_size"] = 0
    for child in treeRoot["children"]:
        treeRoot["group_size"] += child["group_size"]
        treeRoot["matched_size"] += child["matched_size"]

    if treeRoot["group_size"] == 0:
        treeRoot["occurrence_fraction"] = 0
    else :
        treeRoot["occurrence_fraction"] = treeRoot["matched_size"] / treeRoot["group_size"]


if __name__ == '__main__':
    # parsing the arguments (all optional)
    parser = argparse.ArgumentParser(description='merge an ontology with external data')
    parser.add_argument('--ontology', type=str, help='the json ontology file with children',
                        default="../data/GFOP.json")
    parser.add_argument('--in_data', type=str, help='a tab separated file with additional data that is added to the '
                                                    'ontology', default="../examples/caffeic_acid.tsv")
    parser.add_argument('--node_key', type=str, help='the field in the ontology to be compare to the field in the '
                                                     'data file', default="name")
    parser.add_argument('--data_key', type=str,
                        help='the field in the data file to be compared to the field in the ontology',
                        default="group_value")
    parser.add_argument('--out_tree', type=str, help='output file', default="dist/merged_ontology_data.json")
    parser.add_argument('--format', type=bool, help='Format the json output False or True',
                        default=True)

    args = parser.parse_args()

    # is a url - try to download file
    # something like https://raw.githubusercontent.com/robinschmid/GFOPontology/master/data/GFOP.owl
    # important use raw file on github!
    try:
        add_data_to_ontology_file(output=args.out_tree, ontology_file=args.ontology, in_data=args.in_data,
                                  node_key=args.node_key, data_key=args.data_key, format_out_json=args.format)
    except Exception as e:
        # exit with error
        logger.exception(e)
        sys.exit(1)

    # exit with OK
    sys.exit(0)
