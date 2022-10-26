import logging
import glob, re
import json
import pandas as pd
from tqdm import tqdm
import masst_utils
from masst_utils import SpecialMasst

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def create_summary_file(
    parent_directory,
    special_masst: SpecialMasst,
    out_file=None,
    min_matches=1,
    matches_to_binary_presence=True,
) -> pd.DataFrame:
    dfs = []
    node_id = special_masst.tree_node_key
    for file in tqdm(glob.glob(parent_directory + f"*{special_masst.prefix}.json")):
        comp_id = re.search(r"_(CCMSLIB\d+)_", file).group(1)
        df = json_to_dataframe(file, node_key=node_id, min_matches=min_matches)
        if df is None:
            continue

        df[node_id] = df[node_id].astype(str)
        df = create_multi_index(df.drop_duplicates(node_id), special_masst)
        # set detected to 1
        if matches_to_binary_presence:
            df[comp_id] = 1
        else:
            df[comp_id] = df["matched_size"]

        df = df.drop(columns=["matched_size"])
        if len(df) > 0:
            dfs.append(df)

    merged_df = pd.concat(dfs, join="outer", axis=1).fillna(0)
    if out_file:
        merged_df.to_csv(out_file)
    return merged_df


def for_all_children(rows, node, minimum_matches=1, node_key="NCBI", level=1):
    matched = node.get("matched_size", 0)
    if matched >= minimum_matches:
        rows.append(
            {
                "Name": node.get("name", None),
                "Rank": node.get("Rank", None),
                "Level": level,
                node_key: str(node[node_key]),
                "group_size": node["group_size"],
                "matched_size": matched,
            }
        )
        # apply to all children
        if "children" in node:
            for child in node["children"]:
                for_all_children(rows, child, minimum_matches, node_key, level + 1)


def json_to_dataframe(file, node_key="NCBI", min_matches=1):
    rows = []
    with open(file) as json_file:
        treeRoot = json.load(json_file)
        for node in treeRoot["children"]:
            for_all_children(rows, node, min_matches, node_key=node_key)
    if len(rows) == 0:
        return None
    return pd.DataFrame(rows).sort_values(
        ["Rank", "matched_size", "group_size"], ascending=False
    )


def create_multi_index(df: pd.DataFrame, special_masst: SpecialMasst) -> pd.DataFrame:
    return df.set_index(
        [special_masst.tree_node_key, "Name", "Rank", "group_size", "Level"]
    )


if __name__ == "__main__":
    merged_df = create_summary_file(
        parent_directory=r"D:\\git\\microbe_masst\\output\\library\\",
        out_file=r"D:\git\microbe_masst\output\summary_gnps_lib2.csv",
        special_masst=masst_utils.MICROBE_MASST,
        min_matches=2,
        matches_to_binary_presence=True,
    )
