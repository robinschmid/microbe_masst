import logging
import glob, re
import json
import pandas as pd
from tqdm import tqdm
import masst_utils
from masst_utils import SpecialMasst

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def create_quant_summary(
        quant_csv,
        summary_df: pd.DataFrame,
        out_file=None,
        sum_as_binary_presence=False,
) -> pd.DataFrame:
    """
    Creates a table of samples (columns) and MASST matches (like NCBI id) as rows. The intensity in the quant table
    is row normalized to the relative intensity in the row. then multiplied and summed for each feature that had a
    match to a specific MASST metadata entry

    :param quant_csv: imports the quant table from MZmine
    :param summary_df: summary of MASST matches with the feature id as columns
    :param out_file: saves the results to csv
    :return: the final data frame
    """
    summary_df.columns = pd.to_numeric(summary_df.columns)
    quant_df = import_quantdf(quant_csv)
    quant_df = quant_df[quant_df.index.get_level_values("row ID").isin(summary_df.columns)]
    if sum_as_binary_presence:
        quant_df = quant_df.ge(0.01)

    samples = quant_df.columns
    final_df = pd.DataFrame(index=summary_df.index, columns=samples).fillna(0)

    summary_df.columns = pd.to_numeric(summary_df.columns)
    features = summary_df.columns
    for f in tqdm(features):
        try:
            # take first row
            _, quant_row = next(quant_df[quant_df.index.get_level_values("row ID") == f].iterrows())

            # add MASST matches to samples
            for sample in samples:
                final_df[sample] += summary_df[f] * quant_row[sample]
        except:
            logger.warning("This error should not happen, every row ID should be available in the table")
            pass

    if out_file:
        final_df.to_csv(out_file)
    return final_df


def import_quantdf(quant_csv):
    quant_df = pd.read_csv(quant_csv)
    samples = [col[:-10] for col in quant_df.columns if col.endswith(" Peak area")]
    quant_df.rename(columns=
                    dict([(col, col[:-10]) for col in quant_df.columns if col.endswith(" Peak area")]), inplace=True
                    )
    # ensure numeric columns, set index and retain only sample columns
    quant_df["row ID"] = pd.to_numeric(quant_df["row ID"])
    quant_df["row m/z"] = pd.to_numeric(quant_df["row m/z"])
    quant_df["row retention time"] = pd.to_numeric(quant_df["row retention time"])
    quant_df.set_index(["row ID", "row m/z", "row retention time"], inplace=True)
    quant_df = quant_df[samples]
    # relative intensities in rows
    quant_df = quant_df.div(quant_df.max(axis=1), axis=0) * 100
    quant_df.round(3)
    return quant_df


def create_summary_file(
        parent_directory,
        special_masst: SpecialMasst,
        out_file=None,
        min_matches=1,
        matches_to_binary_presence=True,
) -> pd.DataFrame | None:
    dfs = []
    node_id = special_masst.tree_node_key
    for file in tqdm(glob.glob(parent_directory + f"*{special_masst.prefix}.json")):
        comp_id = re.search(r"_(\d+)_"+special_masst.prefix, file).group(1)
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
        # if len(dfs) > 10:
        #     break

    if len(dfs) == 0:
        return None

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


def create_all_summary_files(special_masst:SpecialMasst, masst_directory, quant_csv, out_base_file, min_matches=1):
    out_base_file = "{}_{}".format(out_base_file, special_masst.prefix)
    merged_df = create_summary_file(
        parent_directory=masst_directory,
        out_file= "{}_spectral_matches.csv".format(out_base_file),
        special_masst=special_masst,
        min_matches=min_matches,
        matches_to_binary_presence=True,
    )
    if merged_df is None:
        return

    if quant_csv is not None:
        create_quant_summary(
            quant_csv=quant_csv,
            summary_df=merged_df,
            out_file=out_base_file + "_samples_binary.csv",
            sum_as_binary_presence=True,
        )

        create_quant_summary(
            quant_csv=quant_csv,
            summary_df=merged_df,
            out_file=out_base_file + "_samples_row_normalized.csv",
            sum_as_binary_presence=False,
        )

def create_all_masst_summaries(masst_directory, quant_csv, out_base_file, min_matches=1, ):
    for special_masst in masst_utils.SPECIAL_MASSTS:
        create_all_summary_files(special_masst, masst_directory, quant_csv, out_base_file, min_matches)


if __name__ == "__main__":
    create_all_masst_summaries(
        r"D:\git\microbe_masst\output\sydney\\",
        r"D:\git\microbe_masst\local_files\SEED_Grant_iimn_gnps_quant.csv",
        r"D:\git\microbe_masst\output\mmasst_summary_seed_grant",
        min_matches = 1,
    )