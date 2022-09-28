import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging
from usi_utils import ensure_simple_file_usi

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def clean_filename(name):
    return name.replace("/peak/", "/ccms_peak/").replace(".mzML", "").replace(".mzXML", "").replace("f.MSV", "MSV")


def filename_from_path(pathname):
    return Path(pathname).stem


def massive_id_from_path(pathname):
    path = Path(pathname.replace("f.MSV", "MSV"))
    for part in path.parts:
        if part.startswith("MSV"):
            massive_id = str(part)
            return massive_id


def usi_massiveid_filename(usi):
    """
    :param usi: input universal spectrum identifier
    :return: tuple(input usi, MassIVE ID, sample name without extension)
    """
    split = usi.split(":")
    return usi, split[1], split[2]


def create_counts_file(metadata_file, masst_file, out_tsv_file, meta_col_header="Taxa_NCBI", out_id_col_header='ncbi'):
    if str(metadata_file).endswith(".tsv"):
        metadata_df = pd.read_csv(metadata_file, sep="\t")
    else:
        metadata_df = pd.read_csv(metadata_file)
    masst_file = pd.read_csv(masst_file, sep="\t")

    metadata_df["Filepath"] = metadata_df["Filename"].apply(clean_filename)
    masst_file["filename"] = masst_file["filename"].apply(clean_filename)

    # count matches per ID
    id_matches_dict = dict()

    for index, match_row in masst_file.iterrows():
        # might have multiple rows in the metadata table if multiple IDs
        matching_metadata = metadata_df.loc[metadata_df["Filepath"] == match_row["filename"]]
        for index2, meta_row in matching_metadata.iterrows():
            ncbi_ = meta_row[meta_col_header]
            id_matches_dict[ncbi_] = id_matches_dict.get(ncbi_, 0) + 1

    export_ncbi_counts(id_matches_dict, out_tsv_file, out_id_col_header)


def create_counts_file_from_usi(metadata_file, matches_df: pd.DataFrame, out_tsv_file, meta_col_header="Taxa_NCBI",
                                out_id_col_header='ncbi'):
    if str(metadata_file).endswith(".tsv"):
        metadata_df = pd.read_csv(metadata_file, sep="\t")
    else:
        metadata_df = pd.read_csv(metadata_file)

    # create a usi column that only points to the dataset:file (not scan)
    matches_df["file_usi"] = [ensure_simple_file_usi(usi) for usi in matches_df["USI"]]
    matches_df = matches_df.sort_values(by=["Cosine", "Matching Peaks"], ascending=[False, False]).drop_duplicates(
        "file_usi")
    # join on the file usi
    results_df = pd.concat([matches_df.set_index('file_usi'), metadata_df.set_index('file_usi')], axis=1,
                           join='inner').reset_index()
    grouped = results_df.groupby("Taxa_NCBI")
    results_df = grouped.agg(matched_size=("Taxa_NCBI", "size"), taxa_name=("Taxaname_file", "first"))
    results_df["matches_json"] = grouped['USI', 'Cosine', 'Matching Peaks'].apply(lambda x: x.to_json(
        orient='records'))

    results_df = results_df.reset_index().rename(columns={'Taxa_NCBI': 'ncbi'})
    # export file with ncbi, matched_size,
    results_df.to_csv(out_tsv_file, index=False, sep="\t")


def export_ncbi_counts(id_matches_dict, out_tsv_file, out_id_col_header='ncbi'):
    df = pd.DataFrame().from_dict(id_matches_dict, orient='index', columns=['matched_size'])
    df.reset_index(inplace=True)
    df = df.rename(columns={'index': out_id_col_header})
    df.to_csv(out_tsv_file, index=False, sep="\t")


if __name__ == '__main__':
    # parsing the arguments (all optional)
    parser = argparse.ArgumentParser(description='Merge MASST results with microbeMASST metadata')
    parser.add_argument('--metadata_file', type=str, help='microbe masst metadata',
                        default="../data/microbe_masst_table.csv")
    parser.add_argument('--masst_file', type=str, help='a tab separated file with additional data that is added to '
                                                       'metadata file',
                        default="../examples/phelylglycocholic_acid.tsv")
    parser.add_argument('--out_tsv_file', type=str, help='output file in .tsv format',
                        default="../output/microbe_masst_counts.tsv")

    args = parser.parse_args()

    try:
        create_counts_file(metadata_file=args.metadata_file, masst_file=args.masst_file, out_tsv_file=args.out_tsv_file)
    except Exception as e:
        # exit with error
        logger.exception(e)
        sys.exit(1)

    # exit with OK
    sys.exit(0)
