import sys
import argparse
import logging
import pandas as pd
from usi_utils import create_file_usi_column

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def prepare_check_metadata_file(metadata_file, output_file):
    if metadata_file.endswith(".tsv"):
        df = pd.read_csv(metadata_file, sep="\t")
    else:
        df = pd.read_csv(metadata_file, sep=",")

    # remove white space around values in columns
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    if "file_usi" not in df.columns:
        create_file_usi_column(df)

    # check uniqueness of files
    logger.info("Checking duplicated usi")
    # microbe masst specific
    if "Taxa_NCBI" in df.columns:
        duplicates = df[df.duplicated(["file_usi"], keep=False)].sort_values(
            by=["file_usi", "Taxa_Assigment", "Taxa_NCBI"],
            ascending=[True, False, False],
        )
    else:
        duplicates = df[df.duplicated(["file_usi"], keep=False)].sort_values(
            by=["file_usi"],
            ascending=[True],
        )

    try:
        if len(duplicates) > 0:
            logger.info(duplicates)
            duplicates.to_csv("../data/duplicates.csv", index=False)
        else:
            logger.info("NO DUPLICATES")
    except:
        logger.warning("Cannot export duplicates file")

    df = df.sort_values(
        by=["file_usi", "Taxa_Assigment", "Taxa_NCBI"],
        ascending=[True, False, False],
    ).drop_duplicates(["file_usi"])

    # export final table
    if output_file.endswith(".tsv"):
        df.to_csv(output_file, index=False, sep="\t")
    else:
        df.to_csv(output_file, index=False)


if __name__ == "__main__":
    # parsing the arguments (all optional)
    parser = argparse.ArgumentParser(description="Update metadata file and check")
    parser.add_argument(
        "--metadata_file",
        type=str,
        help="input masst metadata",
        default="../data/table_microbeMASST_manula_fix_final.csv",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        help="output masst metadata",
        default="../data/mmast_masst_metadata.csv",
    )
    # foodmasst
    # parser.add_argument('--metadata_file', type=str, help='input masst metadata',
    #                     default="../data/food_masst_metadata.csv")
    # parser.add_argument('--output_file', type=str, help='output masst metadata',
    #                     default="../data/food_masst_metadata.csv")
    args = parser.parse_args()

    try:
        prepare_check_metadata_file(args.metadata_file, args.output_file)
    except Exception as e:
        # exit with error
        logger.exception(e)
        sys.exit(1)

    # exit with OK
    sys.exit(0)
