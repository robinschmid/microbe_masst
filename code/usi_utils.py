import json
from pathlib import Path
import requests

USI_URL = "https://metabolomics-usi.gnps2.org/json/"


def create_simple_file_usi(filename, dataset):
    filename = Path(filename).stem
    if not filename:
        raise ValueError("Filename is empty")
    if not dataset:
        raise ValueError("Dataset is empty")
    return "mzspec:{}:{}".format(dataset, filename)


def ensure_simple_file_usi(usi):
    # remove scan, check only dataset and filename
    scan = str(usi).rfind(":scan")
    if scan > -1:
        elements = usi[:scan].split(":")
    else:
        elements = usi.split(":")
    return "mzspec:{}:{}".format(elements[1], elements[-1])


def create_file_usi_column(
    df, original_usi_col="USI", dataset_col="MassIVE", filename_col="Filename"
):
    if original_usi_col in df.columns:
        df["file_usi"] = [ensure_simple_file_usi(usi) for usi in df[original_usi_col]]
    else:
        if filename_col not in df.columns:
            raise ValueError("Missing Filename column")
        if dataset_col not in df.columns:
            raise ValueError("Missing MassIVE column for datasets")
        df["file_usi"] = [
            create_simple_file_usi(filename, dataset)
            for filename, dataset in zip(df[filename_col], df[dataset_col])
        ]


def ensure_usi(usi_or_lib_id):
    """

    :param usi_or_lib_id: USI or gnps library ID
    :return: the input usi, the GNPS library usi, or None if not a USI
    """
    if not usi_or_lib_id:
        return None
    input = str(usi_or_lib_id)
    if input.startswith("mzspec:"):
        return input
    elif input.startswith("CCMSLIB"):
        return "mzspec:GNPS:GNPS-LIBRARY:accession:{}".format(input)
    else:
        return None


def get_spectrum(usi: str):
    resp = requests.get(USI_URL, params={"usi1": usi})
    resp.raise_for_status()
    return json.loads(resp.text)
