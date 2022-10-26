import sys
import logging
import pandas as pd
from tqdm import tqdm
import re
import argparse
import pyteomics.mgf
from pathlib import Path

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait

import masst_client

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# activate pandas tqdm progress_apply
tqdm.pandas()


def path_safe(file):
    return re.sub("[^-a-zA-Z0-9_.() ]+", "_", file)


def run_on_usi_list_or_mgf_file(
    in_file,
    out_file_no_extension="../output/fastMASST",
    # only for USI list
    usi_or_lib_id="USI",
    compound_name_header="Compound",
    sep=",",
    # general matching parameters
    precursor_mz_tol=0.05,
    mz_tol=0.05,
    min_cos=0.7,
    min_matched_signals=3,
    analog: bool = False,
    analog_mass_below=150,
    analog_mass_above=200,
    parallel_queries=10,
    skip_existing=False,
):
    """

    :param in_file: the input file either tabular USI,Compound or mgf with spectra
    :param out_file_no_extension: the output file prefix - no format extension
    :param usi_or_lib_id: header for USI or GNPS library ID (not for mgf)
    :param compound_name_header: compound name header (not for mgf)
    :param sep: separator (not for mgf)
    :param precursor_mz_tol: matching precursor mz tolerance in Da
    :param mz_tol: the signal mz tolerance in Da
    :param min_cos: minimum cosine score
    :param min_matched_signals: minimum matched signals
    :param analog: search analogs bool
    :param analog_mass_below: analog search window below precursor mz
    :param analog_mass_above: analog search window above precursor mz
    :param parallel_queries: perform queries in parallel
    :param skip_existing: skip existing files
    :return: success rate between 0-1 (skipped existing files excluded)
    """
    if str(in_file).endswith(".mgf"):
        return run_on_mgf(
            input_file=in_file,
            out_filename_no_ext=out_file_no_extension,
            precursor_mz_tol=precursor_mz_tol,
            mz_tol=mz_tol,
            min_cos=min_cos,
            min_matched_signals=min_matched_signals,
            analog=analog,
            analog_mass_below=analog_mass_below,
            analog_mass_above=analog_mass_above,
            parallel_queries=parallel_queries,
            skip_existing=skip_existing,
        )
    else:
        return run_on_usi_and_id_list(
            input_file=in_file,
            out_filename_no_ext=out_file_no_extension,
            usi_or_lib_id=usi_or_lib_id,
            compound_name_header=compound_name_header,
            sep=sep,
            precursor_mz_tol=precursor_mz_tol,
            mz_tol=mz_tol,
            min_cos=min_cos,
            min_matched_signals=min_matched_signals,
            analog=analog,
            analog_mass_below=analog_mass_below,
            analog_mass_above=analog_mass_above,
            parallel_queries=parallel_queries,
            skip_existing=skip_existing,
        )


def run_on_usi_and_id_list(
        input_file,
        out_filename_no_ext,
        usi_or_lib_id="Output USI",
        compound_name_header="COMPOUND_NAME",
        sep=",",
        precursor_mz_tol=0.05,
        mz_tol=0.02,
        min_cos=0.7,
        min_matched_signals=3,
        analog=False,
        analog_mass_below=150,
        analog_mass_above=200,
        parallel_queries=100,
        skip_existing=False,
):
    jobs_df = pd.read_csv(input_file, sep=sep)
    jobs_df.rename(
        columns={usi_or_lib_id: "input_id", compound_name_header: "Compound"},
        inplace=True,
    )
    jobs_df = jobs_df[jobs_df["input_id"].notnull()]
    jobs_df = jobs_df.astype({"Compound": "string"})
    jobs_df["Compound"] = jobs_df["Compound"].apply(path_safe)

    if skip_existing:
        all_len = len(jobs_df)
        jobs_df["file_path"] = [
            "{}_matches.tsv".format(
                masst_client.common_base_file_name(compound_name, out_filename_no_ext)
            )
            for compound_name in jobs_df["Compound"]
        ]
        jobs_df["finished"] = [Path(file).is_file() for file in jobs_df["file_path"]]
        jobs_df = jobs_df[~jobs_df["finished"]]
        logger.info(
            "Running fast microbe masst on input n={} spectra (total with already finished was {} spectra)".format(
                len(jobs_df), all_len
            )
        )
    else:
        logger.info(
            "Running fast microbe masst on input n={} spectra".format(len(jobs_df))
        )

    with ThreadPoolExecutor(parallel_queries) as executor:
        futures = [
            executor.submit(
                masst_client.query_usi_or_id,
                out_filename_no_ext,
                compound_id,
                name,
                precursor_mz_tol=precursor_mz_tol,
                mz_tol=mz_tol,
                min_cos=min_cos,
                min_matched_signals=min_matched_signals,
                analog=analog,
                analog_mass_below=analog_mass_below,
                analog_mass_above=analog_mass_above,
            )
            for compound_id, name in zip(jobs_df["input_id"], jobs_df["Compound"])
        ]

        wait(futures)
        jobs_df["success"] = [f.result() for f in futures]

    # return success rate
    total_jobs = len(jobs_df)
    return 1 if total_jobs == 0 else len(jobs_df[jobs_df["success"]]) / float(total_jobs)



def run_on_mgf(
        input_file,
        out_filename_no_ext,
        precursor_mz_tol=0.05,
        mz_tol=0.02,
        min_cos=0.7,
        min_matched_signals=3,
        analog=False,
        analog_mass_below=150,
        analog_mass_above=200,
        parallel_queries=100,
        skip_existing=False,
):
    ids, precursor_mzs, precursor_charges, lib_ids = [], [], [], []
    mzs, intensities = [], []

    with pyteomics.mgf.MGF(input_file) as f_in:
        for spectrum_dict in tqdm(f_in):
            abundances = spectrum_dict["intensity array"]
            if len(abundances) >= min_matched_signals:
                # GNPS library mgf has SPECTRUMID for IDs
                specid = spectrum_dict["params"].get("spectrumid", None)
                lib_ids.append(specid)
                specid = "_{}".format(specid) if specid else ""
                # scan number and optional specid
                ids.append(spectrum_dict["params"].get("scans", "") + specid)
                precursor_mzs.append(float(spectrum_dict["params"]["pepmass"][0]))
                if "charge" in spectrum_dict["params"]:
                    precursor_charges.append(int(spectrum_dict["params"]["charge"][0]))
                else:
                    precursor_charges.append(1)

                mzs.append(spectrum_dict["m/z array"])
                intensities.append(abundances)

    jobs_df = pd.DataFrame(
        {
            "Compound": ids,
            "lib_id": lib_ids,
            "precursor_mz": precursor_mzs,
            "precursor_charge": precursor_charges,
            "mzs": mzs,
            "intensities": intensities,
        }
    )

    if skip_existing:
        all_len = len(jobs_df)
        jobs_df["file_path"] = [
            "{}_matches.tsv".format(
                masst_client.common_base_file_name(compound_name, out_filename_no_ext)
            )
            for compound_name in ids
        ]
        jobs_df["finished"] = [Path(file).is_file() for file in jobs_df["file_path"]]
        jobs_df = jobs_df[~jobs_df["finished"]]
        logger.info(
            "Running fast microbe masst on input n={} spectra (total with already finished was {} spectra)".format(
                len(jobs_df), all_len
            )
        )
    else:
        logger.info(
            "Running fast microbe masst on input n={} spectra".format(len(jobs_df))
        )

    total_jobs = len(jobs_df)
    if total_jobs <= 1:
        jobs_df["success"] = [
            masst_client.query_spectrum(
                out_filename_no_ext,
                name,
                prec_mz,
                prec_charge,
                mz_array,
                intensity_array,
                precursor_mz_tol=precursor_mz_tol,
                mz_tol=mz_tol,
                min_cos=min_cos,
                min_matched_signals=min_matched_signals,
                analog=analog,
                analog_mass_below=analog_mass_below,
                analog_mass_above=analog_mass_above,
                lib_id=lib_id,
            )
            for name, lib_id, prec_mz, prec_charge, mz_array, intensity_array in zip(
                jobs_df["Compound"],
                jobs_df["lib_id"],
                jobs_df["precursor_mz"],
                jobs_df["precursor_charge"],
                jobs_df["mzs"],
                jobs_df["intensities"],
            )
        ]
    else:
        with ThreadPoolExecutor(parallel_queries) as executor:
            futures = [
                executor.submit(
                    masst_client.query_spectrum,
                    out_filename_no_ext,
                    name,
                    prec_mz,
                    prec_charge,
                    mz_array,
                    intensity_array,
                    precursor_mz_tol=precursor_mz_tol,
                    mz_tol=mz_tol,
                    min_cos=min_cos,
                    min_matched_signals=min_matched_signals,
                    analog=analog,
                    analog_mass_below=analog_mass_below,
                    analog_mass_above=analog_mass_above,
                    lib_id=lib_id,
                )
                for name, lib_id, prec_mz, prec_charge, mz_array, intensity_array in zip(
                    jobs_df["Compound"],
                    jobs_df["lib_id"],
                    jobs_df["precursor_mz"],
                    jobs_df["precursor_charge"],
                    jobs_df["mzs"],
                    jobs_df["intensities"],
                )
            ]

            wait(futures)
            jobs_df["success"] = [f.result() for f in futures]

    # return success rate
    total_jobs = len(jobs_df)
    return (
        1 if total_jobs == 0 else len(jobs_df[jobs_df["success"]]) / float(total_jobs)
    )


def create_params_label(
        analog,
        analog_mass_above,
        analog_mass_below,
        min_cos,
        min_matched_signals,
        mz_tol,
        precursor_mz_tol,
):
    params_label = (
        "min matched signals: {};  min cosine: {};  precursor m/z tolerance: {};  m/z "
        "tolerance: {};  analog: {}".format(
            min_matched_signals, min_cos, precursor_mz_tol, mz_tol, analog
        )
    )
    if analog:
        params_label += ";  analogs below m/z:{};  analogs above m/z:{}".format(
            analog_mass_below, analog_mass_above
        )
    return params_label


if __name__ == "__main__":
    # parsing the arguments (all optional)
    parser = argparse.ArgumentParser(description="Run fast microbeMASST in batch")
    parser.add_argument(
        "--in_file",
        type=str,
        help="input file either mgf with spectra or table that contains the two columns specified by "
             "usi_or_lib_id and compound_header",
        # default="../examples/example_links.tsv",
        default="../examples/small.mgf",
        # default="../examples/empty.mgf",
    )

    parser.add_argument(
        "--out_file",
        type=str,
        help="output html and other files, name without extension",
        default="../output/fastMASST",
    )

    # only for USI or lib ID file
    parser.add_argument(
        "--usi_or_lib_id",
        type=str,
        help="specify the usi or GNPS library id to search",
        default="USI",
    )
    parser.add_argument(
        "--compound_header",
        type=str,
        help="defines the header of the compound names, make sure to be path safe",
        default="Compound",
    )
    parser.add_argument(
        "--separator",
        type=str,
        help="separator for input file, e.g., \\t for tab",
        default="\t",
    )

    # MASST params
    parser.add_argument(
        "--precursor_mz_tol", type=float, help="precursor mz tolerance", default="0.05"
    )
    parser.add_argument(
        "--mz_tol", type=float, help="mz tolerance to match signals", default="0.05"
    )
    parser.add_argument(
        "--min_cos", type=float, help="Minimum cosine score for a match", default="0.7"
    )
    parser.add_argument(
        "--min_matched_signals", type=int, help="Minimum matched signals", default="4"
    )
    parser.add_argument(
        "--analog",
        type=bool,
        help="Search for analogs within mass window",
        default=False,
    )
    parser.add_argument(
        "--analog_mass_below",
        type=float,
        help="Maximum mass delta for analogs",
        default="150",
    )
    parser.add_argument(
        "--analog_mass_above",
        type=float,
        help="Maximum mass delta for analogs",
        default="200",
    )
    parser.add_argument(
        "--parallel_queries",
        type=int,
        help="the number of async queries. fastMASST step is IO bound so higher number than CPU "
        "speeds up the process",
        default="10",
    )
    parser.add_argument(
        "--skip_existing",
        type=bool,
        help="skip existing already processed entries",
        default=True,
    )

    args = parser.parse_args()

    try:
        success_rate = run_on_usi_list_or_mgf_file(
            in_file=args.in_file,
            out_file_no_extension=args.out_file,
            # only for USI list
            usi_or_lib_id=args.usi_or_lib_id,
            compound_name_header=args.compound_header,
            sep=args.separator,
            # matching parameters
            precursor_mz_tol=args.precursor_mz_tol,
            mz_tol=args.mz_tol,
            min_cos=args.min_cos,
            min_matched_signals=args.min_matched_signals,
            analog=args.analog,
            analog_mass_below=args.analog_mass_below,
            analog_mass_above=args.analog_mass_above,
            parallel_queries=args.parallel_queries,
            skip_existing=args.skip_existing,
        )
        logger.info(
            "Batch microbe MASST success rate (fastMASST query success) was %.3f",
            success_rate,
        )
        if success_rate < 1:
            sys.exit(1)
    except Exception as e:
        # exit with error
        logger.exception(e)
        sys.exit(1)

    # exit with OK
    sys.exit(0)
