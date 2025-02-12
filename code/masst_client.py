import sys
import logging
from tqdm import tqdm
import re
import argparse
from distutils.util import strtobool

import masst_utils
from masst_utils import DataBase
from utils import prepare_paths
from masst_tree import create_enriched_masst_tree
from masst_tree import create_combined_masst_tree
import masst_utils as masst
import usi_utils

MATCH_COLUMNS = ["Delta Mass", "USI", "Cosine", "Matching Peaks", "Status"]

LIB_COLUMNS = [
    "USI",
    "GNPSLibraryAccession",
    "Cosine",
    "Matching Peaks",
    "CompoundName",
    "Adduct",
    "Charge",
]

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# activate pandas tqdm progress_apply
tqdm.pandas()


def process_matches(
    file_name,
    compound_name,
    matches,
    library_matches,
    precursor_mz_tol,
    min_matched_signals,
    analog,
    input_label,
    params_label,
    usi=None,
):
    common_file = common_base_file_name(compound_name, file_name)

    # extract results
    extracted_results = masst.extract_matches_from_masst_results(
        matches,
        precursor_mz_tol,
        min_matched_signals,
        analog,
        limit_to_best_match_in_file=True,
        add_dataset_titles=False,
    )

    if analog:
        matches_df, analog_matches_df, original_matches_df = extracted_results
        # Save analog results separately
        analog_file = "{}_analog_matches.tsv".format(common_file)
        prepare_paths(file=analog_file)
        analog_matches_df[MATCH_COLUMNS + ["rounded_delta"]].to_csv(analog_file, index=False, sep="\t")

        # Here we want to export the unfiltered results
        masst_file = "{}_matches.tsv".format(common_file)
        prepare_paths(file=masst_file)
        matches_df[MATCH_COLUMNS].to_csv(masst_file, index=False, sep="\t")

        # Here we want to export all FasstMASST results
        original_masst_file = "{}_unfiltered_matches.tsv".format(common_file)
        prepare_paths(file=original_masst_file)
        original_matches_df[MATCH_COLUMNS].to_csv(original_masst_file, index=False, sep="\t")
    else:
        matches_df = extracted_results

        # always export match table even with 0 matches to mark that it was successful
        masst_file = "{}_matches.tsv".format(common_file)
        prepare_paths(file=masst_file)
        matches_df[MATCH_COLUMNS].to_csv(masst_file, index=False, sep="\t")

    lib_matches_df = masst.extract_matches_from_masst_results(
        library_matches, precursor_mz_tol, min_matched_signals, analog, False
    )
    if len(lib_matches_df) > 0:
        lib_matches_df[LIB_COLUMNS].to_csv(
            "{}_library.tsv".format(common_file), index=False, sep="\t"
        )

    if "grouped_by_dataset" not in matches:
        logger.debug("Missing datasets")
    # extract matches
    datasets_df = masst.extract_datasets_from_masst_results(matches, matches_df)
    if len(datasets_df) > 0:
        datasets_df.to_csv("{}_datasets.tsv".format(common_file), index=False, sep="\t")

    # add library matches to table
    lib_match_json = lib_matches_df.to_json(orient="records")

    # microbeMASST
    logger.debug("Exporting microbeMASST %s", compound_name)
    create_enriched_masst_tree(
        matches_df,
        masst.MICROBE_MASST,
        common_file=common_file,
        lib_match_json=lib_match_json,
        input_str=input_label,
        parameter_str=params_label,
        usi=usi,
        format_out_json=False,
        compress_out_html=True,
    )

    if analog:
        logger.debug("Exporting microbeMASST analog %s", compound_name)
        create_enriched_masst_tree(
            analog_matches_df,
            masst.MICROBE_MASST,
            common_file=common_file + "_analog",
            lib_match_json=lib_match_json,
            input_str=input_label,
            parameter_str=params_label,
            usi=usi,
            format_out_json=False,
            compress_out_html=True,
        )

    # plantMASST
    logger.debug("Exporting plantMASST %s", compound_name)
    create_enriched_masst_tree(
        matches_df,
        masst.PLANT_MASST,
        common_file=common_file,
        lib_match_json=lib_match_json,
        input_str=input_label,
        parameter_str=params_label,
        usi=usi,
        format_out_json=False,
        compress_out_html=True,
    )

    if analog:
        logger.debug("Exporting plantMASST analog %s", compound_name)
        create_enriched_masst_tree(
            analog_matches_df,
            masst.PLANT_MASST,
            common_file=common_file + "_analog",
            lib_match_json=lib_match_json,
            input_str=input_label,
            parameter_str=params_label,
            usi=usi,
            format_out_json=False,
            compress_out_html=True,
        )

    # tissueMASST
    logger.debug("Exporting tissueMASST %s", compound_name)
    create_enriched_masst_tree(
        matches_df,
        masst.TISSUE_MASST,
        common_file=common_file,
        lib_match_json=lib_match_json,
        input_str=input_label,
        parameter_str=params_label,
        usi=usi,
        format_out_json=False,
        compress_out_html=True,
    )

    if analog:
        logger.debug("Exporting tissueMASST analog %s", compound_name)
        create_enriched_masst_tree(
            analog_matches_df,
            masst.TISSUE_MASST,
            common_file=common_file + "_analog",
            lib_match_json=lib_match_json,
            input_str=input_label,
            parameter_str=params_label,
            usi=usi,
            format_out_json=False,
            compress_out_html=True,
        )

    # foodMASST
    logger.debug("Exporting foodMASST %s", compound_name)
    create_enriched_masst_tree(
        matches_df,
        masst.FOOD_MASST,
        common_file=common_file,
        lib_match_json=lib_match_json,
        input_str=input_label,
        parameter_str=params_label,
        usi=usi,
        format_out_json=False,
        compress_out_html=True,
    )

    if analog:
        logger.debug("Exporting foodMASST analog %s", compound_name)
        create_enriched_masst_tree(
            analog_matches_df,
            masst.FOOD_MASST,
            common_file=common_file + "_analog",
            lib_match_json=lib_match_json,
            input_str=input_label,
            parameter_str=params_label,
            usi=usi,
            format_out_json=False,
            compress_out_html=True,
        )

    # personalCareProductMASST
    logger.debug("Exporting personalCareProductMASST %s", compound_name)
    create_enriched_masst_tree(
        matches_df,
        masst.PERSONALCAREPRODUCT_MASST,
        common_file=common_file,
        lib_match_json=lib_match_json,
        input_str=input_label,
        parameter_str=params_label,
        usi=usi,
        format_out_json=False,
        compress_out_html=True,
    )

    if analog:
        logger.debug("Exporting personalCareProductMASST analog %s", compound_name)
        create_enriched_masst_tree(
            analog_matches_df,
            masst.PERSONALCAREPRODUCT_MASST,
            common_file=common_file + "_analog",
            lib_match_json=lib_match_json,
            input_str=input_label,
            parameter_str=params_label,
            usi=usi,
            format_out_json=False,
            compress_out_html=True,
        )

    # microbiomeMASST
    logger.debug("Exporting microbiomeMASST %s", compound_name)
    create_enriched_masst_tree(
        matches_df,
        masst.MICROBIOME_MASST,
        common_file=common_file,
        lib_match_json=lib_match_json,
        input_str=input_label,
        parameter_str=params_label,
        usi=usi,
        format_out_json=False,
        compress_out_html=True,
    )

    if analog:
        logger.debug("Exporting microbiomeMASST analog %s", compound_name)
        create_enriched_masst_tree(
            analog_matches_df,
            masst.MICROBIOME_MASST,
            common_file=common_file + "_analog",
            lib_match_json=lib_match_json,
            input_str=input_label,
            parameter_str=params_label,
            usi=usi,
            format_out_json=False,
            compress_out_html=True,
        )

    # combined from all
    logger.debug("Exporting combined tree %s", compound_name)
    create_combined_masst_tree(
        matches_df,
        common_file=common_file,
        lib_match_json=lib_match_json,
        input_str=input_label,
        parameter_str=params_label,
        usi=usi,
        format_out_json=False,
        compress_out_html=True,
    )

    if analog:
        logger.debug("Exporting combined tree analog %s", compound_name)
        create_combined_masst_tree(
            analog_matches_df,
            common_file=common_file + "_analog",
            lib_match_json=lib_match_json,
            input_str=input_label,
            parameter_str=params_label,
            usi=usi,
            format_out_json=False,
            compress_out_html=True,
        )

    return matches_df


def common_base_file_name(compound_name, file_name):
    if compound_name:
        return "{}_{}".format(file_name, compound_name.replace(" ", "_"))
    else:
        return "{}".format(file_name)


def query_usi_or_id(
    file_name,
    usi_or_lib_id,
    compound_name,
    precursor_mz_tol=0.05,
    mz_tol=0.02,
    min_cos=0.7,
    min_matched_signals=3,
    analog=False,
    analog_mass_below=130,
    analog_mass_above=200,
    database: str | DataBase = None,
    library: str | DataBase = None,
):
    """
    NOTE: database and library are the fasst database, if None, we fall back on defaults provided by the system, otherwise we can set a string

    :return: True if fastmasst query was successful otherwise False
    """
    # might raise exception for service
    try:
        logger.debug("Query fastMASST id:%s  of %s", usi_or_lib_id, compound_name)

        if database is None:
            database = masst.DataBase.metabolomicspanrepo_index_latest

        matches = masst.fast_masst(
            usi_or_lib_id,
            precursor_mz_tol=precursor_mz_tol,
            mz_tol=mz_tol,
            min_cos=min_cos,
            analog=analog,
            analog_mass_below=analog_mass_below,
            analog_mass_above=analog_mass_above,
            database=database,
        )

        if not matches or "results" not in matches:
            logger.debug(
                "Empty fastMASST response for compound %s with id %s",
                compound_name,
                usi_or_lib_id,
            )
            return False

        if len(matches["results"]) == 0:
            export_empty_masst_results(compound_name, file_name)
            # succeeded with 0 matches
            # currently fastMASST returns empty response without results dictionary
            return True

        if library is None:
            library = masst.DataBase.gnpslibrary

        library_matches = masst.fast_masst(
            usi_or_lib_id,
            precursor_mz_tol=precursor_mz_tol,
            mz_tol=mz_tol,
            min_cos=min_cos,
            analog=False,
            database=library,
        )

        params_label = create_params_label(
            analog,
            analog_mass_above,
            analog_mass_below,
            min_cos,
            min_matched_signals,
            mz_tol,
            precursor_mz_tol,
        )
        input_label = "ID: {};  Descriptor: {}".format(usi_or_lib_id, compound_name)

        process_matches(
            file_name,
            compound_name,
            matches,
            library_matches,
            precursor_mz_tol,
            min_matched_signals,
            analog,
            input_label,
            params_label,
            usi_utils.ensure_usi(usi_or_lib_id),
        )
        return True
    except Exception as e:
        # logger.exception(e)
        return False


def query_spectrum(
    file_name,
    compound_name,
    precursor_mz,
    precursor_charge,
    mzs,
    intensities,
    precursor_mz_tol=0.05,
    mz_tol=0.02,
    min_cos=0.7,
    min_matched_signals=3,
    analog=False,
    analog_mass_below=130,
    analog_mass_above=200,
    lib_id=None,
    database: str | DataBase = None,
    library: str | DataBase = None,
):
    """
    NOTE: database and library are the fasst database, if None, we fall back on defaults provided by the system, otherwise we can set a string

    :return: True if fast masst query was successful otherwise False
    """
    # might raise exception for service
    try:
        if database is None:
            database = masst.DataBase.metabolomicspanrepo_index_latest

        matches, filtered_dps = masst.fast_masst_spectrum(
            mzs=mzs,
            intensities=intensities,
            precursor_mz=precursor_mz,
            precursor_charge=precursor_charge,
            precursor_mz_tol=precursor_mz_tol,
            mz_tol=mz_tol,
            min_cos=min_cos,
            analog=analog,
            analog_mass_below=analog_mass_below,
            analog_mass_above=analog_mass_above,
            database=database,
        )
        if not matches or "results" not in matches:
            # export empty masst results file to signal that service was successful
            logger.debug("Empty fastMASST response for spectrum %s", compound_name)
            return False

        if len(matches["results"]) == 0:
            export_empty_masst_results(compound_name, file_name)
            return True

        if library is None:
            library = masst.DataBase.gnpslibrary

        library_matches, _ = masst.fast_masst_spectrum(
            mzs=mzs,
            intensities=intensities,
            precursor_mz=precursor_mz,
            precursor_charge=precursor_charge,
            precursor_mz_tol=precursor_mz_tol,
            mz_tol=mz_tol,
            min_cos=min_cos,
            analog=False,
            database=library,
        )

        params_label = create_params_label(
            analog,
            analog_mass_above,
            analog_mass_below,
            min_cos,
            min_matched_signals,
            mz_tol,
            precursor_mz_tol,
        )
        input_label = "Descriptor: {};  Precursor m/z: {};  Data points:{}".format(
            compound_name, round(precursor_mz, 5), len(filtered_dps)
        )

        usi = usi_utils.ensure_usi(lib_id)

        process_matches(
            file_name,
            compound_name,
            matches,
            library_matches,
            precursor_mz_tol,
            min_matched_signals,
            analog,
            input_label,
            params_label,
            usi,
        )
        return True
    except Exception as e:
        return False


def export_empty_masst_results(compound_name, file_name):
    try:
        path = "{}_matches.tsv".format(common_base_file_name(compound_name, file_name))
        with open(path, "w") as file:
            file.write("USI	Cosine	Matching Peaks	Status\n")
    except:
        pass


def path_safe(file):
    return re.sub("[^-a-zA-Z0-9_.() ]+", "_", file)


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
        "--usi_or_lib_id",
        type=str,
        help="usi or GNPS library ID",
        default='mzspec:MSV000095331:peak/mzML_3/Full_syncom_T72_drug_mix_2_2.mzML:scan:2176',  # microbiome pos control
        # default="mzspec:MSV000095003:P2_E1_584:scan:2609",  # personal care pos control
        # default="CCMSLIB00005435899",  # moroidin
        # default="CCMSLIB00005883945",  # tryptophan 6
        # default="CCMSLIB00004679239", # commendamide
        # default="CCMSLIB00006582001",
    )
    parser.add_argument(
        "--out_file",
        type=str,
        help="output html and other files, name without extension",
        default="../output/fastMASST",
    )

    # only for USI or lib ID file
    parser.add_argument("--compound_name", type=str, help="compound name", default="")

    # search database
    parser.add_argument(
        "--database",
        type=str,
        help="fasst database for public data",
        default="metabolomicspanrepo_index_latest",
    )
    parser.add_argument(
        "--library", type=str, help="fasst library for reference spectra", default=None
    )

    # MASST params
    parser.add_argument(
        "--precursor_mz_tol", type=float, help="precursor mz tolerance", default="0.05"
    )
    parser.add_argument(
        "--mz_tol", type=float, help="mz tolerance to match signals", default="0.02"
    )
    parser.add_argument(
        "--min_cos", type=float, help="Minimum cosine score for a match", default="0.7"
    )
    parser.add_argument(
        "--min_matched_signals", type=int, help="Minimum matched signals", default="3"
    )
    parser.add_argument(
        "--analog",
        type=lambda x: bool(strtobool(str(x.strip()))),
        help="Search for analogs within mass window",
        default=False,
    )
    parser.add_argument(
        "--analog_mass_below",
        type=int,
        help="Maximum mass delta for analogs",
        default="130",
    )
    parser.add_argument(
        "--analog_mass_above",
        type=int,
        help="Maximum mass delta for analogs",
        default="200",
    )

    args = parser.parse_args()

    try:
        usi_or_lib_id = args.usi_or_lib_id
        compound_name = args.compound_name
        out_file = args.out_file

        # MASST parameters
        precursor_mz_tol = args.precursor_mz_tol
        mz_tol = args.mz_tol
        min_cos = args.min_cos
        min_matched_signals = args.min_matched_signals
        analog = args.analog
        analog_mass_below = args.analog_mass_below
        analog_mass_above = args.analog_mass_above
        database = args.database
        library = args.library

        query_usi_or_id(
            out_file,
            usi_or_lib_id,
            compound_name,
            precursor_mz_tol=precursor_mz_tol,
            mz_tol=mz_tol,
            min_cos=min_cos,
            min_matched_signals=min_matched_signals,
            analog=analog,
            analog_mass_below=analog_mass_below,
            analog_mass_above=analog_mass_above,
            database=database,
            library=library,
        )
    except Exception as e:
        # exit with error
        logger.exception(e)
        sys.exit(1)

    # exit with OK
    sys.exit(0)
