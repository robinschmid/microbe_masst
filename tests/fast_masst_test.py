import pytest

import masst_utils
import usi_utils


def test_fast_masst_spectrum():
    mzs, intensities = zip(
        *[
            (80.9734, 955969.8),
            (81.9816, 542119.2),
            (98.9841, 483893630.0),
            (116.9947, 1605324.2),
            (127.0155, 182958080.0),
            (131.0102, 878951.4),
            (155.0467, 73527150.0),
            (183.0781, 16294011.0),
        ]
    )

    matches, _ = masst_utils.fast_masst_spectrum(
        mzs, intensities, precursor_mz=183.078, precursor_charge=1
    )
    assert len(matches) > 1


def test_fast_masst_usi():
    usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005883671"
    matches = masst_utils.fast_masst(usi)
    assert len(matches) > 1


def test_fast_masst_moroidin():
    usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005435899"
    matches = masst_utils.fast_masst(
        usi,
        database=masst_utils.DataBase.gnpsdata_index_11_25_23,
        mz_tol=0.05,
    )
    assert len(matches) > 1
    # this currently fails when mz_tol is 0.05 - fixed by setting 0.02
    assert len(matches["results"]) > 1


def test_fastmasst_trp():
    usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005883950"
    matches = masst_utils.fast_masst(
        usi,
        precursor_mz_tol=0.05,
        mz_tol=0.05,
    )
    df = masst_utils.extract_matches_from_masst_results(
        matches,
        precursor_mz_tol=0.5,
        min_matched_signals=4,
        analog=False,
        limit_to_best_match_in_file=True,
    )
    assert len(df) > 1
    assert len(df.index) > 50, 000


def test_fastmasst_trp_analog():
    usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005883950"
    matches = masst_utils.fast_masst(usi, precursor_mz_tol=0.05, analog=True)
    df = masst_utils.extract_matches_from_masst_results(
        matches,
        precursor_mz_tol=0.5,
        min_matched_signals=4,
        analog=True,
        limit_to_best_match_in_file=True,
    )
    assert len(df) > 1
    assert len(df.index) > 444, 000


def test_get_spectrum():
    spec = usi_utils.get_spectrum(
        "mzspec:MSV000090162:peak/F2202_NRP_purification_check_mzXML/Sample_10_RB5_01_63425"
        ".mzXML:scan:583"
    )
    assert spec["precursor_mz"] > 0
