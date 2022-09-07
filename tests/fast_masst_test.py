import pytest

import masst_utils


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

    matches,_ = masst_utils.fast_masst_spectrum(
        mzs, intensities, precursor_mz=183.078, precursor_charge=1
    )
    assert len(matches) > 1

def test_fast_masst_usi():
    usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005883671"
    matches = masst_utils.fast_masst(usi)
    assert len(matches) > 1
