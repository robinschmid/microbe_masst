import pytest

import microbe_masst


def test_microbe_masst():
    result = microbe_masst.run_microbe_masst(
        "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005883671"
    )
    assert result is not None
