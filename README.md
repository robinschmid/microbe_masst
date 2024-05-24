[![DOI](https://zenodo.org/badge/492844724.svg)](https://zenodo.org/badge/latestdoi/492844724)

# How to cite?
_microbeMASST: a taxonomically informed mass spectrometry search tool for microbial metabolomics data_
in [Nature Microbiology](https://www.nature.com/articles/s41564-023-01575-9)

# fast microbeMASST
Using MASST or fastMASST to search all public datasets on GNPS/MassIVE, adding extensive metadata onto an ontology for microbes.

## Batch jobs runner
1. Navigate to the jobs.py file add entries to the files list as `("input file", "output_directory/file_prefix)`
2. Run jobs.py
3. Some entries might fail - adding the same input output file combination multiple times to run iteratively with the option: `skip_existing=True`

## Usage

```
# single USI or GNPS library ID
microbe_masst.py
# batch jobs
microbe_masst_batch.py
```
