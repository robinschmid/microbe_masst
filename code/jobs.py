import logging
import sys

import masst_batch_client

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

files = [
    (r"../examples/plant_mass_examples.tsv", "../output/plant/MASST_"),
    (r"../examples/example_links5.csv", "../output/global/example_links6"),
    # (r"../examples/vta_filter.mgf", "../output/examples/vta_"),
]

if __name__ == "__main__":
    for file, out_file in files:
        try:
            logger.info("Starting new job for input: {}".format(file))
            sep = (
                "," if file.endswith("csv") else "\t"
            )  # only used if tabular format not for mgf
            masst_batch_client.run_on_usi_list_or_mgf_file(
                in_file=file,
                out_file_no_extension=out_file,
                min_cos=0.7,
                mz_tol=0.02,
                precursor_mz_tol=0.02,
                min_matched_signals=4,
                parallel_queries=6,
                skip_existing=True,
                analog=False,
                sep=sep,
            )
        except:
            pass
    # exit with OK
    sys.exit(0)
