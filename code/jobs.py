import logging
import sys

import masst_batch_client

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

files = [
    # (r"D:\git\microbe_masst\local_files\piper_vs_microbiome_iimn_gnps.mgf", "../output/piper_vs_microbiome/fastMASST_"),
    # (r"D:\Data\salt_infusion\iimn\salt_infusion_iimn_gnps.mgf", "../output/salt_infusion/fastMASST_"),
    # (r"../examples/plant_mass_examples.tsv", "../output/plant/MASST_"),
    # (r"../examples/wender_plants.csv", "../output/wender/MASST_"),
    # (r"../examples/example_links.tsv", "../output/examples/fastMASST_"),
    # (r"../examples/example_links2.tsv", "../output/examples/fastMASST_"),
    (r"../examples/example_links3.tsv", "../output/examples/fastMASST_"),
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
                min_matched_signals=4,
                parallel_queries=5,
                skip_existing=True,
                sep=sep,
            )
        except:
            pass
    # exit with OK
    sys.exit(0)
