import logging
import sys

import masst_batch_client

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

files = [
    # ("D:\Robin\git\microbe_masst\local_files\ipsita\Dihydroxy_MSCluster.mgf", "../output/fastMASSTDiOHBA"),
    (r"../examples/example_links.tsv", "../output/examples/fastMASST_"),
    (r"../examples/example_links2.tsv", "../output/examples/fastMASST_"),
    (r"../examples/example_links3.tsv", "../output/examples/fastMASST_"),
]

if __name__ == "__main__":
    for file, out_file in files:
        try:
            logger.info("Starting new job for input: {}".format(file))
            sep = "," if file.endswith("csv") else "\t" # only used if tabular format not for mgf
            masst_batch_client.run_on_usi_list_or_mgf_file(
                in_file=file,
                out_file_no_extension=out_file,
                min_matched_signals=4,
                parallel_queries=22,
                skip_existing=True,
                sep=sep
            )
        except:
            pass
    # exit with OK
    sys.exit(0)
