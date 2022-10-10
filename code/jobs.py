import logging
import masst_batch_client

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


if __name__ == "__main__":
    files = [
        # ("D:\Robin\git\microbe_masst\local_files\ipsita\Dihydroxy_MSCluster.mgf", "../output/fastMASSTDiOHBA"),
        # ("D:\Robin\git\microbe_masst\local_files\ipsita\Monohydroxy_MSCluster.mgf", "../output/fastMASSTMonoOHBA"),
        # ("D:\Robin\git\microbe_masst\local_files\ipsita\Tetrahydroxy_MSCluster.mgf", "../output/fastMASSTTetraOHBA"),
        # ("D:\Robin\git\microbe_masst\local_files\ipsita\Trihydroxy_MSCluster.mgf", "../output/fastMASSTTriOHBA"),
        (r"D:\Robin\git\microbe_masst\local_files\casmi_pos_sirius\bifido.mgf", "../output/bifido/fastMASST_"),
        ("D:\Robin\git\microbe_masst\local_files\casmi_pos_sirius\MIND.mgf", "../output/MIND/fastMASST_MIND"),
        (r"D:\Robin\git\microbe_masst\local_files\221005_gwas_rats_fbmn.mgf", "../output/gwas/fastMASST_gwas"),
    ]

    for file, out_file in files:
        try:
            logger.info("Starting new job for input: {}".format(file))
            masst_batch_client.run_on_mgf(
                input_file=file,
                out_filename_no_ext=out_file,
                min_matched_signals=4,
                parallel_queries=2,
                skip_existing=True
            )
        except:
            pass
