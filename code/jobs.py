import logging
import masst_batch_client

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

files = [
    # ("D:\Robin\git\microbe_masst\local_files\ipsita\Dihydroxy_MSCluster.mgf", "../output/fastMASSTDiOHBA"),
    # ("D:\Robin\git\microbe_masst\local_files\ipsita\Monohydroxy_MSCluster.mgf", "../output/fastMASSTMonoOHBA"),
    # ("D:\Robin\git\microbe_masst\local_files\ipsita\Tetrahydroxy_MSCluster.mgf", "../output/fastMASSTTetraOHBA"),
    # ("D:\Robin\git\microbe_masst\local_files\ipsita\Trihydroxy_MSCluster.mgf", "../output/fastMASSTTriOHBA"),
    # BASF
    (r"D:\Robin\git\microbe_masst\local_files\BASF\RP18_neg_microbial_network.mgf",
     "../output/BASF/fastMASST_RP18_neg_"),
    (r"D:\Robin\git\microbe_masst\local_files\BASF\RP18_neg_microbial_network.mgf",
     "../output/BASF/fastMASST_RP18_neg_"),
    # (r"D:\Robin\git\microbe_masst\local_files\BASF\RP18_pos_microbial_network.mgf", "../output/BASF/fastMASST_RP18_pos_"),
    # (r"D:\Robin\git\microbe_masst\local_files\BASF\RP18_pos_microbial_network.mgf", "../output/BASF/fastMASST_RP18_pos_"),
    # (r"D:\Robin\git\microbe_masst\local_files\BASF\HILIC_pos_microbial_network.mgf", "../output/BASF/fastMASST_HILIC_pos_"),
    # (r"D:\Robin\git\microbe_masst\local_files\BASF\HILIC_neg_microbial_network.mgf", "../output/BASF/fastMASST_HILIC_neg_"),
    #
    # (r"D:\Robin\git\microbe_masst\local_files\BEAM\BEAM_FECAL.mgf", "../output/beam/fecal/fastMASST_beam_fecal"),
    # (r"D:\Robin\git\microbe_masst\local_files\BEAM\BEAM_SERUM.mgf", "../output/beam/serum/fastMASST_beam_serum"),
    # (r"D:\Robin\git\microbe_masst\local_files\BEAM\BEAM_FECAL.mgf", "../output/beam/fecal/fastMASST_beam_fecal"),
    # (r"D:\Robin\git\microbe_masst\local_files\BEAM\BEAM_SERUM.mgf", "../output/beam/serum/fastMASST_beam_serum"),
    # (r"D:\Robin\git\microbe_masst\local_files\BEAM\BEAM_FECAL.mgf", "../output/beam/fecal/fastMASST_beam_fecal"),
    # (r"D:\Robin\git\microbe_masst\local_files\BEAM\BEAM_SERUM.mgf", "../output/beam/serum/fastMASST_beam_serum"),
    # (r"D:\Robin\git\microbe_masst\local_files\BEAM\BEAM_FECAL.mgf", "../output/beam/fecal/fastMASST_beam_fecal"),
    # (r"D:\Robin\git\microbe_masst\local_files\BEAM\BEAM_SERUM.mgf", "../output/beam/serum/fastMASST_beam_serum"),
    # (r"D:\Robin\git\microbe_masst\local_files\BEAM\BEAM_FECAL.mgf", "../output/beam/fecal/fastMASST_beam_fecal"),
    # (r"D:\Robin\git\microbe_masst\local_files\BEAM\BEAM_SERUM.mgf", "../output/beam/serum/fastMASST_beam_serum"),
    # (r"D:\Robin\git\microbe_masst\local_files\BEAM\BEAM_FECAL.mgf", "../output/beam/fecal/fastMASST_beam_fecal"),
    # (r"D:\Robin\git\microbe_masst\local_files\BEAM\BEAM_SERUM.mgf", "../output/beam/serum/fastMASST_beam_serum"),
    # (r"D:\Robin\git\microbe_masst\local_files\BEAM\BEAM_FECAL.mgf", "../output/beam/fecal/fastMASST_beam_fecal"),
    # (r"D:\Robin\git\microbe_masst\local_files\BEAM\BEAM_SERUM.mgf", "../output/beam/serum/fastMASST_beam_serum"),
    (r"D:\Robin\git\microbe_masst\local_files\ALL_GNPS_NO_PROPOGATED.mgf",
     "../output/library/fastMASST_library_"),
    # nina space station
    # (r"D:\Robin\git\microbe_masst\local_files\nina_space\20220927_3DMM_rerun_execlusion_one.mgf",
    #  "../output/nina_space/fastMASST_nina_in_space_"),
    # (r"D:\Robin\git\microbe_masst\local_files\nina_space\20220927_3DMM_rerun_execlusion_one.mgf",
    #  "../output/nina_space/fastMASST_nina_in_space_"),
    # (r"D:\Robin\git\microbe_masst\local_files\nina_space\20220927_3DMM_rerun_execlusion_one.mgf",
    #  "../output/nina_space/fastMASST_nina_in_space_"),
    # (r"D:\Robin\git\microbe_masst\local_files\nina_space\20220927_3DMM_rerun_execlusion_one.mgf",
    #  "../output/nina_space/fastMASST_nina_in_space_"),

    # bile second
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Dihydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_2OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Monohydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_1OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Non_hydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_0OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Pentahydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_5OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\tetrahydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_4OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Trihydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_3OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Dihydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_2OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Monohydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_1OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Non_hydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_0OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Pentahydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_5OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\tetrahydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_4OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Trihydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_3OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Dihydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_2OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Monohydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_1OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Non_hydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_0OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Pentahydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_5OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\tetrahydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_4OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Trihydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_3OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Dihydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_2OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Monohydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_1OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Non_hydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_0OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Pentahydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_5OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\tetrahydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_4OH_"),
    # (r"D:\Robin\git\microbe_masst\local_files\ipsita\bile_non_refined\Trihydroxy_non_refined.mgf",
    #  "../output/bile_unrefined/fastMASST_bile_unref_3OH_"),
    #
    # # other
    # (r"D:\Robin\git\microbe_masst\local_files\221005_gwas_rats_fbmn.mgf", "../output/gwas/fastMASST_gwas"),
    # (r"D:\Robin\git\microbe_masst\local_files\casmi_pos_sirius\bifido.mgf", "../output/bifido/fastMASST_"),
    # ("D:\Robin\git\microbe_masst\local_files\casmi_pos_sirius\MIND.mgf", "../output/MIND/fastMASST_MIND"),
]

if __name__ == "__main__":

    for file, out_file in files:
        try:
            logger.info("Starting new job for input: {}".format(file))
            masst_batch_client.run_on_mgf(
                input_file=file,
                out_filename_no_ext=out_file,
                min_matched_signals=4,
                parallel_queries=20,
                skip_existing=True
            )
        except:
            pass
