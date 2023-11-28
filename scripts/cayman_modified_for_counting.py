import logging
import pathlib
import os
import sys
import argparse

# pylint: disable=W0611
from gqlib.db.db_import import SmallDatabaseImporter
from gqlib.profilers import RegionQuantifier
from gqlib.runners.alignment_runner import BwaMemRunner
from gqlib.ui.validation import check_bwa_index, check_input_reads

logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(
        description="Quantify gene counts based on modified Cayman script",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        '--bam_file', '-b', type=str, required=True,
        help='BAM file to count alignments from.'
    )

    parser.add_argument(
        '--bed_file', '-d', type=str, required=True,
        help='GFF file with genes to be counted.'
    )

    args = parser.parse_args()
    bam_file = args.bam_file
    gff_file = args.bed_file
    min_identity = 0.97
    min_seqlen = 50


    # profiler.count_alignments(
    #             bam_file=bam_file,
    #             aln_format="bam",
    #             min_identity=0.97,# default: 0.97
    #             min_seqlen=50, # default: 45 to keep to same as KE
    #         )

    # db_importer = SmallDatabaseImporter(
    #     logger, args.bed_file, single_category="cazy", db_format=args.db_format,
    # )
    # logger.info("Finished loading database.")

    # profiler = RegionQuantifier(
    #     db=db_importer,
    #     out_prefix=args.out_prefix,
    #     ambig_mode="1overN",
    #     reference_type="domain",
    # )

    # aln_runner = BwaMemRunner(
    #     args.cpus_for_alignment,
    #     args.bwa_index,
    #     sample_id=os.path.basename(args.out_prefix),
    # )

    # for input_type, *reads in input_data:

    #     logger.info("Running %s alignment: %s", input_type, ",".join(reads))
    #     proc, call = aln_runner.run(
    #         reads,
    #         single_end_reads=input_type == "single",
    #     )

    #     try:
            

    #     except Exception as err:
    #         if isinstance(err, ValueError) and str(err).strip() == "file does not contain alignment data":
    #             # pylint: disable=W1203
    #             logger.error("Failed to align. This could have different reasons:")
    #             logger.error(f"* Is `{args.aligner}` installed and on the path? Type `bwa mem` and see what happens.")
    #             logger.error("* Syntax errors or missing files. Please try running the aligner call below manually to troubleshoot the problem.")
    #             logger.error("* Alignment stream was interrupted, perhaps due to a memory issue.")
                
    #             logger.error("Aligner call was:")
    #             logger.error("%s", call)
    #             sys.exit(1)

    #         logger.error("Encountered problems digesting the alignment stream:")
    #         logger.error("%s", err)
    #         logger.error("Aligner call was:")
    #         logger.error("%s", call)
    #         logger.error("Shutting down.")
    #         sys.exit(1)

    # profiler.finalise(restrict_reports=("raw", "rpkm",))


if __name__ == "__main__":
    main()
