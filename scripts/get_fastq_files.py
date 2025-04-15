import os
import sys
import glob
import logging
import shutil

import argparse
import pandas as pd

logger = logging.getLogger("get_fastq")


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--build_dir", "-b", help="Project directory", required=True)
    parser.add_argument(
        "--fastq_dir", "-f", help="Location of fastq files", required=True
    )
    parser.add_argument(
        "--seq_type",
        help="Machine that generated fastq files. (Affects file names)",
        choices=["MiSeq", "HiSeq", "NovaSeq", "DRAGEN"],
        required=True,
    )
    parser.add_argument(
        "--verbose",
        "-v",
        help="Whether to print a bunch of output",
        action="store_true",
        default=False,
    )

    return parser


def format_wells(row, well_field="pcr_well"):
    well = row[well_field]
    col = well[0]
    row = well[1:]
    return "{}{:02d}".format(col, int(row))


def make_file_names(row, seq_type):
    if seq_type == "HiSeq":
        return "[{fc_lane}]_{fc_name}.[{fc_lane}].{index_1}_{index_2}.unmapped.{read_type}".format(  # HiSeq
            fc_name=row["flowcell_name"],
            fc_lane=row["flowcell_lane"],
            index_1=row["IndexBarcode1"],
            index_2=row["IndexBarcode2"],
            read_type="*",
        )

    elif seq_type == "NovaSeq":
        return "[{fc_lane}]_{fc_name}.[{fc_lane}].{index_1}_{index_2_rc}.unmapped.{read_type}".format(  # NovaSeq
            fc_name=row["flowcell_name"],
            fc_lane=row["flowcell_lane"],
            index_1=row["IndexBarcode1"],
            index_2_rc=reverse_complement(row["IndexBarcode2"]),
            read_type="*",
        )

    elif seq_type == "DRAGEN":
        return "{fc_name}_[{fc_lane}]_{PoolTubeBarcode}_{index_1}-{index_2_rc}_{SampleNum}_{fc_L_lane}_{read_type}".format(  # NextSeq DRAGEN
            fc_name=row["flowcell_name"],
            fc_lane=row["flowcell_lane"],
            index_1=row["IndexBarcode1"],
            index_2_rc=row["IndexBarcode2"],
            read_type="*",
            PoolTubeBarcode="*",
            SampleNum="*",
            fc_L_lane="*",
        )

    elif seq_type == "MiSeq":
        return "{pcr_plate}_{pcr_well}_*.fastq.gz".format(  # MiSeq
            pcr_plate=row["pcr_plate"], pcr_well=row["pcr_well"]
        )
    else:
        raise ValueError("unknown sequencer type")


def reverse_complement(sequence):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join([complement[base] for base in reversed(sequence)])


def main(args):
    project_info = pd.read_csv(os.path.join(args.build_dir, "sample_meta.csv"))
    project_info["pcr_well"] = project_info.apply(format_wells, axis=1)

    proj_files = project_info.apply(
        lambda row: make_file_names(row, args.seq_type), axis=1
    )

    paths = proj_files.apply(lambda f: glob.glob(os.path.join(args.fastq_dir, f)))
    paths = [file for row in paths for file in row]

    proj_dir = args.build_dir
    fastq_out = os.path.join(proj_dir, "fastq/")

    logger.info("dir: {}".format(fastq_out))
    if os.path.isdir(fastq_out):
        logger.info("directory exists")
    else:
        os.umask(0o002) # Change masking permission
        os.makedirs(fastq_out)

    logger.debug(fastq_out)
    logger.debug(paths)
    for f in paths:
        # fl.append(glob.glob(f))
        dst = shutil.copy(f, fastq_out)
        if dst:  # only works on python3
            logger.info(f + " --> " + dst)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=level)
    logger.info("args:  {}".format(args))

    main(args)
