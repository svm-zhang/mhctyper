import argparse

from tinyscibio import parse_path


def parse_cmd() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam",
        metavar="FILE",
        type=parse_path,
        required=True,
        help="specify path to BAM file",
    )
    parser.add_argument(
        "--freq",
        metavar="FILE",
        type=parse_path,
        required=True,
        help="specify path to HLA frequency file",
    )
    parser.add_argument(
        "--race",
        metavar="STR",
        type=str,
        default="Unknown",
        help="specify race [default: Unknown]",
    )
    parser.add_argument(
        "--out",
        metavar="FILE",
        type=parse_path,
        required=True,
        help="specify path to output hlatyping result file",
    )
    parser.add_argument(
        "--min_ecnt",
        metavar="INT",
        type=int,
        default=999,
        help="specify minimum # of mm events (999)",
    )
    parser.add_argument(
        "--nproc",
        metavar="INT",
        type=int,
        default=8,
        help="specify # processes to use (8)",
    )
    return parser
