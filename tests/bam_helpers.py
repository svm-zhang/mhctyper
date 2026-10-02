"""Tiny BAM fixtures; sequences and qualities use stored alignment orientation."""

from pathlib import Path

import pysam


def make_pair(qname, reference_id=0, *, secondary=False):
    """Make a four-base, all-match FR pair; callers can alter individual records."""
    records = []
    for flag, start, mate_start, span in [(99, 10, 30, 24), (147, 30, 10, -24)]:
        read = pysam.AlignedSegment()
        read.query_name = qname
        read.query_sequence = "AAAA"
        read.query_qualities = [30, 30, 30, 30]
        read.flag = flag | (256 if secondary else 0)
        read.reference_id = reference_id
        read.reference_start = start
        read.mapping_quality = 0 if secondary else 60
        read.cigarstring = "4M"
        read.next_reference_id = reference_id
        read.next_reference_start = mate_start
        read.template_length = span
        read.set_tag("MD", "4")
        read.set_tag("NM", 0)
        read.set_tag("RG", "lane_1")
        records.append(read)
    return records


def write_bam(path: Path, references, records):
    """Write records in coordinate order and index without external executables."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": 100} for name in references],
        "RG": [{"ID": "lane_1", "SM": "sample_1"}],
    }
    ordered_records = sorted(
        records, key=lambda read: (read.reference_id, read.reference_start)
    )
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for record in ordered_records:
            bam.write(record)
    pysam.index(str(path))
    return path
