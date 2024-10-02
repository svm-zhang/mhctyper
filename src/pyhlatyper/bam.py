import itertools
import re
from collections.abc import MutableSequence
from dataclasses import dataclass, field
from typing import Any, Generator

import pysam
from libscibio import parse_path


@dataclass
class BAMetadata:
    fspath: str
    sort_by: str = field(init=False, default="")
    references: list[dict[str, int]] = field(init=False, default_factory=list)
    read_groups: list[dict[str, str]] = field(init=False, default_factory=list)

    def __post_init__(self) -> None:
        """Parse metadata out of given BAM file"""
        with pysam.AlignmentFile(self.fspath, "rb") as bamf:
            header = bamf.header.to_dict()
            self._parse_read_groups(header)
            self._parse_sort_by(header)
            self._parse_references(header)

    def _parse_read_groups(self, header: dict[str, Any]) -> None:
        """Parse read groups from the header"""
        self.read_groups = header.get("RG", [])

    def _parse_sort_by(self, header: dict[str, Any]) -> None:
        """Parse sort_by information from the header"""
        hd = header.get("HD", {})
        self.sort_by = hd.get("SO", "")
        if self.sort_by == "coordinate":
            self._check_bai()

    def _parse_references(self, header: dict[str, Any]) -> None:
        """Parse references from the header"""
        sqs = header.get("SQ", [])
        if not sqs:
            raise IndexError("No sequence information found in the header")

        self.references = [
            {s["SN"]: int(s["LN"])} for s in sqs if "SN" in s and "LN" in s
        ]

    def _check_bai(self) -> None:
        """Check if index file exists for the coordinate-sorted BAM"""
        bam = parse_path(self.fspath)
        bai = bam.parent / f"{bam.name}.bai"
        if not bai.exists():
            bai = bam.with_suffix(".bai")
            if not bai.exists():
                raise FileNotFoundError(
                    f"Cannot find the index file for the given BAM {bam}"
                )

    def __repr__(self) -> str:
        """Return BAMetadata object representation"""
        return (
            f"BAM file: {self.fspath}\n"
            f"Sort by: {self.sort_by}\n"
            f"# references: {len(self.references)}\n"
            f"# read groups: {len(self.read_groups)}\n"
        )

    def seqnames(self) -> list[str]:
        return [k for r in self.references for k in r.keys()]


@dataclass(slots=True)
class AlignmentRecord:
    rname: str
    mapq: int
    start: int
    end: int
    is_proper: bool
    is_mapped: bool
    is_secondary: bool
    is_duplicate: bool
    n_sc: int
    n_indels: int
    n_ecnt: int
    qname: str
    seq: str
    qual: MutableSequence[int]


# TODO: perhaps return a dictionary is better
# def parse_cigar(cigar: str) -> list[tuple[int, str]]:
#     cigar_iter = itertools.groupby(cigar, lambda k: k.isdigit())
#     cigar_list = []
#     for _, n in cigar_iter:
#         op = int("".join(n)), "".join(next(cigar_iter)[1])
#         cigar_list.append(op)
#
#     return cigar_list


# def parse_md(md_str: str) -> list[str]:
# TODO: check if given md string is empty
#     md_iter = itertools.groupby(
#         md_str, lambda k: k.isalpha() or not k.isalnum()
#     )
#     return ["".join(group) for c, group in md_iter if not c or group]


def count_soft_clip_bases(cigar: str) -> int:
    pattern = re.compile(
        "^(?:(?P<ls>[0-9]+)S)?(?:[0-9]+[MIDNHP=X])+(?:(?P<rs>[0-9]+)S)?$"
    )
    m = pattern.search(cigar)
    if m is None:
        return 0
    n_sc = 0
    if m.group("ls") is not None:
        n_sc += int(m.group("ls"))
    if m.group("rs") is not None:
        n_sc += int(m.group("rs"))

    return n_sc


def count_indel_events(cigar: str) -> int:
    """Count the number of Is and Ds event in the given cigar string"""
    cigar_iter = itertools.groupby(cigar, lambda k: k.isalpha())
    aln_events = ["".join(grpv) for grpk, grpv in cigar_iter if grpk]
    return len([e for e in aln_events if e in ["I", "D"]])


def count_mismatch_events(md: str) -> int:
    md_iter = itertools.groupby(md, lambda k: k.isalpha() or not k.isalnum())
    aln_events = ["".join(grpv) for grpk, grpv in md_iter if grpk]
    return len([e for e in aln_events if e.isalpha()])


def parse_alignments(
    bam_fspath: str,
    rname: str,
) -> Generator[AlignmentRecord, None, None]:
    with pysam.AlignmentFile(bam_fspath, "rb") as bamf:
        for aln in bamf.fetch(contig=rname):
            # Skip alignment record with qc_fail and supplementary marked
            if aln.is_qcfail or aln.is_supplementary:
                continue

            n_sc = 0
            n_ecnt = 0
            n_indels = -1
            if aln.cigarstring is not None:
                # TODO: cigar string parsed twice, wasted
                n_indels = count_indel_events(cigar=aln.cigarstring)
                n_ecnt += n_indels
                n_sc = count_soft_clip_bases(cigar=aln.cigarstring)
            n_mms = -1
            if aln.has_tag("MD"):
                n_mms = count_mismatch_events(md=str(aln.get_tag("MD")))
                n_ecnt += n_mms

            yield AlignmentRecord(
                rname=aln.reference_name
                if aln.reference_name is not None
                else "",
                mapq=aln.mapping_quality,
                start=aln.reference_start,
                end=aln.reference_end if aln.reference_end is not None else -1,
                is_proper=aln.is_proper_pair,
                is_mapped=not aln.is_unmapped,
                is_secondary=aln.is_secondary,
                is_duplicate=aln.is_duplicate,
                n_sc=n_sc,
                n_indels=n_indels,
                n_ecnt=n_ecnt,
                qname=aln.query_name if aln.query_name is not None else "",
                seq=aln.query_sequence
                if aln.query_sequence is not None
                else "",
                qual=aln.query_qualities
                if aln.query_qualities is not None
                else [],
            )
