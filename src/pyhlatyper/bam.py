import itertools
import re
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any, Union

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


def parse_cigar(
    cigar: str,
) -> list[tuple[int, str]]:
    if not cigar:
        raise ValueError("Cannot parse empty CIGAR string")

    cigar_iter = itertools.groupby(cigar, lambda k: k.isdigit())
    cigar_parsed = [
        (int("".join(n)), "".join(next(cigar_iter)[1])) for _, n in cigar_iter
    ]
    if not cigar_parsed:
        raise ValueError(
            f"CIGAR string {cigar} failed to be parsed. Empty list returned"
        )
    return cigar_parsed


def parse_md(md: str) -> list[str]:
    if not md:
        raise ValueError("Cannot parse empty MD string")

    md_iter = itertools.groupby(md, lambda k: k.isalpha() or not k.isalnum())
    md_parsed = ["".join(group) for c, group in md_iter if not c or group]
    if not md_parsed:
        raise ValueError(
            f"MD string {md} failed to be parsed. Empty list returned"
        )
    return md_parsed


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


def count_unaligned_events(
    cigar: Union[str, Sequence[tuple[int, str]]],
) -> int:
    if isinstance(cigar, str):
        cigar = parse_cigar(cigar)
    elif isinstance(cigar, list):
        pass
    else:
        raise TypeError(
            "cigar parameter must be either str or list[tuple[int, str]] type"
        )
    return len([op for _, op in cigar if op in ["I", "D", "S"]])


def count_indel_events(cigar: str) -> int:
    """Count the number of Is and Ds event in the given cigar string"""
    cigar_iter = itertools.groupby(cigar, lambda k: k.isalpha())
    aln_events = ["".join(grpv) for grpk, grpv in cigar_iter if grpk]
    return len([e for e in aln_events if e in ["I", "D"]])


def count_mismatch_events(md: Union[str, Sequence[str]]) -> int:
    if isinstance(md, str):
        md = parse_md(md)
    elif isinstance(md, list):
        pass
    else:
        raise TypeError(
            "md parameter must be either str or list[tuple[int, str]] type"
        )
    # md_iter = itertools.groupby(md, lambda k: k.isalpha() or not k.isalnum())
    # aln_events = ["".join(grpv) for grpk, grpv in md_iter if grpk]
    return len([e for e in md if e.isalpha()])
