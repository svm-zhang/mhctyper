import re
import warnings
from dataclasses import dataclass, field, fields

__VALID_LOCI: list[str] = [
    "A",
    "B",
    "C",
    "E",
    "F",
    "G",
    "H",
    "J",
    "K",
    "L",
    "N",
    "P",
    "R",
    "S",
    "T",
    "U",
    "V",
    "W",
    "X",
    "Y",
    "Z",
    "DRA",
    "DRB1",
    "DRB2",
    "DRB3",
    "DRB4",
    "DRB5",
    "DRB6",
    "DRB7",
    "DRB8",
    "DRB9",
    "DQA1",
    "DQB1",
    "DQA2",
    "DQB2",
    "DQB3",
    "DOA",
    "DOB",
    "DMA",
    "DMB",
    "DPA1",
    "DPB1",
    "DPA2",
    "DPB2",
    "DPA3",
]


@dataclass
class HLAllelePattern:
    prefix: str = field(default="(:?HLA|hla)[-_]")
    locus: str = field(default="[a-zA-Z]+[0-9]?")
    digit_fields: str = field(init=False)
    digit_field_sep: list[str] = field(default_factory=lambda: [":", "_"])
    expr_suffix: list[str] = field(
        default_factory=lambda: ["N", "L", "S", "C", "A", "Q"]
    )
    resolution: int = field(default=4)
    sep: str = field(default="[\\*_]")

    def __post_init__(self) -> None:
        self._validate_resolution()
        # the first line gives "(?:[0-9]+[:_])"
        # the second line gives {0, resolution-1}
        # the last 2 lines give [0-9NLSCAQnlscaq]+
        self.digit_fields = (
            f"(?:[0-9]+[{''.join(self.digit_field_sep)}])"
            f"{{0,{self.resolution-1}}}"
            f"[0-9{''.join(self.expr_suffix)}"
            f"{''.join([c.lower() for c in self.expr_suffix])}]+"
        )
        self._validate_fields()
        self._add_name_to_group_match()

    def _validate_resolution(self) -> None:
        """Make sure resolution is set between 1 and 4"""
        if self.resolution < 1:
            self.resolution = 1
        if self.resolution > 4:
            self.resolution = 4

    def _validate_fields(self) -> None:
        for f in fields(self):  # f is field object of this class
            v = getattr(self, f.name)
            if isinstance(v, str):
                if not v or not v.strip():
                    raise ValueError(
                        f"Property {f.name} cannot be empty or "
                        "just whitespaces. Whitespaces at beginning and "
                        "end of the pattern string are stripped"
                    )

    def _add_name_to_group_match(self) -> None:
        for f in fields(self):
            v = getattr(self, f.name)
            if isinstance(v, str):
                setattr(self, f.name, f"(?P<{f.name}>{v.strip()})")


@dataclass
class HLAllele:
    """
    HLAllele object to represent a HLA allele

    A HLA allele = "{prefix}{locus}{sep}{fields}"
    """

    prefix: str
    locus: str
    digit_field: str
    sep: str

    def __str__(self) -> str:
        return f"{self.prefix}{self.locus}{self.sep}{self.digit_field}"


def _check_if_allele_empty(allele: str) -> None:
    """Check if given allele string is empty"""
    if not allele or not allele.strip():
        raise ValueError(
            "allele cannot be empty or just whitespaces. "
            "Whitespaces at beginning and end of the pattern "
            "string are stripped"
        )


def reduce_resolution(allele: str, ap: HLAllelePattern) -> str:
    """Reduce resolution of input HLA allele"""
    pattern = f"{ap.prefix}?{ap.locus}{ap.sep}{ap.digit_fields}"

    pattern_c = re.compile(pattern)
    m = pattern_c.search(allele)
    if m is None:
        raise ValueError("Failed to decompose the given allele")

    if m.group() == allele:
        warnings.warn(
            f"No resolution reduced for input {allele=}", RuntimeWarning
        )
    return m.group()


def decompose(allele: str, ap: HLAllelePattern) -> HLAllele:
    """Decompose given HLA allele string"""
    _check_if_allele_empty(allele)

    pattern = f"^{ap.prefix}?{ap.locus}{ap.sep}{ap.digit_fields}$"

    pattern_c = re.compile(pattern)
    m = pattern_c.search(allele)
    # quit when no match found
    if m is None:
        raise ValueError("Failed to decompose the given allele")

    prefix = m.group("prefix") if m.group("prefix") is not None else ""

    locus = m.group("locus")
    # make sure identified locus is valid
    if locus is not None and locus.upper() not in __VALID_LOCI:
        raise ValueError

    return HLAllele(
        prefix=prefix,
        locus=m.group("locus"),
        digit_field=m.group("digit_fields"),
        sep=m.group("sep"),
    )
