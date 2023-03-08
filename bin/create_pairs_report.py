#!/usr/bin/env python
"""Create pairs report."""

from collections import defaultdict
from typing import List, Tuple

import pandas as pd
import panel as pn
import typer
import hvplot.pandas  # noqa

pn.extension()

PAIR_TYPES = {
    "W": "walk",
    "N": "null",
    "X": "corrupt",
    "M": "multi",
    "R": "rescued",
    "U": "unique",
    "D": "duplicate",
}
# https://github.com/4dn-dcic/pairsqc/blob/master/pairsqc.py
ORI_NAMES = dict(zip(["+-", "-+", "++", "--"], ["Inner", "Outer", "Right", "Left"]))


# %%
def _parse_totals_table(data=List[Tuple[str, str]]):
    """Parse totals table."""
    res = []
    total = 0
    for key, val in data:
        key = key.strip()
        if key == "total":
            section = "all"
            total = int(val)
        elif key in ("total_unmapped", "total_single_sided_mapped", "total_mapped"):
            section = "mapping"
        elif key in ("total_dups", "total_nodups"):
            section = "duplicates"
        elif key in ("cis", "trans"):
            section = "cis/trans"
        elif key.startswith("cis_"):
            section = "distance"
        else:
            raise ValueError(f"#{key}#")

        res.append((section, key, int(val)))
    df = pd.DataFrame(res, columns=["Section", "Type", "Count"])
    df["Perc. of Total"] = df["Count"] / total * 100.0
    df["Perc. of Section"] = df.groupby("Section")["Count"].transform(
        lambda x: 100.0 * x / x.sum()
    )
    return df


def _parse_pair_types(data=List[Tuple[str, str]]):
    """Parse pair types."""
    res = []
    for code, val in data:
        left, right = code[0], code[1]
        label = f"{PAIR_TYPES[left]}-{PAIR_TYPES[right]}"
        res.append((code, left, right, label, int(val)))
    df = pd.DataFrame(res, columns=["code", "left", "right", "label", "pairs"])
    df["perc"] = 100.0 * df["pairs"] / df["pairs"].sum()
    return df


def _parse_chrom_freq(data=List[Tuple[str, str]]):
    """Parse chrom freq."""
    res = []
    for code, val in data:
        chr1, chr2 = code.split("/")
        res.append((chr1, chr2, int(val)))

    df = (
        pd.DataFrame(res, columns=["chrom1", "chrom2", "count"])
        .set_index(["chrom1", "chrom2"])
        .sort_index()
        .unstack(fill_value=0)
    )
    df = df.xs("count", axis=1)
    return df


def _parse_summary(data=List[Tuple[str, str]]):
    """Parse summary."""
    res = []
    for key, val in data:
        res.append({"statistic": key, "value": float(val)})
    return pd.DataFrame(res)


def _parse_dist_freq(data=List[Tuple[str, str]]):
    """Parse dist freq."""
    res = []
    for key, val in data:
        interval, ori = key.split("/")
        interval = interval.strip()
        if interval.endswith("+"):
            bin_left = bin_right = interval[:-1]
        else:
            bin_left, bin_right = interval.split("-")
        res.append(
            (int(bin_left), int(bin_right), ori, ORI_NAMES[ori] + f" ({ori})", int(val))
        )
    res = pd.DataFrame(
        res, columns=["bin_left", "bin_right", "ori", "ori_name", "count"]
    )
    return res


def read_pairs_stats(path):
    """Read Pairs stats."""
    _data = defaultdict(list)
    with open(path) as f:
        for i in f:
            if "/" not in i:
                table = "totals"
            else:
                table, i = i.split("/", 1)
            _data[table].append(tuple(i.strip().split("\t")))
    totals = _parse_totals_table(_data["totals"])
    pair_types = _parse_pair_types(_data["pair_types"])
    chrom_freq = _parse_chrom_freq(_data["chrom_freq"])
    summary = _parse_summary(_data["summary"])
    dist_freq = _parse_dist_freq(_data["dist_freq"])
    return totals, pair_types, chrom_freq, summary, dist_freq


def main(pair_stats, report_html, show_chroms=None):
    """Entry point."""
    totals, pair_types, chrom_freq, summary, dist_freq = read_pairs_stats(
        pair_stats)
    totals_pane = pn.Column(
        pn.Row(
            pn.pane.DataFrame(totals.set_index(
                ["Section", "Type"]), width=600),
            totals.query("Section == 'mapping'").hvplot.bar(
                x="Section",
                y="Perc. of Total",
                by="Type",
                hover_cols=["Count", "Perc. of Total"],
                stacked=True,
                width=400,
                title="Mapping Rate",
            ),
        ),
        totals.query("Section == 'distance'").hvplot.bar(
            x="Type", y="Perc. of Section",
            title="Genomic Distance Distribution"
        ),
    )

    pair_type_pane = pn.Column(
        pair_types.hvplot.bar(
            x="label", y="perc", hover_cols=["pairs"], title="Pair Types"
        ),
        pn.pane.DataFrame(pair_types, width=600),
    )

    if not show_chroms:
        show_chroms = chrom_freq.columns

    chrom_freq = chrom_freq.reindex(columns=show_chroms, fill_value=0)
    chrom_contact_pane = pn.Row(
        chrom_freq.loc[show_chroms, show_chroms].hvplot.heatmap(
            width=600,
            height=600,
            colorbar=False,
            rot=45,
            colormap="viridis",
            title="Contact Count",
        ),
        chrom_freq.loc[show_chroms, show_chroms]
        .pipe(lambda x: x.div(x.sum(axis=0), axis=1))
        .hvplot.heatmap(
            width=600,
            height=600,
            colorbar=False,
            rot=45,
            colormap="viridis",
            title="Contact Proportion (normalized by Chromosome)",
        ),
    )

    distance_pane = pn.Row(
        dist_freq.hvplot.line(
            x="bin_right", by="ori_name", y="count", logx=True)
    )

    report = pn.Tabs(
        ("Pairs", totals_pane),
        ("Pair Types", pair_type_pane),
        ("Chrom Contacts", chrom_contact_pane),
        ("Distance", distance_pane),
    )
    report.save(report_html)


if __name__ == "__main__":
    typer.run(main)
