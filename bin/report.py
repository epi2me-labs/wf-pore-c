#!/usr/bin/env python
"""Create workflow report."""

import argparse
import json

from aplanat.components import fastcat
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
import pandas


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("report", help="Report output file")
    parser.add_argument("summaries", nargs='+', help="Read summary file.")
    parser.add_argument(
        "--metadata", default='metadata.json',
        help="sample metadata")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    args = parser.parse_args()

    report = WFReport(
        "Workflow Template Sequencing report", "wf-template",
        revision=args.revision, commit=args.commit)

    with open(args.metadata) as metadata:
        sample_details = [
            {
                'sample': d['sample_id'],
                'type': d['type'],
                'barcode': d['barcode']
            } for d in json.load(metadata)
        ]

    report.add_section(
        section=fastcat.full_report(args.summaries))

    section = report.add_section()
    section.markdown('## Samples')
    section.table(pandas.DataFrame(sample_details))

    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))

    # write report
    report.write(args.report)


if __name__ == "__main__":
    main()
