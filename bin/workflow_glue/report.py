"""Create workflow report."""
import json

from aplanat.components import fastcat
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
import pandas

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = WFReport(
        "Workflow Template Sequencing report", "wf-template",
        revision=args.revision, commit=args.commit)

    with open(args.metadata) as metadata:
        sample_details = [
            {
                'sample': d['alias'],
                'type': d['type'],
                'barcode': d['barcode']
            } for d in json.load(metadata)
        ]

    if args.stats:
        report.add_section(section=fastcat.full_report(args.stats))

    section = report.add_section()
    section.markdown('## Samples')
    section.table(pandas.DataFrame(sample_details))

    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))

    # write report
    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument("--stats", nargs='*', help="Fastcat per-read stats file(s).")
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
    return parser
