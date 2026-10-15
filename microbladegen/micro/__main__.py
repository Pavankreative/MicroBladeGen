"""``python -m microbladegen.micro`` - list, summarise and export the bundled micro cases."""
import argparse
import sys
import warnings

from microbladegen.micro.cases import list_cases, build_case
from microbladegen.micro.summary import case_summary
from microbladegen.micro.export_cfd import export_case


def main(argv=None) -> int:
    warnings.filterwarnings("ignore")
    p = argparse.ArgumentParser(prog="python -m microbladegen.micro", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="list bundled cases")
    s = sub.add_parser("summary", help="build a case and print its summary tables")
    s.add_argument("case", help="case name or path to a case .json")
    e = sub.add_parser("export", help="build a case and write a CFD-ready folder")
    e.add_argument("case")
    e.add_argument("outdir")
    a = p.parse_args(argv)
    if a.cmd == "list":
        print("\n".join(list_cases()))
    elif a.cmd == "summary":
        print(case_summary(build_case(a.case)))
    else:
        man = export_case(build_case(a.case), a.outdir)
        print(f"Wrote {a.outdir}/ (manifest.json, sections/, STL, curves, summaries); periodic sector {man['periodic_sector_deg']:.2f} deg")
    return 0


if __name__ == "__main__":
    sys.exit(main())
