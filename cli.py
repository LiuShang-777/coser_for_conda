from __future__ import annotations

import argparse
import sys
from typing import List, Optional


def _parse_index_col(value: str):
    if value.lower() in {"none", "false", "-1"}:
        return None
    try:
        return int(value)
    except ValueError:
        return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="coser",
        description="Calculate CoSEr coupling scores for two traits.",
    )
    subparsers = parser.add_subparsers(dest="mode")

    def add_common_args(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument("-i", "--input", required=True, help="Input phenotype CSV/TSV file.")
        subparser.add_argument("-x", "--x-trait", required=True, help="Trait column for x axis.")
        subparser.add_argument("-y", "--y-trait", required=True, help="Trait column for y axis.")
        subparser.add_argument("-o", "--output", required=True, help="Output CSV file.")
        subparser.add_argument("--sep", default=",", help="Input delimiter. Default: ','.")
        subparser.add_argument("--index-col", default="0", help="Index column name/number, or 'none'. Default: 0.")
        subparser.add_argument("--drop-first-trait-col", action="store_true", help="Drop the first trait column before analysis.")
        subparser.add_argument("--no-scale", action="store_true", help="Do not apply MinMax scaling.")
        subparser.add_argument("--no-fillna", action="store_true", help="Do not fill missing values with column means.")

    angle_parser = subparsers.add_parser("angle", help="Angle-based CoSEr score: 1 - sin(theta).")
    add_common_args(angle_parser)

    distance_parser = subparsers.add_parser("distance", help="Distance-based CoSEr score.")
    add_common_args(distance_parser)
    distance_parser.add_argument(
        "--no-normalize-distance",
        action="store_true",
        help="Output raw distances without normalized distance score.",
    )

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.mode is None:
        parser.print_help()
        return 1

    try:
        if args.mode == "angle":
            from .core import calculate_angle_coupling_score

            calculate_angle_coupling_score(
                pheno_file=args.input,
                x_trait=args.x_trait,
                y_trait=args.y_trait,
                output_file=args.output,
                sep=args.sep,
                index_col=_parse_index_col(args.index_col),
                drop_first_trait_col=args.drop_first_trait_col,
                scale_data=not args.no_scale,
                fillna=not args.no_fillna,
            )
        elif args.mode == "distance":
            from .core import calculate_distance_coupling_score

            calculate_distance_coupling_score(
                pheno_file=args.input,
                x_trait=args.x_trait,
                y_trait=args.y_trait,
                output_file=args.output,
                sep=args.sep,
                index_col=_parse_index_col(args.index_col),
                drop_first_trait_col=args.drop_first_trait_col,
                scale_data=not args.no_scale,
                fillna=not args.no_fillna,
                normalize_distance=not args.no_normalize_distance,
            )
        else:
            parser.error(f"Unknown mode: {args.mode}")
    except Exception as exc:
        print(f"coser: error: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
