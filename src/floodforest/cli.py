"""Command line entry point.

Built so far: fetch, label, transects, profiles, table, the validation
sampler, the validation-edge KML converter, and qa_centerline.
The later analysis steps are not.
"""

from __future__ import annotations

import argparse

STEPS = ["fetch", "label", "transects", "profiles", "table", "stats",
         "validate", "placebo", "report", "all", "qa_centerline", "dw_gap_check",
         "edges"]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="floodforest")
    parser.add_argument("step", choices=STEPS)
    parser.add_argument("--event", required=True,
                        choices=["melamchi2021", "lende2025", "bhotekoshi2026"])
    parser.add_argument("--only", action="append", default=None,
                        choices=["s1_pre", "copernicus_dem", "s1_near",
                                 "s1_geometry"],
                        help="Refresh these products and leave the others in place.")
    parser.add_argument("--kml", default=None,
                        help="Placemark KML for edges. Names are station_id_side. "
                             "Description lines: image_date: YYYY-MM-DD, "
                             "confidence: 1, 2 or 3, note: text.")
    args = parser.parse_args(argv)
    if args.only and args.step != "fetch":
        raise SystemExit("--only applies to fetch.")
    if args.kml and args.step != "edges":
        raise SystemExit("--kml applies to edges.")
    if args.step == "fetch":
        from .ee_data import fetch_event
        fetch_event(args.event, only=args.only)
        return 0
    if args.step == "label":
        from .label import run_label
        run_label(args.event)
        return 0
    if args.step == "transects":
        from .centerline import run_transects
        run_transects(args.event)
        return 0
    if args.step == "profiles":
        from .profiles import run_profiles
        run_profiles(args.event)
        return 0
    if args.step == "table":
        from .table import run_table
        run_table(args.event)
        return 0
    if args.step == "dw_gap_check":
        from .ee_data import dynamic_world_tile_dates
        dynamic_world_tile_dates(args.event)
        return 0
    if args.step == "qa_centerline":
        from .qa import qa_centerline
        qa_centerline(args.event)
        return 0
    if args.step == "validate":
        from .validate import run_validation_sample
        run_validation_sample(args.event)
        return 0
    if args.step == "edges":
        if args.kml:
            from .validation_kml import import_edge_kml
            import_edge_kml(args.event, args.kml)
            return 0
        from .validate import import_validation_edges
        import_validation_edges(args.event)
        return 0
    raise SystemExit(
        f"Step '{args.step}' for {args.event} is not built yet."
    )


if __name__ == "__main__":
    main()
