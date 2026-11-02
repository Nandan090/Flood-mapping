## Clone

Clone the repository and enter its root directory:

````sh

- Windows 10/11 with PowerShell

## Clone and install
## Setup
Open PowerShell and run:
Create and activate a virtual environment, then install the project in editable mode.

Linux/macOS:

```sh
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
````

Windows PowerShell:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e .
```

If PowerShell blocks environment activation, allow local scripts for your user, then activate again:

````powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
.\.venv\Scripts\Activate.ps1
cd Flood-mapping
py -3.11 -m venv .venv
Verify the installation:
python -m pip install --upgrade pip
```sh
python -m floodforest.cli --help
````

$env:PYTHONPATH = "$PWD\src"
Earth Engine is required to fetch remote satellite and terrain data, but not to run the local test suite. With the virtual environment active, authenticate once:

```sh
earthengine authenticate
```

````powershell
Register a Google Cloud project for noncommercial Earth Engine use, then set its project ID in `config.yaml` under `ee_project`. Do not add credentials or tokens to the repository. Before running `fetch`, make sure `data/<event>/centerline.geojson` exists and has been verified in QGIS; the configured centerline path for each event is listed in `config.yaml`.


Earth Engine is needed to fetch remote satellite and terrain data, but not to run the local test suite.
With the virtual environment active, run from the project root:
1. Register a Google Cloud project for noncommercial Earth Engine use and make a note of its project ID.
2. With the virtual environment active, authenticate once:

   ```powershell
   earthengine authenticate
````

3. Open `config.yaml` and set `ee_project` to that registered project ID. Do not add credentials or tokens to the repository.
4. Before fetching, confirm the event's centerline file exists and has been checked in QGIS with `verified_by_user: true`. The Melamchi event is the project's current fetch target.

## Test the installation

From the project folder, with `.venv` active and `PYTHONPATH` set:

```powershell
pytest
```

These tests use small fixtures and made-up data; they do not download data or produce study results. Slow statistical tests are excluded from a normal run and can be selected with `pytest -m slow`.

## Run the implemented pipeline

The current remote-data fetch is implemented for Melamchi 2021 only. Run its steps in order:

```powershell
python -m floodforest.cli fetch --event melamchi2021
python -m floodforest.cli label --event melamchi2021
python -m floodforest.cli transects --event melamchi2021
python -m floodforest.cli profiles --event melamchi2021
python -m floodforest.cli table --event melamchi2021
```

`fetch` requires a working Earth Engine login, a valid project ID, and the verified centerline. The other steps use the cached downloads and intermediate products from the earlier steps. Each step can be rerun independently after its prerequisites have completed.

Other implemented utilities include:

```powershell
python -m floodforest.cli qa_centerline --event melamchi2021
python -m floodforest.cli validate --event melamchi2021
python -m floodforest.cli edges --event melamchi2021 --kml path\to\validation_edges.kml
```

The validation sampler creates locations and transects for manual flood-edge checking. `edges` imports digitized edges from a KML file; see `data/melamchi2021/validation/how_to_digitise.md` for the expected format and workflow. `dw_gap_check` is an additional Earth Engine diagnostic for the configured Melamchi tile.

## Run the placebo (no-flood control)

The placebo (spec 8.9 and 14.2) runs the same pipeline on two pre-event years with no notable flood. It is built for Melamchi 2021:

```powershell
python -m floodforest.cli placebo --event melamchi2021
```

What it does:

- **Checks the period first.** It refuses to run unless `flood_check.status` under the event's `placebo` block in `config.yaml` is `confirmed`. The confirmed interval must also cover the placebo image windows. For Melamchi the interval is 3 Aug 2018 – 15 Dec 2019; see `results/placebo_period_check.md`.
- **Runs the unchanged pipeline.** It runs `fetch`, `label`, `transects`, `profiles` and `table` on a derived event, `melamchi2021_placebo`. That event uses the real event's centerline, reach, manual exclusions and thresholds, with the placebo image windows. Its fake event date is 2019-06-15.
- **Keeps its files apart.** Everything goes to `data/melamchi2021_placebo/` and `outputs/melamchi2021_placebo/`. The real event's files are only read.
- **Reuses downloads.** The first run downloads the placebo rasters, which needs Earth Engine. Later runs reuse them; add `--refetch` to download again.
- **Writes a report.** `outputs/melamchi2021_placebo/placebo_report/PLACEBO_REPORT.md` compares the placebo with the real event. Alongside it are the CSV tables, the figures and contact sheets of everything the placebo label detected.
- **Sorts detections for review.** Each detected patch is placed in a diagnostic category: channel margin, floodplain or terrace, or hillslope. The rules are in `placebo_diagnostics` in `config.yaml`. These categories only help a person review the patches. They never decide what counts as flood impact.

**Forest exposure.**

- The primary placebo forest is built from the placebo "pre" window, Oct–Dec 2018, as 14.2 requires.
- Where Dynamic World has no 2018 image, WorldCover 2020 stands in. This is a placebo-only setting (`worldcover_standin`); the affected rows are flagged, and a sensitivity run leaves them out.
- The real event's Oct–Dec 2020 forest is shown as a separate sensitivity and is never changed.

A single placebo step can be rerun with the derived key, for example `python -m floodforest.cli label --event melamchi2021_placebo`.

The forest-effect half of the placebo test (β1 not distinguishable from zero) needs the `stats` step, which is not built yet.

## What a successful run produces

The pipeline caches downloaded inputs and intermediate files under `data/` and writes event results under `outputs/`. For the implemented Melamchi steps, expect:

- `data/melamchi2021/rasters/`: downloaded satellite and terrain GeoTIFFs.
- `data/melamchi2021/tables/profile_geometry.csv` and `profile_outcomes.csv`: transect geometry, controls, and measured flood-reach outcomes.
- `outputs/melamchi2021/provenance.json`: data sources, scene dates, settings, and processing details.
- `outputs/melamchi2021/label/`: combined and source-specific impact-label rasters, plus QA map sheets.
- `outputs/melamchi2021/transects/`: transect GeoPackage and review sheets.
- `outputs/melamchi2021/transects.csv`, `transects_labels.csv`, `transects.gpkg`, and `transects_columns.csv`: the assembled per-bank analysis table, label-specific outcomes, spatial layers, and column descriptions.

The impact labels distinguish impacted, not impacted, and unknown pixels. The table contains one row per station and riverbank, including reach distance/height, overtopping, forest exposure, terrain controls, and exclusion or censoring details where available.

## Current implementation status

The intended full workflow also includes statistical models and HTML reports. Those steps are not implemented yet: although the CLI lists `stats`, `report`, and `all`, invoking them currently reports that the step is not built. The `placebo` step is built (see above), but without the model part. Fetching is likewise currently restricted to `melamchi2021`; the other event keys are planned/configured but are not yet available for the full run. Do not expect `report.html` or `summary_report.html` from the current implemented pipeline.

For the full analytical specification, event settings, and study limitations, see `CURSOR_BUILD_SPEC.md`, `config.yaml`, and `info.md`.
