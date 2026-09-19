

"""
Data Acquisition Script
Study: Quantifying the buffering role of Himalayan forest cover against
       secondary slope failure following the August 26, 2026 Langtang
       glacial-collapse flood, Trishuli River corridor, Nepal.

Event summary (for reference / methods section):
  - Trigger: ice-rock avalanche / partial glacier collapse, Langtang Lirung,
    Langtang National Park (NOT rainfall — Nepal DHM reported no significant
    antecedent rainfall in Rasuwa immediately before the event).
  - Date: 26 August 2026, ~08:37 local time (02:52:10 UTC).
  - Path: Lhende Khola -> Bhote Koshi -> Trishuli River, ~72-100 km runout.
  - Affected districts: Rasuwa, Nuwakot, Dhading (Bagmati Province).

Author: [Nandan K C]
"""

import ee
import geemap


# 0. INITIALIZE
# ---------------------------------------------------------------------------
ee.Authenticate()
ee.Initialize(project='')  


# ---------------------------------------------------------------------------
# 1. STUDY AREA
# ---------------------------------------------------------------------------
# Original bbox (85.30, 28.20, 85.55, 28.32) covers the upper Rasuwa /
# Langtang collapse zone. The debris/flood corridor travelled ~72-100 km
# downstream (south) along the Trishuli, so the AOI is extended southward
# to fully capture the runout path and secondary slope failures in
# Nuwakot and Dhading districts. Adjust further once you've checked
# imagery coverage.
aoi = ee.Geometry.Rectangle([85.15, 27.75, 85.60, 28.35])

# Keep the original tight box too, in case you want a "source zone only"
# subset analysis (glacier collapse + upper valley).
# Source-zone box refined around the independently reported source centroid
# (85.52701 E, 28.28858 N), with enough buffer to comfortably contain the
# full mapped changed-surface envelope (up to ~1.8 km2, max dimension ~1.5 km).
aoi_source_zone = ee.Geometry.Rectangle([85.505, 28.275, 85.550, 28.305])

# ---------------------------------------------------------------------------
# 2. DATE WINDOWS
# ---------------------------------------------------------------------------
# Pre-event: clear dry-ish window before the collapse.
PRE_START, PRE_END = '2026-07-15', '2026-08-20'

# Post-event: cloud cover will be heavy for days after (monsoon season +
# disturbed atmosphere). Widen the post window and rely more on Sentinel-1
# (cloud-penetrating) for the immediate aftermath; use a later, clearer
# Sentinel-2 window for optical confirmation of the scar/runout extent.
POST_START_S1, POST_END_S1 = '2026-08-26', '2026-09-10'   # immediate, SAR
POST_START_S2, POST_END_S2 = '2026-09-01', '2026-10-15'   # clearer optical

# ---------------------------------------------------------------------------
# 3. SENTINEL-2 OPTICAL (cloud masked)
# ---------------------------------------------------------------------------
def mask_s2_clouds(image):
    qa = image.select('QA60')
    cloud_bit, cirrus_bit = 1 << 10, 1 << 11
    mask = qa.bitwiseAnd(cloud_bit).eq(0).And(qa.bitwiseAnd(cirrus_bit).eq(0))
    return image.updateMask(mask).divide(10000)

import sys
from datetime import datetime, timedelta

def get_s2_composite_robust(start, end, aoi, cloud_thresh=30, label="window",
                             max_extra_days=60, extra_day_step=15):
    """
    Tries the given date range/cloud threshold first. If zero images are
    found, progressively widens the end date (in extra_day_step increments,
    up to max_extra_days) and relaxes the cloud threshold, printing progress
    at every attempt. Raises only if nothing works even at the widest/loosest
    settings, so you get useful imagery instead of a silent crash.
    """
    attempts = [(start, end, cloud_thresh)]
    end_dt = datetime.strptime(end, '%Y-%m-%d')
    for extra in range(extra_day_step, max_extra_days + 1, extra_day_step):
        new_end = (end_dt + timedelta(days=extra)).strftime('%Y-%m-%d')
        attempts.append((start, new_end, min(cloud_thresh + 20, 90)))

    for s, e, ct in attempts:
        collection = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                      .filterBounds(aoi)
                      .filterDate(s, e)
                      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', ct))
                      .map(mask_s2_clouds))
        size = collection.size().getInfo()
        print(f"[{label}] tried {s} to {e}, cloud<{ct}%: {size} image(s) found")
        sys.stdout.flush()
        if size > 0:
            return collection.median()

    raise RuntimeError(
        f"[{label}] No Sentinel-2 imagery found for {start}..{end} even after "
        f"widening to {attempts[-1][1]} and cloud_thresh<={attempts[-1][2]}%. "
        f"The AOI/date combination may need S2_CLOUD_PROBABILITY-based masking "
        f"or a manually chosen clear scene — let me know and I'll rewrite this."
    )

s2_pre = get_s2_composite_robust(PRE_START, PRE_END, aoi, cloud_thresh=30, label="PRE")
s2_post = get_s2_composite_robust(POST_START_S2, POST_END_S2, aoi, cloud_thresh=40, label="POST")

# NDVI (forest/vegetation density proxy)
ndvi_pre = s2_pre.normalizedDifference(['B8', 'B4']).rename('NDVI_pre')
ndvi_post = s2_post.normalizedDifference(['B8', 'B4']).rename('NDVI_post')
ndvi_loss = ndvi_pre.subtract(ndvi_post).rename('NDVI_loss')

# Fractional Vegetation Cover (FVC) from NDVI (Gutman & Ignatov method)
NDVI_SOIL, NDVI_VEG = 0.05, 0.85  # typical bare soil / full veg bounds; tune to region
fvc_pre = (ndvi_pre.subtract(NDVI_SOIL)
           .divide(NDVI_VEG - NDVI_SOIL)
           .clamp(0, 1)
           .rename('FVC_pre'))

# NDWI for flood/water extent confirmation (optical, cloud-permitting)
ndwi_post = s2_post.normalizedDifference(['B3', 'B8']).rename('NDWI_post')

# ---------------------------------------------------------------------------
# 4. SENTINEL-1 SAR (primary disturbance/flood detector — cloud independent)
# ---------------------------------------------------------------------------
s1_collection = (ee.ImageCollection('COPERNICUS/S1_GRD')
                  .filterBounds(aoi)
                  .filter(ee.Filter.eq('instrumentMode', 'IW'))
                  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
                  .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING')))  # keep orbit consistent

s1_pre = s1_collection.filterDate(PRE_START, PRE_END).select('VV').median()
s1_post = s1_collection.filterDate(POST_START_S1, POST_END_S1).select('VV').median()

# Backscatter change: sharp drops = new water/debris smoothing;
# sharp increases = new bare/rough surface (scar, deposited debris)
s1_diff = s1_post.subtract(s1_pre).rename('VV_diff')

# Log-ratio version (often more robust for SAR change detection than simple diff)
s1_log_ratio = s1_post.divide(s1_pre).log().rename('VV_log_ratio')

# ---------------------------------------------------------------------------
# 5. TERRAIN (DEM-derived)
# ---------------------------------------------------------------------------
dem_raw = (ee.ImageCollection('COPERNICUS/DEM/GLO30_2024_1')
           .filterBounds(aoi)
           .select('DEM')
           .mosaic()
           .rename('elevation')
           .clip(aoi))

# Copernicus GLO-30 has known small voids in steep/shadowed/glaciated terrain
# (very relevant here -- Langtang is high-relief, glacier-covered). Filling
# gaps in the DEM itself (before deriving slope/aspect/hillshade) is better
# than filling those derived products afterward, which would smear artificial
# ridgelines. A small focal-mean kernel fills isolated void pixels without
# meaningfully distorting real terrain.
dem_filled = dem_raw.focal_mean(radius=2, kernelType='square', units='pixels')
dem = dem_raw.unmask(dem_filled)  # use original values where present, fill only voids

slope = ee.Terrain.slope(dem).rename('slope')
aspect = ee.Terrain.aspect(dem).rename('aspect')
hillshade = ee.Terrain.hillshade(dem).rename('hillshade')

# ---------------------------------------------------------------------------
# 6. FLOW-CORRIDOR DISTANCE (NEW — key feature for this event framing)
# ---------------------------------------------------------------------------
# Distance from the actual Lhende Khola -> Trishuli flow path, not generic
# stream proximity. Use HydroSHEDS/JRC river network or digitize the
# known corridor centerline manually if the network dataset doesn't
# resolve the Lhende Khola at this scale.
#
# Approach: use JRC Global Surface Water occurrence to derive a river mask,
# or (recommended) manually digitize the corridor centerline in QGIS from
# the post-event imagery, upload as an asset, and compute distance from it.
# Placeholder using HydroSHEDS flow accumulation as a proxy:
flow_acc = ee.Image('WWF/HydroSHEDS/15ACC').select('b1').rename('flow_accumulation')
river_mask = flow_acc.gt(1000)  # threshold tuned to isolate major channels; adjust
dist_from_river = (river_mask.fastDistanceTransform(256)
                    .sqrt()
                    .multiply(ee.Image.pixelArea().sqrt())  # convert to meters
                    .rename('dist_from_channel'))

# TWI must be computed locally (SAGA GIS / richdem / pysheds) after DEM export —
# GEE has no native TWI function. See 02_twi_local.py.

# ---------------------------------------------------------------------------
# 7. RAINFALL (included as a CONTROL variable to demonstrate it was NOT
#    the trigger — expect near-zero feature importance; this is a citable
#    finding, not a data-quality problem)
# ---------------------------------------------------------------------------
rainfall_pre_event = (ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
                       .filterBounds(aoi)
                       .filterDate('2026-08-19', '2026-08-26')
                       .sum()
                       .rename('rainfall_7day_pre'))

# ---------------------------------------------------------------------------
# 8. STACK FEATURES AND EXPORT
# ---------------------------------------------------------------------------
# NOTE: geemap.ee_export_image() uses GEE's direct-download endpoint, which
# is capped at 50 MB per request. Your AOI x 13 bands x 10m resolution
# produces a multi-GB request, so we use ee.batch.Export.image.toDrive()
# instead — this has no size cap and runs as an async task on Google's
# servers. You then download the finished .tif from your Google Drive.

feature_stack = (ndvi_pre
                  .addBands([ndvi_post, ndvi_loss, fvc_pre, ndwi_post,
                             s1_diff, s1_log_ratio,
                             dem, slope, aspect, hillshade,
                             dist_from_river, rainfall_pre_event])
                  .toFloat())  # GEE export requires all bands share one dtype;
                               # our stack mixes Float32/Float64 sources, so
                               # cast everything to Float32 explicitly.

def export_to_drive(image, description, region, folder='nepal_flood_research', scale=10):
    task = ee.batch.Export.image.toDrive(
        image=image,
        description=description,
        folder=folder,
        fileNamePrefix=description,
        region=region,
        scale=scale,
        maxPixels=1e10,  # AOI is large; raise the default pixel cap
        fileFormat='GeoTIFF'
    )
    task.start()
    print(f"Started export task: '{description}' -> Google Drive folder '{folder}/'")
    return task

task_corridor = export_to_drive(feature_stack, 'feature_stack_corridor', aoi)
task_source_zone = export_to_drive(feature_stack, 'feature_stack_source_zone', aoi_source_zone)

print("""
Export tasks submitted to Earth Engine's servers (this runs remotely, not on
your machine). Track progress at: https://code.earthengine.google.com/tasks

Once each task shows COMPLETED there, the .tif files will appear in:
  Google Drive > nepal_flood_research/

From there, download them to your local machine, e.g.:
  ~/Desktop/Research paper/feature_stack_corridor.tif
  ~/Desktop/Research paper/feature_stack_source_zone.tif

Tasks can take anywhere from a few minutes to ~30+ minutes depending on
AOI size and server load -- this is normal for exports this large.
""")