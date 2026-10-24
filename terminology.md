Here's the guide boiled down, in my own words, grouped by its main parts. Examples that aren't cited are mine, not the guide's.

## What EO is

Earth observation means collecting data about Earth's surface and atmosphere with remote sensors, usually as digital images. Planes give the finest detail but are expensive, drones are cheaper but limited by regulation and payload, and satellites give reliable, repeatable global coverage.

## Passive sensors (they only "listen")

These record reflected sunlight or heat emitted by the Earth, so most depend on daylight and are blocked by clouds.

| Term                     | Meaning                                                                 | Typical use                                                          |
| ------------------------ | ----------------------------------------------------------------------- | -------------------------------------------------------------------- |
| **Panchromatic**         | One wide band, shown in greyscale; very sharp                           | Mapping buildings and roads                                          |
| **Thermal infrared**     | Measures emitted heat, so it works at night                             | Sea temperature, urban heat, fires                                   |
| **Multispectral**        | Several narrow bands (e.g. red, green, blue, near-infrared); less sharp | Crop and forest health, floods, burn scars                           |
| **False colour**         | Showing invisible bands as visible colours                              | Healthy vegetation reflects near-infrared strongly, so it stands out |
| **Pan-sharpening**       | Merging a sharp greyscale image with a coarser colour one               | Sharp colour imagery, as in online maps                              |
| **Hyperspectral**        | Hundreds of contiguous bands, giving a full spectrum per pixel          | Identifying minerals, materials, pollution                           |
| **Microwave radiometry** | Passive microwave measurement                                           | Water vapour, soil moisture, ice                                     |

## Active sensors (they send a signal and measure the echo)

These work day or night.

| Term                               | Meaning                                 | Typical use                                                                     |
| ---------------------------------- | --------------------------------------- | ------------------------------------------------------------------------------- |
| **SAR** (synthetic aperture radar) | Radar imaging that sees through cloud   | Ship and oil-spill detection, sea ice, forests, soil moisture                   |
| **InSAR** (interferometry)         | Comparing the phase of two radar images | Ground movement down to centimetres or millimetres, e.g. subsidence or glaciers |
| **Lidar**                          | Same idea with laser light              | Terrain height, tree height, clouds, wind                                       |
| **Radar altimetry**                | Timing a pulse straight down            | Sea level, wave height, ice sheets                                              |
| **GNSS-R**                         | Using reflected GPS/Galileo signals     | Sea wind, soil humidity, water levels                                           |
| **Scatterometry**                  | Radar measuring sea-surface roughness   | Ocean wind for weather forecasts                                                |

## Other measurements

- **Atmospheric chemistry:** each gas absorbs and emits at characteristic wavelengths, which act as a fingerprint. Used for ozone, CO₂, methane, and air quality.
- **Nadir vs limb viewing:** looking straight down gives good horizontal detail; looking sideways through the atmosphere gives good vertical detail.
- **Gravity and the geoid:** satellites sense tiny gravity variations. Used for ocean currents, ice-sheet thickness, and freshwater movement.
- **3D models:** built from stereo images, SAR, or lidar. A **DSM** includes buildings and trees; a **DTM** is bare ground; **DEM** is used loosely for either.
- **Processing levels:** raw data is corrected for sensor errors and atmosphere, and **ortho-rectified** (geometric distortion removed) so it lines up with a map.

## Parameters for choosing imagery

- **Spatial resolution (GSD):** the ground distance one pixel covers. Above 300 m is low, 300–30 m medium, 30–5 m high, under 5 m very high. Weather needs low; precision farming needs high.
- **Swath:** the width of the strip a sensor images. Finer resolution usually means a narrower swath, so more images and more cost.
- **Revisit time:** how long until the same place can be seen again. Flood response needs it short; crop monitoring needs it regular.
- **Orbits:** **sun-synchronous** orbits (about 700 km up) pass at the same local solar time, so lighting is consistent. **Geostationary** satellites (about 36,000 km) stare at one area constantly but at coarse resolution, which suits weather.
- **Agility and constellations:** satellites that can tilt, or several working together, shorten the revisit time.
- **Bit depth:** how many brightness levels are recorded per pixel. 8-bit gives 256, 12-bit gives 4,096.
- **Off-nadir angle:** how far from straight down the image was taken. Large angles make tall buildings appear to lean.
- **Sun elevation and cloud cover:** a low sun means long shadows, and clouds block optical images.

## Getting and paying for data

- **Archive vs tasking:** buying existing images is cheapest; ordering a new acquisition costs more and can take weeks.
- **Free vs commercial:** Copernicus Sentinel and Landsat data are free and open; very-high-resolution imagery is mostly commercial. Price depends on resolution, area, urgency, and processing.
- **Trends:** many small, cheap satellites instead of a few big ones, faster revisit, and selling insights rather than raw images. The guide's example is counting cars in a mall car park to report how busy it is.

The guide is from 2020, so some access links are stale. As far as I know, the Copernicus Open Access Hub it mentions has been replaced by the Copernicus Data Space Ecosystem, which is where you'd get free Sentinel data now.

Sources:

- [Newcomers Earth Observation Guide (ESA Space Solutions)](https://business.esa.int/newcomers-earth-observation-guide)
