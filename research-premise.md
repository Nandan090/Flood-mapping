## The goal in plain words

The paper asks: **when a huge debris flood comes down a Himalayan valley, does forest along the riverbanks reduce how far and how high the flood spreads sideways?**

Forests are widely said to protect against floods, but the evidence comes mostly from small events. This needs to be tested on a very large, sediment-heavy floods in Nepal, using satellite data only.

## What's being tested

**Hypothesis:** riverbanks that were forested before the flood show less sideways reach of the flow than similar banks without forest.

**Three floods, same method on each:**

- Bhotekoshi 2026 (the main test, about 72 km of river)
- Lende Khola 2025 (same river, smaller flood)
- Melamchi 2021 (a different valley, as an independent check)

**The method (transects):** draw a line across the valley every 50–100 m, at right angles to the river. For each line and each bank, record three kinds of numbers:

| Type         | What it is                                              | When it's measured |
| ------------ | ------------------------------------------------------- | ------------------ |
| **Outcome**  | How far and how high the debris reached                 | After the flood    |
| **Exposure** | How much of the bank strip (50/100/200 m) was forest    | Before the flood   |
| **Controls** | Channel width, slope, valley width, bend, upstream area | From terrain       |

Then a statistical model asks whether more forest goes with less reach once the controls are held equal. That gives roughly 700 lines and 1,400 rows for the main event.

**The four traps the paper is designed around:**

1. Forest grows on higher ground the flood wouldn't reach anyway, so forest and non-forest must be compared at the same height and distance from the river.
2. Damage can't be measured as vegetation loss, because only forest can lose vegetation. The flood footprint must come from radar and bare-soil signals instead.
3. Cloud covered the flood days, and the gorge is narrow, so 10 m pixels and radar are needed.
4. Neighbouring pixels aren't independent evidence, hence transects rather than pixels.

## Data you'd likely retrieve

| Need                         | Dataset                                             | Notes                                                                             |
| ---------------------------- | --------------------------------------------------- | --------------------------------------------------------------------------------- |
| Flood footprint (radar)      | Sentinel-1 SAR                                      | Before/after backscatter change; sees through cloud                               |
| Flood footprint (optical)    | Sentinel-2                                          | Brightness and bare-soil change, Oct–Dec composites; NDVI for pre-event greenness |
| Terrain                      | Copernicus DEM (30 m) or similar                    | Height above channel, slope, valley width                                         |
| River line and upstream area | HydroSHEDS / MERIT Hydro, or derived from the DEM   | Needed to place transects                                                         |
| Forest cover and past loss   | ESA WorldCover, Hansen Global Forest Change         | Pre-event only                                                                    |
| Canopy height                | GEDI-based global canopy height maps                | For the mean tree height column                                                   |
| Validation                   | Published extent maps (ICIMOD, Vantor open imagery) | To check your own flood mapping                                                   |

**Where the API work is:** Google Earth Engine (Python API) holds almost all of this in one place and is the fastest route. The alternative is the Copernicus Data Space Ecosystem with a STAC API. I haven't checked current availability of each dataset, so confirm before you commit.

**Python stack:** `geopandas` and `shapely` (transect geometry), `rasterio`/`xarray` (image arrays), `statsmodels` (the regression), plus QGIS for visually checking results.

## What to focus your reading on

1. **SAR basics and change detection**, including layover and shadow in steep terrain. This is the hardest part and the paper leans on it.
2. **Sentinel-2 cloud masking and compositing**, which builds on the bands and indices we covered.
3. **Terrain analysis:** DEMs and "height above nearest drainage" (HAND), which is essentially the paper's main outcome measure.
4. **Generating river centerlines and perpendicular transects**, a pure geometry and programming problem that suits you well.
5. **Lighter touch on the statistics:** confounding, matching, block bootstrap, E-values. Your friend probably owns this, but you need to know what the columns feed into.

## Things to ask your friend

- **Does code already exist?** The text refers to "the E-value already in the pipeline", "the script" and a README. Either there's a codebase you should see, or those passages describe code nobody has written yet.
- **Which DEM gives reach height?** Free global DEMs are 30 m, which is coarse for a gorge about 30 m wide. This is the biggest technical risk I see, since reach height is the primary outcome.
- **What's settled?** Belt width and height-versus-distance are still listed as open choices, and the introduction and title sections are empty.

I haven't verified the event details in the table (dates, volumes, casualties); I've taken them as written.

Absolutely. These concepts fit together into one workflow: **detect where a flood/debris event happened → measure the terrain and vegetation → compare affected and unaffected areas → determine whether the observed changes are actually attributable to the event.**

---

# 1. SAR basics and change detection

**SAR** = **Synthetic Aperture Radar**.

Unlike an ordinary satellite photograph, SAR does **not need sunlight** and can generally see through clouds. It sends microwave radar pulses toward the ground and measures the signal that comes back.

![Image](https://images.openai.com/static-rsc-4/EUMxDkELO3LG1EM2f2DANFUg1pR-2PR1RA4OuNxLOs3UuLBiQcUmeT889ZLcNuzv0_ZFSAaaUTAxYu0S1z4pZ5wUaN-HJXLkzHI4IokPpv4Lc1ERCzzEtZMKlTkya9DLIbqE2xkZpOzHVc0SjHrmA_PXaP-nlIlSfy5D5J-vBUOIKCo89LTWkFeHAjy_V3T4?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/_7CYcl6wr6eZJ24Aq4g07ER0BYegsb1O1ZW-ufnjlrpTZiINwVlA5mVGrz4cCTanHCn1BBg24ckSHaBO58oUNSWKEajY823qkCTMUj4Uuv12KERNH4ydiqbwh04y-Grfn35KT3CVuba1SBcNTJmutsxYyP9VhRBNdHti7dDRt2jEhFnE003kKWwhhTQWbF5a?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/bL8HCBhNxlV_0P1TfkdFLlotxeUof_QtALnxRQKYR4ofDonmPfXR-8luveiO9vphHd7eFNrX_J6IN2bVBbBA6WhguMLxtCdF-j8O_JqrVVlywcX_11TwA_3dih8AtDlOwHlqBje3x5cM7dXZGp9tqhkeABby-ANtQLPI5BA6RjEINeBAWSTZrzSygtGdeS-b?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/EutzZSTXlskR16ApnQdlQVxst6o7tzDBt001cBXlyfDcD1E4yku9rEFATIlMRK7dOxIEoTeFA758a8WBJBF4NbOfq0MNkOVk1obvrnQ-5UIY6rUTYJIJmHxsjfFdou-GuxlmbuPGNwDAjEa9hySGdW43qb6olaeiE7AEImAQQjbWNySSxUd7dFueSenYL3Dk?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/vGlD17NgSc0ZYpay_Jvog0oqqoZT65Ziw35PuFDkDbeVasjI1_9eLmi3TsMGoliuuZmKCXMvMXnRlw4spOAAOoDlJst93cJ8jB8KzaTpKVR53lqKfjExWPTs8GkHJzFc0OAtaPInD1geDCTmmfpMfW5CGb85yJvwOeWwM6QA9yiLbx7FGZvQKTCJ4rVDdGUg?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/DLUT6v8gGggUcrqz7nUy1AVIhWI5nfDi0dUkbE0daShSWQpyE4WX6pMhc3qCdxMNlda4-8eQV3aOhs0Zg_E-r1eAlB5ycf4EVKzOPle7RCKyRBHEJqkBxXj_hRL-aYx-C2NeAlu5hLFHRdRdmBtK79SzlxF2BMo_kVe7rEPEuMxWNuO1SDXG0nvwamG1UkJV?purpose=fullsize)

### The basic idea

Imagine shining a flashlight at a mountain.

With an optical camera, you see the mountain's visible surface.

With SAR, you're effectively asking:

> **How strongly does this surface reflect radar energy back toward the satellite?**

Different surfaces return different amounts of radar energy.

| Surface               | Typical SAR response     |
| --------------------- | ------------------------ |
| 🌊 Smooth water       | Very dark                |
| 🌱 Vegetation         | Variable, often moderate |
| 🏙️ Buildings          | Very bright              |
| 🪨 Rough ground/rocks | Moderate to bright       |

The important thing is that SAR measures **backscatter**, often represented as **dB (decibels)**.

---

## SAR change detection

Suppose you have:

**Before event → SAR image → After event → SAR image**

You can compare the two.

For example:

```text
BEFORE                         AFTER

🌿🌿🌿🌿🌿                    🪨🪨🪨🪨🪨
🌿🌿🌿🌿🌿        →           🪨🌊🪨🌊🪨
~~~~~~~~~~~                  ~~~~~~~~~~~
   River                         River
```

If the radar response changes substantially, that can indicate that something on the ground changed.

A simple change measure might be:

$$
\Delta SAR = SAR_{after} - SAR_{before}
$$

Large differences can indicate:

- flooding
- erosion
- deposited sediment
- vegetation damage
- landslides
- changes in surface roughness

### But here's the difficult part: steep terrain

SAR doesn't see the landscape from directly above.

It looks at it **from the side**.

That creates geometric distortions.

---

# 2. Layover and shadow

These are extremely important when studying **mountain rivers and valleys**.

Imagine a satellite looking toward a mountain:

```text
                 🛰️
                /
               /
          ____/  ← mountain
        /    /
       /    /
______/____/____________
```

Because the radar is looking from an angle, the apparent position of parts of the mountain can be distorted.

## Layover

**Layover** happens when the top of a steep slope is closer to the satellite than the bottom.

The radar signal from the top can arrive **before** the signal from the bottom.

So the radar image effectively squashes the slope.

```text
Real mountain:

       /\
      /  \
     /    \
____/      \____


SAR can appear more like:

      /|
     / |
____/  |
```

In extreme cases, features can appear in the **wrong apparent position**.

### Why this matters

Suppose your river is next to a very steep valley wall.

A change in SAR backscatter might not represent an actual flood/debris change.

It might simply be caused by the **geometry of the radar observation**.

---

## Radar shadow

The opposite problem is **radar shadow**.

If a mountain is extremely steep, the radar signal may not reach the surface behind it.

```text
             🛰️
            ↘
             ↘
              /\
             /  \████████
            /    \███████
___________/______\███████
                   ↑
               radar shadow
```

The area behind the mountain receives little or no radar illumination.

It therefore appears **very dark**.

The key lesson is:

> **Dark in SAR does not necessarily mean water or no vegetation. It could simply be radar shadow.**

---

## Why this is especially important for the paper

If the study examines **debris floods in steep river valleys**, SAR is attractive because it works despite clouds and bad weather.

But the researchers have to distinguish:

**real landscape change**

from

**apparent change caused by SAR geometry.**

That is why DEMs and terrain information become important later.

---

# 3. Sentinel-2 cloud masking and compositing

Now we're switching from **radar** to **optical imagery**.

![Image](https://images.openai.com/static-rsc-4/xHtJdlxUouxBMrjZo4JDpgufn2xo43Yg9Wd0BTObVG9v8DncvKfZNBFInAF1V_B29BtkYzXGhJ2yIPT7E1NgbyudwtBozxMPATydPOsCLOBCBBtiyTIomj4gQCnLfn6idCpGRMXuA8XyMAv1ckVpMrP7vTjEjZ4lwx0qRb2XsVX_boC3cNEYecEci1G2jNNe?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/RBojz458tlvnP9y2AXuKoNne188_UTmBB83NTLzR6rTEQWpR-dk4SrmMp_j2PfZD7byHrSWFj4oyTS9Cv5_JAtu5Q7TlFqAWmOigfz6HlybLP_sgf8UKv-MkANs7bpaVki5uuyg2Om-fBGtl1dAjypuXat6G_FrVCdQ2SrUg-8wBOC0U4pzBaexaGuflv4em?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/NZy8iHKhtFosNRNBHl4xU-SQPaG_jJUG-zqE1wV5n5ILKK1Gyid8l-bR0uSih8-lesCaJjDssTDBctO4ARfQLqwL1OkhUEk2DPnCf8bATDxneFpi6MPsmAGls38fpLrFPS355W2QQAWgor1KvuYDhTTcyuvIsQXc2YErc2zxVqCcLOK3DT9MOGQHHGiP8mkT?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/aBoX8HKxek0ABApQgJkq0RnSB4adMvYnImh8PUT5v_mB89xp3rPoXyCglp2yHfqD4n85oDuuDSzy_N5DND8vwBhiVczlqXIpux-oY8G3AJbqfkWLkGuJXowguDyhV6IYF4I0J4Jk3i8Et3qrGoVM1ssEtWx3-MzYb865IV_eUklanDWcTkvThjNUuQMTmnmO?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/qMfZXd-4Xzxi8NfuoyP9mY5eZ2yzz2BqMcY1AcIJxDbE-Yvlb7vLLIfjjrNoi9EzYumW_tOUN8XynP0rC_1FVg2_PqbW2k6ZnT2JPY_9bQNRUj3Fv7YQ77ND8tG-te463eZjRfjwmGmwST6-DmV6TeJ703t9qBAkRt1EsiU4ESOs_F5Pb39l4KypwJiLl2qP?purpose=fullsize)

**Sentinel-2** is an optical satellite mission. It records reflected sunlight in multiple wavelengths.

Among the important bands are:

- **Blue**
- **Green**
- **Red**
- **Near-infrared (NIR)**
- Several **shortwave infrared (SWIR)** bands

These allow you to calculate things such as **NDVI**.

---

## The cloud problem

Optical satellites can't see the ground properly through clouds.

Imagine collecting images:

```text
June 1       ☀️
June 6       ☁️☁️
June 11      ☁️
June 16      ☀️
June 21      ☁️
June 26      ☀️
```

You don't want clouds contaminating your vegetation measurements.

So you perform **cloud masking**.

### Cloud masking

The satellite data contains information that allows you to identify pixels likely to be:

☁️ clouds
🌫️ cloud shadows
❄️ snow
🌊 water, depending on the application

You create a mask:

```text
Image

🌿 🌿 ☁️ 🌿
🌿 ☁️ ☁️ 🌿
🌿 🌿 🌿 🌿

       ↓

Mask clouds

🌿 🌿 ❌ 🌿
🌿 ❌ ❌ 🌿
🌿 🌿 🌿 🌿
```

The cloudy pixels are excluded from the analysis.

---

# 4. What is compositing?

Even after cloud masking, you might have several usable images.

Instead of choosing one image, you can combine observations from multiple dates into a **composite**.

For example:

```text
Image 1: 🌿 🌿 ❌ 🌿
Image 2: 🌿 ❌ 🌿 🌿
Image 3: ❌ 🌿 🌿 🌿

              ↓

Composite: 🌿 🌿 🌿 🌿
```

Each pixel is constructed from the best/representative available observation.

A common approach is a **median composite**.

For each pixel:

$$
Composite = median(x_1,x_2,x_3,...)
$$

This reduces the influence of unusual observations.

### Why do this?

If you're comparing vegetation **before vs after a flood**, you want:

**Before composite**

versus

**After composite**

rather than:

**cloudy image before**

versus

**partly cloudy image after**.

That makes the comparison much more reliable.

---

# 5. DEMs and terrain analysis

A **DEM** is a **Digital Elevation Model**.

Think of it as a giant grid where every pixel has an elevation value.

```text
500  510  525  540
480  490  505  530
450  460  475  500
420  430  445  470
```

Those numbers represent elevation above a reference level.

From a DEM you can calculate:

- elevation
- slope
- aspect
- drainage
- river networks
- valley structure
- terrain ruggedness

![Image](https://images.openai.com/static-rsc-4/H-qJTeq29cA7-ocG6mM_iyyaCVHfesK64W-btCiQKoV4YbANcG9VhY4KZ5CgUkYEb24lmAMdhRuP85CzXywRy6bdhWlc0DkPND9omE3CA_GY0y5bU2FfJOAdODy8MIrKFVWz3IQGvJo-lxZi-crcwNuC2OC7HrTQ1WRUSuyoQJ2Vcufq7LtY0zQnImoxL5DO?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/rxFiKvUdvDkjyAK0-OE6rX9FJ_O3uRIKUgmaEDWPwfLwpvXr7_AVQ-6GeSakCW5Hy2fFXk3BU-n2dwoVfPq8VIOT4XRpirCQdmRH4aVmbCFOUI0daueeXHWkMOin-pfWBDLwYE9tVb5vEQNYW8DVpFbJYfDngM64pnzO4Hh0UVtzpVA59SA7ctK2X5tsEzHr?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/RcdmFd6OvPafFyLvf7AWbArrZD7jtaiWDL5zeYoJ3jMAzCqBWZylNJRV8PzdwCXSnbAknMsVFKml3YY2Iqe_DgiR16AGIr_mEr7oOQX818V35FJSp2-NTNOIHyhWNlF-n38kmCNZiG8_xkxwavNEI_g_Ajun7xGI3HsiE-45i5FW8lkMe1PPQwebUIprB6m6?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/cr0TKBHebYjhsTypU-ITFl230KTxxrZ5fAM4AoJ7MEQTTDNR2CrAlVJfWHAL0-l1msKkexhi330OiMgAGUUW4SLfsnYu31qOMuDpiWX3nsB6xVGtWs9bAEPLIjxM0GxtPVt6ji_RFvDqmhOTQkUCKNYs6DVTs48J2KXrFcwaWBGknN65gz3XMyhYDb6BM98M?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/3j_z0-L9aM00Owk0rvDjhpPrnqzjeThjhaI5KgPkTYOO-153dEPSaEYbrGVoM9KkSiF0p9g9B_N_KZGmiukj97WRPWsw2DwgQkoZ5cMGR892WGRZZTKtO_EYvzvNDJRfdZVriftRP4Cg0M3cwADyM5tNaupLRpkZl06Sm6Bj4XgMlgaEhVA4v2IPdKxyVxaq?purpose=fullsize)

---

# 6. HAND — Height Above Nearest Drainage

This is probably the most important terrain concept in the paper.

**HAND = Height Above Nearest Drainage**

Instead of simply asking:

> "How high above sea level is this location?"

HAND asks:

> **"How high is this location compared with the nearest river/drainage channel?"**

That's much more useful for understanding **flood exposure**.

---

### Example

Imagine this cross-section:

```text
       hill
        /\
       /  \
      /    \       🌳
     /      \      ↑
    /        \     │ HAND = 20 m
   /          \    │
__/____________\___│____
       river 🌊
```

The hill might be **500 m above sea level**.

But that's not necessarily useful for knowing whether floodwater can reach it.

HAND says:

> This location is **20 m above the nearby river channel**.

That's much more relevant to flooding.

---

## Why HAND is useful for debris floods

Suppose a debris flood travels down a river.

Locations close to the river and only a few metres above it are potentially much more exposed than locations high on the valley slope.

So you can classify the landscape:

```text
HAND

0–2 m       🔴 very close to channel
2–5 m       🟠
5–10 m      🟡
10–20 m     🟢
>20 m       🔵
```

The paper can therefore ask something like:

> **Does the effect of the debris flood change depending on how high the location is above the nearest drainage channel?**

That is a powerful way to connect **terrain** with **observed environmental change**.

---

# 7. River centerlines and perpendicular transects

This part is much more geometric.

First, you need a **river centerline**.

Imagine the river as a long winding line:

```text
       ~~~~~
      ~
     ~
    ~
     ~~~~
        ~~~
```

The **centerline** represents the middle of the river channel.

Then you create **transects** perpendicular to the river.

![Image](https://images.openai.com/static-rsc-4/8jxOPjn7vRGRYoJ9sPu5phkYo5yAEbVlVswSyoYi4D6tLsYwC3ddYLPSsUDDAlSBWaVQYzv7knqV5VRm7Pta9OSLdLveVFDOZjEYaYBhqNg6ypLqd44O4zwIEgcUFv13efFHYReDoC5bw2D2jGhvwk354Tlnx9w1Yl4h3-zgVhOSdbDVDIABlwI3p4RClNBQ?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/U7rMJWN2Dcylv90HTGLT7oZbhQ9K_Tsd5bv6PQHsXofwlsXc72fXZJkNXq2HScSiaYlCrOZ5tOuEpNhrVadBf66xMLfCU3GAdJg5Kh0jSu4URJcZyK9QZE7C2vANDZ9MSA9b-5S9tyHy8pZbPfM9NnOKL8yjRPlFVTGP_0X06xcXm5W8l8fk_SUJ3kFCNIb3?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/AD0AGU9AAFwaZBao4mbxvW1wuShNfmsHPmpqrecMDHPf-xhhan6Dquc4tZUnlZ49PHzfkGc_pS6B7pDAt04-ijvB4fVLRYAMV_juI5-JkwJkNgcBQlhltj0LElIqtq7zu_DdNg1EO-PaRTTOYuQyTKjpd7noR0Bx5xy7WZ65ovliCFIqW9C-po3w9QI8fkhz?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/2ShV-JdePwWCwkj5H1YEdE6ljvEgPVw3Y-qmsADhZRawlUzZmJyJm0MFMKqTZ24pbPXwiDNuOf34GAfcfnzB1FIbkNpqX8N2ee73fS2LZf-NJnF8EGUvDj6FACxIWFvAxXqZ6wu7-pCXuV91b13VpK9657jGpE_23vnYiHzzKG7fCD7oipxsdCOKyx-rgJaL?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/_dFmXdyA7L9AUWE_cKSJnH56Fxvk4iUTPnQDuLrRfRZHXyEd3OXW7movkeRXXyInsPiwSxw9x7zyRDDyA9YQ3QQQwv0XfHg_vCF9IXz2NZYBBI1Xjo4jDOQKLYRTvvyCoqj4ucYLPfTRIXyNwmCPIeZh431Ei9QarxojCnOnl3MD4lFxRD8nE9-5tZ10eiTP?purpose=fullsize)

![Image](https://images.openai.com/static-rsc-4/wTcKLdAvAnGbQKm9U_xuq1PMTggdVQdIJsq5wyG2mY7NphqpljnKKe7IelSSpdTaGfE4QeNgt8-HD8V728XiAaQ1utHRK7FHemV7DRFy6JgcjGJyYzGZo_ye5CoYT5OQTDoT7UqrlP8FyBmUjwI7zNgM2Ek_w-hufAwdTGXbUEPc5oQSx3PKOlvWRaBPPZ-v?purpose=fullsize)

For example:

```text
             |
             |  Transect
             |
      -------+-------
             |
             |
              \
               \
                \  River centerline
                 \
                  \-------
                    |
                    |
                    |
```

Each transect crosses the river and extends into the surrounding landscape.

---

## Why perpendicular?

Because you're trying to measure changes **across the river valley**, rather than simply moving downstream.

For example:

```text
Valley side                    Valley side
      \                              /
       \                            /
        \        🌳🌳              /
         \       🌱🌱             /
          \______🌊_____________/
                 ↑
              river

       ←──── transect ────→
```

You could measure along each transect:

- NDVI
- HAND
- distance from river
- SAR change
- vegetation damage

Then you have a standardized way of comparing locations.

---

# 8. The programming problem

If you were coding this, the basic workflow could be:

```text
River geometry
      ↓
Create centerline
      ↓
Sample points along centerline
      ↓
Calculate direction of river
      ↓
Rotate direction by 90°
      ↓
Create perpendicular transect
      ↓
Sample raster data along transect
      ↓
Get NDVI / HAND / SAR values
```

The interesting part is that the river is **curved**.

So you can't simply draw all transects at the same angle.

You need to calculate the **local direction of the river** at each point.

Mathematically, if the river direction is a vector:

$$
\mathbf{v}=(x,y)
$$

a perpendicular direction can be:

$$
\mathbf{v}_{\perp}=(-y,x)
$$

That's the geometry behind generating perpendicular transects.

---

# 9. The statistics — the lighter part

Now we get to the question:

> **How do we know the flood actually caused the observed change?**

This is where the statistics come in.

---

## Confounding

A **confounder** is another factor that could explain your result.

Imagine:

**Flooded area → vegetation decreases**

You might conclude:

> "The flood destroyed the vegetation."

But perhaps the flooded area also happened to be:

- at a much lower elevation
- closer to roads
- more heavily farmed
- naturally less vegetated

Those factors could influence vegetation independently of the flood.

That's **confounding**.

```text
Flood exposure ─────→ Vegetation change
       ↑
       │
   elevation
```

You want to separate the flood effect from other differences.

---

# 10. Matching

**Matching** tries to create a fair comparison.

Suppose you have:

### Flood-affected location

- HAND = 4 m
- slope = 10°
- vegetation = forest
- elevation = 300 m

You find an unaffected location with similar characteristics.

Then compare them.

```text
Affected                 Control
   📍                       📍
HAND 4 m                  HAND 4 m
slope 10°                 slope 11°
forest                    forest
```

Now the two locations are more comparable.

The idea is:

> **Make the affected and unaffected areas as similar as possible except for the flood exposure.**

---

# 11. Block bootstrap

This sounds complicated, but the basic idea is manageable.

Ordinary statistical methods often assume observations are independent.

But geographical observations aren't necessarily independent.

Two pixels right next to each other are probably more similar than two pixels 50 km apart.

That's called **spatial dependence**.

So instead of randomly resampling individual pixels, you can resample **blocks**.

```text
┌────┬────┬────┬────┐
│ A  │ A  │ B  │ B  │
├────┼────┼────┼────┤
│ A  │ A  │ B  │ B  │
├────┼────┼────┼────┤
│ C  │ C  │ D  │ D  │
├────┼────┼────┼────┤
│ C  │ C  │ D  │ D  │
└────┴────┴────┴────┘
```

Rather than randomly selecting individual pixels, you resample groups like **A, B, C, D**.

This preserves some of the spatial structure.

---

# 12. E-values

An **E-value** is used in causal inference to assess how strong an **unmeasured confounding factor** would need to be to explain away an observed association.

In simple terms:

> **"How powerful would some unknown factor have to be to completely explain our result?"**

Suppose a study finds a strong relationship between flood exposure and vegetation damage.

An E-value asks:

> Could an unmeasured variable realistically be strong enough to explain this entire relationship?

A **larger E-value** generally means the result would require a stronger hidden confounder to disappear.

It's therefore a way of testing the **robustness of a causal conclusion**.

---

# Putting the whole paper together

This is the most important part to understand.

The concepts aren't isolated. They form a chain:

```text
                    DEBRIS FLOOD
                         │
                         ↓
              ┌────────────────────┐
              │ SAR change detection│
              └─────────┬──────────┘
                        │
                Did the landscape
                    change?
                        │
              ┌─────────┴──────────┐
              ↓                    ↓
         Sentinel-2              SAR
       vegetation/NDVI       radar response
              │                    │
              └─────────┬──────────┘
                        ↓
                 RIVER GEOMETRY
                        │
                 centerline
                        ↓
                  transects
                        │
                        ↓
                  TERRAIN DATA
                        │
                       DEM
                        ↓
                      HAND
                        │
                        ↓
          "How high above the river?"
                        │
                        ↓
                STATISTICAL TEST
                        │
        ┌───────────────┼──────────────┐
        ↓               ↓              ↓
    Matching       Block bootstrap   E-value
        │               │              │
        └───────────────┴──────────────┘
                        ↓
              How confident are we
              that the flood caused
                 the observed change?
```

### The easiest way to remember it

**SAR** → _Did the surface change?_
**Sentinel-2 / NDVI** → _Did vegetation change?_
**DEM** → _What does the terrain look like?_
**HAND** → _How high is this location above the river?_
**Transects** → _How does everything change as we move away from the river?_
**Matching/statistics** → _Can we reasonably attribute the difference to the flood?_

And the particularly important SAR warning is:

> **In steep valleys, a SAR "change" isn't automatically a real ground change — layover and shadow can create misleading signals.**

That is why the **SAR + DEM + terrain analysis** combination is so important in this type of study.
