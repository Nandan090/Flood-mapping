# Placebo period flood check: melamchi2021

Status: **confirmed by the user on 2026-11-02.** After reviewing the satellite checks and
the high-resolution Wayback imagery, the user confirmed that no notable flood occurred
in the placebo interval, 3 Aug 2018 – 15 Dec 2019 (spec 14.2). This confirmation covers
*notable* floods. Section 7 lists what the checks cannot rule out, including local floods
and bank erosion.

Who reviewed what:

- *The placebo pair (Oct–Dec 2018 → Oct–Dec 2019):* reviewed by the user.
- *The 2018-monsoon check (Mar–May 2018 → Oct–Dec 2018):*
  - Claude, in this session.
  - Claude, independently in a separate session, from the PNG chips and sheets only.
  - The user, who reviewed the results before confirming.
- *Wayback imagery at four patches:* reviewed by Claude and then by the user. Patch 26
  remains unchecked at high resolution.
Nothing in this check changes the impact labels, the placebo method or the
real-event forest baseline (Sentinel-2, 1 Oct–15 Dec 2020).

## 1. The placebo period (spec 14.2, `config.yaml`)

| Input | Window |
|---|---|
| Sentinel-2 placebo "pre" (forest, channel, controls) | 2018-10-01 to 2018-12-15 |
| Sentinel-2 placebo "post" | 2019-10-01 to 2019-12-15 |
| Sentinel-1 placebo "pre" (baseline) | 2018-08-03 to 2018-09-15 (end exclusive) |
| Sentinel-1 placebo "post" | 2019-08-03 to 2019-09-15 (end exclusive) |
| Fake event date | 2019-06-15 |

"October 2018 to October 2019" is shorthand. The optical composites cover
October–mid-December of each year, and the radar baseline starts on 3 August 2018.
The full span that must be free of floods is therefore **3 August 2018 to
15 December 2019**. A flood inside it would appear in the placebo as change, or
would sit inside the placebo baseline.

## 2. Desk review (user)

- The sources reviewed show no major Melamchi River debris flood in 2019 comparable
  to June 2021.
- The Nepal Disaster Report 2019 mainly covers 2017–2018, so it cannot show on its
  own that 2019 was flood-free.
- Local flooding, erosion and landslides in the wider Melamchi–Helambu area cannot be
  ruled out from the reports. That is why the satellite checks below were run.

## 3. Method

**Imagery.** Sentinel-2 composites built the same way as the study's: Cloud Score+
`cs_cdf ≥ 0.60`, per-band median, and at least 2 clear observations per pixel, all on the
10 m study grid:

- Mar–May 2018
- Oct–Dec 2018
- Mar–May 2019
- Oct–Dec 2019

Sentinel-1 per-orbit median VV and VH (dB) for:

- Mar–May 2018
- Aug–Sep 2018
- Oct–Nov 2018
- Aug–Sep 2019
- Oct–Nov 2019

Every scene over the corridor from August 2018 to December 2019 is listed in
`placebo_check/scene_list.csv`, with the share of the corridor that was clear. The scenes
used in each composite are in `placebo_check/manifest*.json`.

**Screen.** For each pair of composites, the study's optical label rules (5.2, 14.18,
14.19 a, connectivity filter 5.4) were applied with the earlier composite as "pre"
and the later one as "post". The result marks new bare ground on ground that was
vegetated at the start. No radar label was computed; the radar images were only
inspected by eye.

**Patches listed.**

- *Which patches:* patches of at least 500 m² that are either connected to the channel
  or within 300 m of the river line.
- *Patch numbers:* patches are numbered by connection first, then by size.
- *Recorded for each patch:* river km, bank, area, distance to the line, coordinates, and
  height above nearby channel pixels. Height is computed from FABDEM, which is 30 m, so it
  is rough.

**Screening classes.** These are screening choices, not study thresholds. They were set
after a first look at a few chips.

- *Hillslope (landslide or gully):* more than 20 m above the nearby channel.
- *Channel margin (flood-compatible):* connected to the channel, within 30 m of the river
  line, and within 10 m of channel height.
- *Floodplain or terrace field / building:* everything else.

"Flood-compatible" means only that a flood could have produced the change.
Ordinary monsoon high water, bar reworking and bank erosion look the same at 10 m.

## 4. Check 1: Oct–Dec 2018 → Oct–Dec 2019 (the placebo pair)

**Screen results.**

- *Connected to the channel:* 13 patches, 1.47 ha in total over 35.5 km. The largest is
  0.26 ha.
- *Not connected, within 300 m of the line:* 30 patches, 3.26 ha.
- *Active channel:* 193.5 ha before and 155.0 ha after, which is 20% narrower. It did not
  widen in any river km.
- *Not assessable:* 4.5% of corridor pixels are unknown or dark (already bare before, or
  shadowed).

| Class | Patches | Area (ha) |
|---|---:|---:|
| Channel margin (flood-compatible) | 5 | 0.86 |
| Floodplain or terrace field / building | 22 | 2.51 |
| Hillslope (landslide or gully) | 16 | 1.36 |

**Channel-margin patches.** The table lists the five channel-margin patches, plus patch 3,
which the user reviewed and which the rules place in the field class.

| Patch | km | Bank | Area (m²) | Bare by May 2019 | Review |
|---:|---:|---|---:|---:|---|
| 1 | 0.29 | right | 2,600 | 58% | Claude: edge of a bare fan that was already there; change mostly before the monsoon |
| 2 | 16.89 | left | 2,300 | 96% | Claude: fields bare in the dry season; not a deposit |
| 3 | 17.37 | left | 2,100 | 0% | User: no clear flood signal |
| 4 | 5.51 | right | 1,600 | 12% | User: minor ambiguity, no convincing flood deposit |
| 5 | 20.90 | right | 1,100 | 9% | User: no clear flood signal; narrow channel-edge change |
| 6 | 20.74 | right | 1,000 | 30% | User: minor ambiguity, no convincing flood deposit |

**User's visual review.**

- No convincing evidence of a widespread flood-related disturbance.
- The overview sheets show no continuous new sediment belt, no widespread bank stripping
  and no major channel widening.
- The Sentinel-1 comparisons show no obvious coherent new surface-change signature.

## 5. Check 2: Mar–May 2018 → Oct–Dec 2018 (the 2018 monsoon)

This check was added to cover the 2018 monsoon. That stretch falls inside the radar
baseline (Aug–Sep 2018) and before the optical "pre" composite.

**It compares different seasons, which the study avoids (4.1).** Fields bare after the
autumn harvest look like new bare ground. Spring NDVI is also lower, which changes the
"unknown" mask. This check therefore over-flags cropland by design, and its patch counts
are not comparable with Check 1.

**Screen results.**

- *Connected to the channel:* 67 patches, 14.21 ha.
- *Not connected, within 300 m of the line:* 156 patches, 37.89 ha.
- *Active channel:* 262.7 ha in spring and 193.5 ha in autumn. It widened only in km 7,
  km 18 and km 25, each by 0.1 ha or less.
- *Not assessable:* 12.1% of corridor pixels are unknown or dark.

| Class | Patches | Area (ha) |
|---|---:|---:|
| Channel margin (flood-compatible) | 15 | 3.36 |
| Floodplain or terrace field / building | 118 | 20.21 |
| Hillslope (landslide or gully) | 90 | 28.53 |

**Claude's visual review.** All 15 channel-margin patches were reviewed, plus patches 21
and 26.

| Group | Patches | What the chips show |
|---|---|---|
| Flood-compatible, minor channel-edge or bar change | 15 (km 11.7 R), 16 (km 28.0 R), 26 (km 3.3 L), 36 (km 29.7 R) | Narrow strips or a single gravel bar along the active channel turned brighter after the monsoon. Each is at most 0.27 ha. This fits ordinary monsoon high water or bar reworking. No sediment spreads onto the banks. |
| Agricultural | 3, 8, 9, 10, 34, 41, 45, 47, 64 (km 20–33) | Blocks on terraces and valley-floor fields, bare after harvest. Several lie near the 2021-traced river line but away from the 2018 channel (see limitations). |
| Unclear: field, gravel extraction or construction | 29, 31, 59 (km 33.7–34.4 L, near Melamchi Bazaar and the confluence) | Plots at the edge of the wide gravel plain turned bright. The gravel plain itself did not spread. |
| Hillslope | 21 (km 1.5 R) | A gully on the hillside, about 160 m above the river, extended up its head. Not a river flood. |

**Overview sheets.** Through the gorge (km 0–12) only a few short channel-edge strips
appear. From km 15 down, the red patches are scattered blocks on terraces, not a
continuous belt along the river. Hillslope scars appear at km 8 and km 12.6, and are
not connected to the channel.

**Independent review.** A separate Claude session reviewed the same chips and sheets.
It saw only the images, not the tables or counts. Its verdict:
**no convincing evidence of significant flood disturbance.** It found:

- no connected sediment belt
- no fans or bank stripping
- small patches, mostly on farmland or at channel margins
- no coherent radar change at any outline

| Patch | This session | Independent review (confidence) |
|---:|---|---|
| 15 | Minor channel-edge change | Ordinary seasonal river change (medium) |
| 16 | Minor channel-edge or bar change | Ordinary seasonal river change (low); the bank part could be small bank erosion |
| 26 | Minor channel-edge change | Uncertain (low); gorge, radar distorted by layover and shadow |
| 36 | Minor channel-edge change | Agricultural (medium) |
| 29 | Unclear: field, gravel extraction or construction | Ordinary seasonal change at the gravel-plain margin (low) |
| 31, 59 | Unclear: field, gravel extraction or construction | Uncertain (low); farmland or human-modified ground |
| 21 | Hillslope gully | Landslide or gully (medium) |
| 3 | Agricultural | Agricultural (medium) |

**Where the reviews agree and differ.**

- Neither review found a flood-type deposit. The "flood-compatible" class in this report
  (section 3) covers the same changes the independent review calls ordinary seasonal river
  change: gravel exposure and water level at the channel margin.
- The two reviews differ on patch 36 (channel edge or field) and patch 29 (field or
  gravel-plain margin). Neither reading points to a flood.
- Both reviews left 31 and 59, near Melamchi Bazaar, unresolved at 10 m. The high-resolution
  check below resolves them as a gravel-processing or construction yard and fields.

**High-resolution check (Esri World Imagery Wayback, done 2026-11-02).** The capture
date was looked up at each patch in all 198 Wayback releases, using the metadata layer at
map level 17. The two dates the independent review suggested (2017-02-27 and 2019-05-15)
are Wayback *release* dates. The captures they show at these patches are 1 Feb 2016 and
1 Nov 2018.

Images used: 0.5 m captures from **1 Feb 2016**, **1 Nov 2018**, **7 May 2019** and
**3–19 Apr 2020**, all DigitalGlobe or Maxar.

- *2016 → Nov 2018* spans the 2016, 2017 and 2018 monsoons, so it cannot isolate 2018 on
  its own.
- *Nov 2018 → May 2019 → Apr 2020* brackets the 2019 monsoon. That is inside the placebo
  interval and before the 2021 flood.

| Patch | km, bank | What the 0.5 m images show | Class (confidence) |
|---:|---|---|---|
| 31 | 33.82 L | A levelled pad becomes an industrial yard with machinery and a plant (May 2019) and stockpiles (Apr 2020). It grows steadily from 2016 to 2020. | Construction or gravel processing (high) |
| 59 | 33.69 L | Terraced fields on every date. Boundaries are intact and plots alternate between green and harvested. A crusher yard is growing just to the south. | Agricultural (high) |
| 29 | 34.36 L | Scrub crossed by a small sandy side channel or tributary and dirt tracks. The tracks expand in May 2019 and the area revegetates by Apr 2020. No sediment sheet spreads over the scrub. | Small side channel and tracks; a minor local sand deposit in that channel cannot be excluded (medium) |
| 16 | 27.98 R | A smooth tan pad between the road and the river appears by Nov 2018, replacing 2016 vegetation. It partly revegetates by Apr 2020. It is flat, sharp-edged and joined to the road, and confined between the road and the channel. | Uncertain: roadside clearing or stockpile pad, or sand laid on the bar (low) |
| 26 | 3.29 L | No Wayback capture between Dec 2015 and Dec 2021 | Not checkable in Wayback |

At patches 16, 29, 31 and 59, no new deposit appears between Nov 2018 and Apr 2020, which
covers the 2019 monsoon. The images are kept in `outputs/melamchi2021/placebo_check/wayback/`
(`compare_patchNN.png` and `zoom_patchNN.png`). That folder is git-ignored, and the images
are Esri content.

**Point for the main study (from the independent review).** In this spring-to-autumn
comparison, the false "new bare" patches are concentrated on farmland in km 12–35. They are
nearly absent from the forested upper gorge. If cropland change of this kind leaks into the
real flood label, it adds false impact to non-forest banks. That biases the study toward
finding protection (spec 14.17). This comparison mixes seasons, so some of the effect is
expected. The same-season placebo (14.2) is the proper test of how much leaks.

**Status: confirmed by the user on 2026-11-02, together with the whole placebo interval.**

## 6. Conclusion (confirmed by the user, 2026-11-02)

The user confirmed that no notable flood occurred between 3 Aug 2018 and 15 Dec 2019.
The evidence behind that confirmation:


- In both checks, the river-connected change is small, scattered and mostly
  agricultural.
- The flood-compatible patches are narrow channel-edge or bar changes of at most 0.27 ha
  each.
- No continuous new sediment belt, widespread bank stripping or channel widening was
  found between March 2018 and December 2019.
- This supports using 2018–2019 as the placebo period.
- It does **not** show that no flood occurred. Local floods, high flows and bank erosion
  below the limits in section 7 may have gone undetected.

## 7. Limitations

- **Resolution.** Sentinel-2 is 10 m. Bank erosion or deposits narrower than about one
  or two pixels can be missed.
- **Composites.** Each period is a seasonal median. A short-lived flood mark that
  revegetated, or was cleared, before the next composite would not show.
- **Single clear monsoon images were not used.** Clear images exist for 2019-06-07,
  2019-07-02, 2019-08-06, 2019-08-26 and 2019-09-20, and for 2018-09-18.
- **Check 2 compares spring with autumn.** Harvest, crop and greenness differences are
  mixed into it.
- **River line.** The river line was traced on post-2021 imagery. In the lower valley it
  departs from the 2018 channel, so "distance to line" is approximate and some field
  patches fall in the channel-margin class.
- **Heights** come from FABDEM, which is 30 m. They are rough, especially in the narrow
  gorge.
- **Screening rules.** The classes are screening rules set after a first look at a few
  chips. Claude's visual classifications are judgement calls, and the 2018 monsoon ones
  still need the user's confirmation. Both Claude reviews used 10 m images. The
  independent one saw only the PNGs.
- **Radar** was only inspected by eye, in the chips (orbit 121). No radar change label was
  computed.
- **High resolution.** Wayback imagery was checked only at patches 16, 29, 31 and 59. Its
  captures are single dates, and there is none between Feb 2016 and Nov 2018 at these
  points. Patch 26 has no capture between Dec 2015 and Dec 2021.
- **Not assessed.** Pixels that are unknown or dark could not be checked: 4.5% of the
  corridor in Check 1 and 12.1% in Check 2.
- **Desk review** found no documented 2019 flood, but the main report used covers mainly
  2017–2018.

## 8. Files (in `results/placebo_check/`)

- **Both checks** (`2018_2019/` and `monsoon2018/`):
  - `change_patches.csv`: every listed patch, with its class, coordinates and a Google
    Maps link
  - `channel_change_per_km.csv`
  - `screen_summary.json`
  - `sheet_new_bare_*.png`: the whole river in 2-km panels, new bare ground in red
  - `chips/`: the patches discussed above
- **Scene records:** `scene_list.csv` and `manifest*.json`.

The full set of chips and composites is in `outputs/melamchi2021/placebo_check/`,
which git ignores.
