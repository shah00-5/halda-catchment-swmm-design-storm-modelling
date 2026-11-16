# Event-Based Runoff Modelling of a Halda River Canal Catchment, Chattogram

Design-storm runoff modelling of a 2.23 km² canal catchment beside the Halda River using EPA SWMM 5.2, ArcGIS, Python and NASA IMERG rainfall.

## Overview

Low-lying parts of Chattogram drain through small canals into rivers such as the Halda. This project asks:

1. How much rain falls in 2- to 50-year storms at the site?
2. What peak flow and runoff volume does the catchment produce?
3. How does the canal respond, and how sensitive are the results to storm shape, subcatchment width, infiltration and Halda river level?

The model simulates single design storms. It does not simulate continuous long-term flow.

## Workflow

```
GIS (ArcMap 10.8)  →  Python (rainfall analysis)  →  EPA SWMM 5.2
subcatchments,        Gumbel frequency analysis,     runoff + dynamic-wave
imperviousness, DEM   24-h alternating-block storms  canal routing
```

## Study Area

- East bank of the Halda River, Chattogram
- Catchment area 2.23 km², divided into 5 subcatchments (S1–S5) draining to canal junctions J1–J5
- Canal 2,545 m long, about 35 m wide at the top, with its outlet at 22.5195°N, 91.8495°E
- About 30% impervious (Landsat 8 classification)

## Data and Methods

- **Rainfall:** NASA GPM IMERG V07 Final Run (0.1° pixel), daily 2000–2025 plus half-hourly data for eight extreme events
- **Frequency analysis:** Gumbel distribution fitted to 24 annual daily maxima (2001–2024)
- **Storm shape:** alternating-block 24-hour storms in 30-minute steps for T2, T5, T10, T25 and T50, plus a high-intensity sensitivity case (H)
- **GIS inputs:** Landsat 8 maximum-likelihood classification for imperviousness, SRTM DEM for canal levels
- **SWMM:** Horton infiltration, trapezoidal canal (23 m bottom width, 2:1 side slopes, 3 m depth, assumed), dynamic-wave routing

## Key Results

| Storm | 24-h rain (mm) | Peak at outfall (m³/s) | Runoff volume (million L) |
|---|---|---|---|
| T2 | 154 | 7.10 | 190 |
| T10 | 208 | 11.49 | 296 |
| T50 | 255 | 15.54 | 393 |
| H10 (high intensity) | 208 | 18.16 | 300 |

- The T10 run balances well (continuity error −0.004% for runoff, −0.07% for routing). A rational-method check gives 12.6 m³/s against 11.49 m³/s from the model.
- Impervious land (about 30% of the area) produces about 47% of the T10 runoff volume.
- In the assumed canal section, depth stays below 0.8 m with no flooding. A **Halda stage near 5 m**, not canal capacity, is the main risk: it brings the lower canal to ground level at J1.
- The peak is most sensitive to storm shape and subcatchment width (±20%), and barely to infiltration or river level.

## Repository Structure

```
├── report/    # Full report (.docx)
└── README.md
```

Suggested additions: `python/` (Gumbel fit and storm generation scripts), `swmm/` (`.inp` model file), `gis/` (subcatchment and classification layers), `figures/`.

## Limitations

- IMERG is a satellite estimate on an ~11 km pixel and was not checked against rain-gauge data.
- Canal geometry, subcatchment width, outfall stage and Horton parameters are assumptions, and the model is not calibrated. The results show relative behaviour, not exact flood levels.
- The land-cover classification was checked visually only.

**Next steps:** survey the canal and ground levels, obtain gauge rainfall and Halda stage records, and calibrate against a measured event.

## Tools

EPA SWMM 5.2 · ArcMap 10.8 · Python · NASA IMERG · Landsat 8 · SRTM

## Author

Kaif Al Shahriar Qais Shurid
