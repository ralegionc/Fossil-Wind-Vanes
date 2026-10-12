# Analysis plan (written 2026-10-12, before any comparison of crests with predicted winds)

## Question

Do giant dunes still point along the winds of the Last Glacial Maximum, while the small dunes on top of them have turned to face today's winds?

## Units

Gunn et al. (2022) tiles: 2,093 squares of 0.28125 deg (about 31 km) inside 38 dune fields. Fields are the clusters for all uncertainty estimates.

## Measurements

- **Observed crests.** Directional power spectra of the Copernicus 30 m DEM. Large dunes: band around Gunn's large wavelength. Short dunes: everything below 0.6 times the smaller of Gunn's wavelengths, down to 120 m. Oblique short sets: spectral peaks in the short band at least 15 deg from the large crest, with prominence of at least 0.5 times the profile median. The 15 deg wedge removes harmonics of the large dunes.
- **Predicted crests.** Sand flux from ERA5 10 m winds (6,000 random hours, 1991 to 2020) with a threshold friction velocity of 0.25 m/s. Two modes from Courrech du Pont et al. (2014) with speed-up 1.6: bed instability (BI) and fingering (F).
- **Past climates.** Delta method: each model's monthly mean change in surface wind (LGM or mid-Holocene minus its pre-industrial control) is added to each ERA5 hour of that month. Ensemble = multi-model mean change. Ten LGM models.

## Filters (fixed in advance)

- Large crest must be clear: peak-to-median contrast at least 2.
- Star dunes are dropped (no single crest direction).
- Tiles whose modern prediction is unstable (bootstrap sd over 10 deg) are dropped.

## Mode for each scale (primary)

- Large linear dunes: F. Large barchanoid, transverse and unlabelled dunes: BI.
- Short oblique sets: BI.
- Robustness: best of the two modes, applied identically to every epoch, so no epoch gains extra freedom.

## Tests

- **H1, giant dunes remember the LGM.** Misfit = axial angle between observed and predicted crest (0 to 90 deg). Statistic: mean over tiles of (misfit to modern minus misfit to LGM). Positive means LGM fits better. 95% CI from a bootstrap over fields. Reported for all tiles and for informative tiles, where the LGM and modern predictions differ by more than 10 deg.
- **H2, short dunes follow today's winds.** The same statistic for oblique short sets. The prediction is that it is zero or negative.
- **H3, the gap between scales is the climate change.** Across tiles with an oblique short set, correlate the signed angle (short crest minus large crest) with the signed predicted rotation (modern prediction minus LGM prediction). Field-bootstrap CI.
- **Placebo.** Give each field the predicted LGM rotation of a different, randomly chosen field: each tile's LGM crest becomes its modern prediction turned by the other field's median predicted rotation (1,000 permutations). This keeps the size and spread of real predicted rotations but breaks their link to place. The H1 statistic from the real assignment is compared with this null distribution.
- **Mid-Holocene** is reported alongside as an intermediate climate, with no directional prediction.

## What would count as failure

- H1 statistic with a CI that includes zero, or one that is no better than the placebo, means giant dunes do not carry a detectable LGM signal in this test.
- If modern winds fit the short oblique sets no better than chance (mean misfit near 45 deg), the method cannot see the modern signal either, and the H1 result is uninformative rather than negative.

## Addendum 1 (2026-10-12, written after the primary null and the power check, before any of the tests below were run)

The primary test was underpowered: the monthly-mean delta method turns predicted crests by a median of only 3.5 deg. Two follow-ups, fixed here before running them.

### A. Daily model winds (storminess, not just mean shifts)

- Models with daily `uas`, `vas` for both LGM and pre-industrial: AWI-ESM-1-1-LR and MPI-ESM1-2-LR. Last 30 simulated years of each run.
- At the model grid cell nearest each tile, build flux roses from daily-mean 10 m winds for the LGM and pre-industrial runs. Daily means are weaker than hourly peaks, so the threshold is lowered: primary 4 m/s, sensitivity 3 and 5 m/s.
- Predicted rotation = crest from the LGM rose minus crest from the pre-industrial rose (signed, axial), per mode. The LGM prediction is the modern ERA5 prediction turned by that rotation. Ensemble = mean of the two models' rotations.
- Same H1 statistic, filters, mode rule, field bootstrap and informative subset (predictions differ by more than 10 deg) as the primary test.

### B. Tiles where the models agree on a large change

- Selection uses model output only, never observed crests.
- Monthly ensemble (10 models): tiles where the median signed rotation (type mode) is at least 10 deg in size and at least 8 of 10 models rotate the same way.
- Daily ensemble (2 models): tiles where both models rotate the same way by at least 10 deg.
- H1 statistic with field bootstrap on each subset. If a subset has fewer than 5 fields, it is reported but not interpreted.

### Failure and success

- Success for the memory hypothesis: H1 gain positive with a field-bootstrap CI above zero in A (primary threshold) or in B, and not reversed in the other.
- A gain with CI below zero in these high-signal subsets would count as evidence that giant crests follow modern rather than LGM winds.

### Deviation (2026-10-12, before any MPI-ESM1-2-LR daily result was computed)

DKRZ, the only server holding MPI-ESM1-2-LR daily LGM winds, slowed to about 40 kB/s. That run uses its last 20 simulated years (one file) instead of 30, with `uas` trimmed to the same 20 years. Its pre-industrial run and both AWI-ESM-1-1-LR runs keep 30 years.

## Addendum 2 (2026-10-12, written before either test below was run)

The climate-model tests are limited by disagreement between models. Two tests that need no ice-age model.

### C. Does the misfit to today's winds grow with turnover time?

- Turnover time per tile: T = H x wavelength / (2 x DP), in years. H and wavelength are Gunn's large-dune height and wavelength; DP is the gross sand flux from the ERA5 hours, with the same flux law as the predictions, in m2/yr.
- Outcome: angular misfit of the giant crest to the modern prediction (type mode), same filters as before.
- Model: misfit = field effect + b x log10(T) + c x (RDP/DP) + error. 95% CI for b from a field bootstrap. Memory predicts b > 0 (slower dunes sit further from today's winds).
- Sensitivity: T from this project's own DEM relief and peak wavelength; best-of-two modes; no directionality term; linear dunes only.
- Known-answer check: the fifth of tiles with the shortest T must fit today's winds with mean misfit below 35 deg. If not, the test is not interpreted.

### D. Young dunes from Sentinel-2

- Base tiles. Central 10 km window of the Sentinel-2 L2A red band (10 m). Two scenes per tile with under 5% cloud, one from June to July and one from December to January (2021 to 2024), the clearest of each, so the sun comes from different directions.
- Spectra as for the DEM, rotated from UTM grid north to true north. The direction profiles of the two scenes are averaged.
- Young-dune crest: main crest in the 40 to 150 m band.
- Check 1 (illumination): the Sentinel-2 crest in the giant-dune band must match the DEM giant crest with a median absolute difference below 10 deg, or Sentinel-2 directions are not used.
- Check 2 (young dunes follow today's winds): mean misfit of young crests to the modern bed-instability prediction must be below 35 deg, or the young-dune arm is declared uninformative again.
- If both checks pass: (i) paired within tiles, is the giant crest further from today's winds than the young crest? Best-of-two modes for both scales so neither gets more freedom; field bootstrap. Memory predicts yes. (ii) Young dunes are added to test C as the fastest scale.

### Amendment to D (2026-10-12, after a two-tile pilot of the imagery and before any comparison with winds)

In the pilot, the 40 to 150 m band of the Simpson tile peaked exactly on the giant crest (156.5 deg in both): sharp crest lines put energy at short wavelengths in the giant dunes' own direction, the same harmonic trap as in the DEM. The young-dune crest is therefore the strongest peak in the 40 to 150 m band at least 15 deg from the giant crest, with prominence at least 0.5 times the profile median, as for the DEM short sets. The unmasked main crest is reported only as a secondary result, flagged as contaminated by giant-crest edges. The share of tiles with an oblique young set is reported.

## Addendum 3 (2026-10-12, written before any of the tests below were run)

### E. Reorientation time for long dunes

Test C used the time to rebuild a dune's cross-section. Turning a dune means moving sand along its whole length: rotating a dune of length L by an angle a moves about H x L^2 x a / 4 of sand, supplied at a rate DP x L, so the time to turn it by about one radian is T_rot = H x L / (4 x DP).

- L is Gunn's large-dune width (the long axis of the planform, capped near the 31 km tile size). For short dunes T_rot is close to the cross-section time; for long linear dunes it can be tens of times larger.
- Same model, filters, covariate and field bootstrap as test C, with log10(T_rot) in place of log10(T). Memory predicts a positive slope.
- Also reported: the share of tiles with T_rot above 21,000 years, and the mean misfit difference between tiles above and below 21,000 years within fields that have both (field bootstrap). Memory predicts that the slow tiles fit today's winds worse.
- Known-answer check: as in C, the fastest fifth must fit today's winds below 35 deg.

### F. Young dunes in winds bent by the giant dunes

Giant dunes steer the wind near the ground toward their crest line, most strongly when the wind meets the crest at 45 deg and not at all when it blows along or straight across it. Each ERA5 hour is turned by delta = -k x sin(2 x (wind direction - giant crest)) before the flux rose for the young dunes is built.

- k is chosen on a grid (0 to 45 deg, step 5) to minimise the young dunes' bed-instability misfit, using half the dune fields, and evaluated on the other half. The halves are then swapped (two-fold cross-fitting, fields assigned at random with seed 0), and only held-out tiles are scored.
- Check: the held-out young misfit to the bent-wind bed-instability prediction must fall below 35 deg. If it does, run test D(i) on held-out tiles: giant crest against regional winds and young crest against bent winds, best-of-two modes for both, field bootstrap. Memory predicts that the giant crest misfit is larger. If the check fails, the young-dune arm stays uninformative.
- Physical sanity: a chosen k of 0 means no bending helps; a chosen k at the top of the grid (45 deg) is reported as implausible.

### G. Weighting the models by skill

- Skill: D_m = root-mean-square difference between a model's pre-industrial monthly mean wind vectors (uas, vas) and ERA5's (from the sampled hours), over all base tiles and 12 months.
- Spatial null: the same with each model's tile values shuffled between tiles. A model whose D_m is not below its null's 5th percentile (200 shuffles) gets zero weight.
- Weights: w_m = exp(-(D_m / s)^2) with s the median D_m over the models, renormalised. Skill-weighted LGM change = weighted mean of model changes.
- Tests with the weighted prediction: the H1 statistic (all tiles and informative), the memory-fraction slope, and the turn compared with each model's.
- Across the ten models: Spearman correlation between skill (D_m) and each model's memory-fraction slope and H1 gain. Reported with n = 10 and not over-read.
