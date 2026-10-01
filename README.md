# CRISPRi HLC Trial3 D28 (fixed) — B2M knockdown analysis

**Revision note:** this version adds compensation (spillover matrix from
single-stain controls), an FSC-H:FSC-A singlet gate, an SSC-A upper cap on
the cell gate, a time/acquisition-stability QC check, MFI reporting, a
per-well gating-hierarchy table, numbered plots with Why/How/How-to-read/
Result text throughout, and splits Arm 1 into dox-off/dox-on groups (see
"Doxycycline induction correction" below) — all of which change the Arm 1
knockdown numbers from the previous version of this notebook. Every chart
and label throughout uses reader-friendly names instead of well IDs, and
no replicate wells that differ by a real experimental condition (dox
induction, transduction reagent) are pooled together without being
labeled as distinct groups.

Analysis of `gs://20240617-landerlab-crispri-project/Flow_Data_CRISPRi/20260929_HLC_Trial3_D28_fixed`
(fixed panel only), the sibling dataset to
[crispri-hlc-trial3-d27-fixed-kd](https://github.com/ageorge-sketch/crispri-hlc-trial3-d27-fixed-kd)
one day later. Unlike D27 (CD81 knockdown, 6 arms), this dataset targets
**B2M** across 3 arms comparing different guide/effector delivery modes.
Raw FCS files are parsed directly (no FlowJo CSV export); debris/doublets are
excluded on FSC-A/FSC-Width; knockdown is the percentage-point shift in "%
cells below a fixed gate" on the raw B2M antibody channel between the
non-targeting guide (NT-GFP) and the B2M-targeting guide (B2M-GFP), within the
guide/effector-delivery-gated population.

## Corrected sample sheet

The sample sheet originally uploaded to the GCS folder contained data-entry
errors that made this look like a CD81 arm sharing D27's exact two-virus
design: a stray `AA173 BFP CD81`/`AA173 BFP ORK` guide-virus assignment on
wells D7–D10, a resulting BV421-channel collision between that BFP guide
marker and the in-well ASGR1-BV421 antibody stain, and stray Doxycycline
Induction flags on two hepatocyte-marker wells. The user re-uploaded a
corrected sample sheet to the same GCS path mid-session; this notebook was
rebuilt against the corrected file (re-listed the bucket, confirmed via the
GCS object's `updated` timestamp that a newer file was in place, and re-loaded
fresh rather than reusing any cached copy). After correction: no CD81 guide
appears anywhere in this dataset (all 3 arms target B2M), and there is no
BV421 collision (ASGR1-BV421 is the only BV421 stain in any arm well).

## Doxycycline induction correction

The sample sheet's `Doxycycline Induction` column reads "No" for every well
in Arm 1 (B1-B6), including B5 and B6. B5 and B6 in fact received
doxycycline induction and B3/B4 did not — the sample sheet does not capture
this distinction. WTC11 carries a doxycycline-inducible CRISPRi effector, so
this is Arm 1's main experimental variable: the `B2M-GFP (guide)` wells are
split into **no dox** (B3/B4, effector off) and **+dox** (B5/B6, effector
on) groups rather than averaged together, since pooling them would dilute
the dox-dependent knockdown signal. The arm's headline knockdown number is
the +dox group vs. control.

## Exclusions

- **K562 wells** (`02-Well-A5/A6/A7`) — included only to confirm ASGR1
  antibody specificity, not part of the knockdown analysis (same precedent as
  D27).
- **D5/D6** (guide-only, no-effector singleton wells, n=1 each) — GFP
  compensation controls only, not used for replicate-based statistics.

## Arms

All 3 arms compare `NT-GFP` (non-targeting guide, control) vs. `B2M-GFP`
(B2M-targeting guide), read out on a B2M antibody channel, gated on
guide/effector-delivery markers (channel assignments verified empirically
against single-stain control wells, not assumed from naming conventions):

| Arm | Cell line | Design | Readout channel | n (mock/guide) | Knockdown, 1st pctile (pp) | Knockdown, 50th pctile (pp) |
|---|---|---|---|---|---|---|
| 1, no dox | WTC11 | Transient `B2M-GFP`/`NT-GFP` virus only, no separate effector virus, dox-inducible effector OFF | B2M-APC → `APC-A` | 2 / 2 | 0.2 | 4.9 |
| 1, +dox | WTC11 | Same as above, dox-inducible effector ON | B2M-APC → `APC-A` | 2 / 2 | 5.2 | 6.0 |
| 2 | 17_3 | Stably-integrated guide (`NT-GFP`/`B2M-GFP`) + transient `AA239 dCas9-KRAB Thy1.1` effector virus | B2M-PE-Dazzle-594 → `B610-ECD-A` | 6 / 6 | 16.2 | 28.6 |
| 3 | 17_3 | Fully-transient two-virus: guide virus (`NT-GFP`/`B2M-GFP`) + `AA239` effector virus, no integration | B2M-PE-Dazzle-594 → `B610-ECD-A` | 2 / 2 | 6.6 | 24.8 |

Arm 1's +dox knockdown (5.2pp / 6.0pp) is larger than its no-dox knockdown
(0.2pp / 4.9pp) at the 1st-percentile metric, consistent with a
dox-dependent CRISPRi effect, though both are small relative to arms 2 and 3.

Both metrics per the standing convention: **1st-percentile-of-mock** (red
dashed line, tail effect) and **50th-percentile-of-mock** (purple dashed
line, population-level shift), reported as the percentage-point shift of the
guide group vs. its own arm's mock group.

Thy1.1-APC effector-delivery gate (arms 2/3) uses well **D11** (fully
triple-stained: B2M-PE-Dazzle/Thy1.1-APC/ASGR1-BV421, no guide, no virus) as
an FMO-style stained-but-uninfected reference, per the FMO-gating convention.

Arm 3 has only n=2 replicates per group (flagged in-notebook) — less robust
than arm 2's n=6. Arm 1 (WTC11) shows markedly weaker knockdown than the
17_3 arms; no FMO-style reference was available for arm 1 (WTC11 has no
effector-only compensation well analogous to D11), so it falls back to a
fully-unstained naive reference, same fallback pattern as D27.

## ASGR1-BV421 stratification

Unlike D27 (which required separate Albumin/ASGR1 co-stain wells), ASGR1-BV421
is co-stained directly in every arm's own wells here, so stratification needs
no extra reference wells. ASGR1+ knockdown stays close to the unstratified
numbers in arms 2 and 3; Arm 1's ASGR1+ stratum has too few cells per
replicate well to read reliably (see table and caveats below):

| Arm | Stratum | 1st pctile (pp) | 50th pctile (pp) |
|---|---|---|---|
| 1 (WTC11), no dox | All cells | 0.2 | 4.9 |
| 1 (WTC11), no dox | ASGR1+ | both replicate wells below the 300-cell floor, no value computed | both replicate wells below the 300-cell floor, no value computed |
| 1 (WTC11), +dox | All cells | 5.2 | 6.0 |
| 1 (WTC11), +dox | ASGR1+ | 1.1 (one of two replicate wells below the 300-cell floor; low-confidence) | -12.2 (one of two replicate wells below the 300-cell floor; low-confidence) |
| 2 (17_3, integrated) | All cells | 16.2 | 28.6 |
| 2 (17_3, integrated) | ASGR1+ | 15.9 | 28.8 |
| 3 (17_3, transient) | All cells | 6.6 | 24.8 |
| 3 (17_3, transient) | ASGR1+ | 7.3 | 27.5 |

This dataset does not contain the previously-flagged non-specific ASGR1-FITC
clone at all (ASGR1 here is stained with BV421 or PE only). Albumin-AF647 and
ASGR1-PE exist only on a separate set of hepatocyte co-stain reference wells
(never co-stained together with the knockdown-arm readout antibodies), so
they're shown as an **informational-only** calibration section and are not
used to stratify any arm.

## Arm 2 delivery-reagent comparison

Arm 2's 12 wells (C1–C8, D1–D4) also compare 3 lentiviral delivery reagents
for the `AA239` effector-virus transduction, per the sample sheet's
`Transduction Reagent` column — not surfaced until specifically requested:
abbreviated in every chart/label as **PB**, **PS**, and **LB** respectively,
4 wells each (2 NT-GFP + 2 B2M-GFP per reagent). A bar chart of %
Thy1.1-APC+ (effector-infected) per well sits at the top of Arm 2's tab in
the infection-gating section (separate from the knockdown-metric section
below it) and is fully reactive to the same FMO-calibrated Thy1.1-APC gate
slider used for that arm's knockdown gating. Each bar is labeled by reagent,
guide condition, and replicate number (e.g. "PB: NT-GFP rep1") and grouped
into 3 reagent blocks separated by dotted lines — no two wells from
different reagents ever share a label, since "rep1" under one reagent is a
different physical well than "rep1" under another, not a matched replicate.
Arm 2's distribution and scatter knockdown plots also split each guide
condition by reagent for the same reason, rather than pooling all 3 reagents
into one trace. At the baseline gate, all 3 reagents land in a similar
~75–81% infected range (Polybrene/Protamine Sulfate slightly ahead of
Lentiboost, ~78–80% vs ~75–76%) — no dramatic reagent-driven confound.

## Debris/doublet QC scatter fix

The pooled debris (FSC-A vs SSC-A) and doublet (FSC-A vs FSC-Width) scatter
plots initially rendered as empty-looking. Root cause: the pooled *histogram*
QC cells always clipped their arrays to a sane range before plotting, but the
*scatter* cells never applied that same clip. ~12% of SSC-A events in this
dataset reach up to ~9.16e7 (vs. a real population sitting in 0–3e6), so
Plotly's autoranged axes stretched to fit those extreme outliers and
compressed the real ~88% of valid events into an imperceptible sliver near
the origin — the plots weren't empty, they were rendered at the wrong scale.
Fixed by clipping both scatter cells to the same ranges as their companion
histograms (FSC-A/SSC-A to `(-50_000, 3_000_000)`, FSC-Width to `(0, 8000)`)
and pinning explicit axis ranges so a future dataset's outliers can't blow
the scale out again. Verified with actual point counts (112,500 points in
each scatter, matching the pre-fix well/cap count, confirming no wells were
silently dropped by the new clip) rather than just checking for cell errors.

## Notebook

`notebook.py` — PEP 723 sandbox header with pinned deps (numpy, pandas,
plotly, requests, xlrd). Open with:

```bash
uvx marimo edit --sandbox notebook.py
```

GCS access uses the same bearer-token pattern as the D27 pipeline: paste the
output of `gcloud auth print-access-token` into the password-masked
`token_input` widget (expires ~1hr; re-paste if you hit 401s). No secrets are
stored in the notebook itself.

See `CONVENTIONS.md` (reused verbatim from the D27 repo) for the standing
chart/notebook conventions this analysis follows.
