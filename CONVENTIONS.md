# Chart & Notebook Conventions

Standing conventions used throughout `notebook.py`, so future edits (by anyone, or by an
AI assistant) stay consistent with the rest of the analysis. Apply these by default to any
new infection/knockdown flow chart or micro-project added to this notebook, not just on
request.

- **Keep replicates separate, never pool/average into one bar or trace.** Every bar chart
  shows one bar per replicate well (Well ID / Sample ID in the tick label), and every
  histogram shows one clickable trace per well/sample rather than a pooled or
  group-averaged series. If multiple wells share a condition, show them as adjacent
  bars/traces rather than collapsing with `.mean()`.

- **Pair every infection/knockdown histogram with a companion FlowJo-style scatter plot,
  stacked below it (not side-by-side).** X-axis = the fluorescence channel (matching the
  histogram), Y-axis = SSC-A. One clickable/toggleable trace per sample. Downsample to
  ~2000 events per sample (a small gray footnote notes the cap, not per-trace "N shown"
  text) to keep the JSON payload small. Debris/doublets must already be excluded via the
  same live-gating used for the histogram. Reuse the shared `scatter_gate`/
  `scatter_gate_grid` helpers (paired with `interactive_hist`/`interactive_hist_grid`)
  rather than one-off scatter code per chart. Give the scatter's y-axis (SSC-A) explicit
  headroom above the tallest plotted point (`* 1.15`) rather than relying on default
  autorange — otherwise a pileup of cells at the channel's detector ceiling reads as a
  line hugging the plot's top edge.

- **All Plotly figures should have clickable/toggleable legends, one trace per
  population/sample** (click to hide, double-click to isolate).

- **Every knockdown/infection %-bar must show the numeric % value directly on the bar, not
  just the cell count.** Label each bar `{pct:.1f}%\n(n={n})` (via `ax.bar_label` or a
  per-bar `ax.text` at the bar's height).

- **Bar chart y-axis: floor of at least 50%, scale up only if real values exceed it.**
  `ax.set_ylim(0, max(50, values.max() * 1.25))` on every knockdown/infection bar chart —
  never let matplotlib autoscale tightly around small values, which makes a negligible
  effect look dramatic.

- **Histogram/scatter x-axis (or y-axis for scatter): extend past the pooled data's actual
  range so low-abundance positive tails aren't clipped**, using a `headroom=True` flag that
  multiplies the percentile-based upper bound by 1.5 — not an arbitrary log-fold multiplier
  disconnected from the data.

- **Every infection/knockdown gating chart needs a visible dashed threshold line**,
  matching the actual gate used to compute the %. On a log-scale Plotly axis, pass the raw
  (non-log-transformed) threshold value to `add_vline`/`add_hline` — Plotly auto-transforms
  shape positions on log axes same as trace data, so pre-computing `log10(threshold)`
  places the line off-screen. Always set `layer="above"` and an explicit `line_width` on
  threshold shapes, since default styling can render underneath semi-transparent fill
  traces in multi-panel grids.

- **Median-of-signal bar charts are not used** — prefer full distributions (histograms) or
  %-above-threshold summaries over median-signal summaries.

- **Keep multi-panel Plotly titles short and non-duplicated.** Don't repeat the
  group/condition name in both the shared title and each subplot title — say it once.

- **Metadata and infection-gating-summary tables are tucked into a collapsed
  `mo.accordion`**, not shown inline, to keep each micro-project section scannable. This
  applies per micro-project and per sub-experiment (nested accordions).

- **Exclude any gated population with fewer than 300 cells (`MIN_CELLS_PER_GATE`) from
  %-based summaries (knockdown, infection).** Drop it from the chart and add a small
  markdown note below listing which sample(s) were excluded and their actual n (via the
  shared `excluded_note(df, ...)` helper). If every sample in a chart falls below the
  floor, show a placeholder instead of an empty/misleading chart. Small-n gated
  populations (double-positive/multi-marker gates, infection-gated knockdown subsets,
  RFP+/GFP+-only gates especially) are noisy enough to be actively misleading, not just
  imprecise.

  **Exception: never drop control wells for low n.** Non-targeting (NT-GFP, NT-RNP), ORK
  (mock), unstained, or parental/no-guide wells are exempt from the 300-cell floor and
  should always be shown even below it (`is_control_label(text)` helper) — a control
  naturally having few cells pass an infection/guide-positive gate is expected, not a
  data-quality problem, and excluding it would hide a still-relevant reference point.

- **Summary tables (the final %-above-background table per micro-project) are tidied for
  readability, not shown as raw dataframes.** Round percentage columns to 1 decimal place,
  rename internal column names to human-readable headers, sort rows by metric/condition/
  well, and render via `mo.ui.table` (sortable/searchable) — via the shared
  `tidy_summary_table(df)` helper.

- **Overall notebook organization**, the standard template for every micro-project here and
  any future addition:
  - One top-level `mo.accordion` with one entry per micro-project/experiment (see the
    in-notebook "How to use this notebook" note at the top for navigation help).
  - Each micro-project starts with a short markdown "Experiment layout" + "Samples"
    overview note (2 short paragraphs, no well-level detail) before any tables or charts,
    plus any caveats about control quality specific to that micro-project.
  - If a micro-project bundles multiple distinct sub-experiments (e.g., different viral
    delivery constructs), split those into their own nested `mo.accordion` sections, each
    with its own small header + 1-paragraph description.
  - Order within each section: overview note → metadata/infection-summary accordion →
    flow distribution charts (histogram + scatter) → knockdown/infection bar charts →
    summary table.
