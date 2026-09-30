# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "numpy==2.4.6",
#     "pandas==3.0.3",
#     "plotly==6.9.0",
#     "requests==2.34.2",
#     "xlrd==2.0.2",
# ]
# ///

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium", auto_download=["html"])


@app.cell
def _(mo):
    mo.md("""
    # CRISPRi HLC Trial3 D28 (fixed) -- B2M knockdown analysis

    Three independent B2M knockdown arms from the fixed HLC Trial3 D28 sample
    sheet (`20260929_HLC_Trial3_D28_fixed`, corrected version), each using a
    different guide-delivery mode and/or cell line. Raw FCS files are parsed
    directly (no FlowJo CSV export), debris/doublets are excluded on FSC-A/
    FSC-Width, and knockdown is computed as the percentage-point shift in "% of
    cells below a fixed gate" on the raw (uncorrected) B2M antibody channel
    between a non-targeting guide (NT-GFP) and the B2M-targeting guide
    (B2M-GFP), within the guide/effector-delivery-gated population.

    **Note on this dataset vs. the sibling D27 (fixed) dataset:** D28 targets
    B2M (not CD81), has no CD81 antibody stain at all, and its original sample
    sheet contained data-entry errors (a stray CD81-guide virus, a BV421
    channel collision, and stray Dox-induction flags) that the user corrected
    by re-uploading the sample sheet mid-session -- this notebook was built
    against the corrected file. K562 wells (ASGR1-antibody-specificity checks
    only) and 2 guide-only/no-effector singleton wells (GFP compensation
    controls only) are excluded from the knockdown analysis, per the user.
    """)
    return


@app.cell(hide_code=True)
def _(arm1_summary, arm2_summary, arm3_summary, go, mo, pd):
    _summary_rows = [arm1_summary, arm2_summary, arm3_summary]
    _summary_df = pd.DataFrame(_summary_rows)
    _summary_df["knockdown_pp_1st"] = _summary_df["knockdown_pp_1st"].round(1)
    _summary_df["knockdown_pp_50th"] = _summary_df["knockdown_pp_50th"].round(1)
    _summary_df = _summary_df.rename(columns={
        "arm": "Arm", "cell_line": "Cell line", "mock_group": "Mock group", "guide_group": "Guide-active group",
        "knockdown_pp_1st": "Knockdown, 1st pctile (pp)", "knockdown_pp_50th": "Knockdown, 50th pctile (pp)",
        "flag": "Caveat",
    })
    top_summary_table = mo.ui.table(_summary_df, selection=None)

    _colors = {"1st": "red", "50th": "purple"}
    _k1 = "knockdown_pp_1st"
    _k50 = "knockdown_pp_50th"
    top_summary_fig = go.Figure()
    top_summary_fig.add_trace(go.Bar(
        x=[r["arm"] for r in _summary_rows], y=[r[_k1] for r in _summary_rows],
        name="1st pctile (tail effect)", marker_color=_colors["1st"],
        text=[f"{r[_k1]:.1f}pp" for r in _summary_rows], textposition="outside",
    ))
    top_summary_fig.add_trace(go.Bar(
        x=[r["arm"] for r in _summary_rows], y=[r[_k50] for r in _summary_rows],
        name="50th pctile (population shift)", marker_color=_colors["50th"],
        text=[f"{r[_k50]:.1f}pp" for r in _summary_rows], textposition="outside",
    ))
    _all_vals = [r[_k1] for r in _summary_rows] + [r[_k50] for r in _summary_rows]
    _ymax = max(50, max(v for v in _all_vals if v == v) * 1.25) if any(v == v for v in _all_vals) else 50
    top_summary_fig.update_layout(
        title="Bottom line across all 3 arms: knockdown (pp vs. mock), 1st vs. 50th percentile gate",
        barmode="group", yaxis_title="Knockdown (percentage points)", yaxis_range=[0, _ymax],
        height=420, margin=dict(t=60),
        legend=dict(itemclick="toggle", itemdoubleclick="toggleothers"),
    )

    top_summary = mo.vstack([
        mo.md(
            "### Bottom line across all 3 arms\n"
            "Primary guide-active group vs. this arm's own mock/control group, at both the "
            "1st-percentile-of-mock (tail effect) and 50th-percentile-of-mock (population-level shift) gates. "
            "See each arm's tab below for full per-replicate detail, ASGR1+ stratification, and caveats."
        ),
        top_summary_fig,
        top_summary_table,
    ])
    top_summary
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ### How to read each panel below
    Each tab is one construct/cell-line arm: a short description of the delivery
    biology, well metadata, a peak-normalized histogram + FlowJo-style scatter of
    the B2M readout channel (gated on guide/effector-delivery markers, dashed red
    line = 1st-percentile-of-mock gate, dashed purple line = 50th-percentile-of-mock
    gate), a per-replicate bar chart for each metric, and a tidy summary table.
    Nested "All cells"/"ASGR1+" sub-tabs stratify by ASGR1-BV421 where available.
    """)
    return


@app.cell(hide_code=True)
def _(arm1_content, arm2_content, arm3_content, mo):
    mo.ui.tabs({
        "1. WTC11 transient B2M-GFP virus": arm1_content,
        "2. 17_3 stably-integrated guide + AA239 effector": arm2_content,
        "3. 17_3 fully-transient two-virus + AA239 effector": arm3_content,
    })
    return


@app.cell(hide_code=True)
def _(mo, pctile_gate):
    # One calibration slider per distinct (channel, reference-well) pair actually
    # used somewhere in the pipeline as a guide/effector-marker gate. Channel
    # assignments were empirically verified against single-stain control wells
    # (see exploration notes): B2M-APC->APC-A, B2M-PE-Dazzle-594->B610-ECD-A,
    # Thy1.1-APC->APC-A, ASGR1-BV421->BV421-A, NT-GFP/B2M-GFP marker->FITC-A.
    # D11 (fully triple-stained: B2M-PE-Dazzle/Thy1.1-APC/ASGR1-BV421, no guide,
    # no virus) is used as an FMO-style stained-but-uninfected reference for the
    # Thy1.1-APC effector-delivery gate, shared by both 17_3 two-virus/integrated
    # arms, per the FMO-gating convention.
    _infection_gate_defaults = {
        "FITC-A__A3": pctile_gate("A3", "FITC-A", 99),
        "FITC-A__A7": pctile_gate("A7", "FITC-A", 99),
        "APC-A__D11_fmo": pctile_gate("D11", "APC-A", 99),
        "BV421-A__A3": pctile_gate("A3", "BV421-A", 99),
        "BV421-A__A7": pctile_gate("A7", "BV421-A", 99),
    }

    def _widget_label(key):
        channel, ref = key.split("__")
        if ref.endswith("_fmo"):
            return f"{channel} infection-gate threshold (FMO ref = well {ref.replace('_fmo', '')}, stained/uninfected)"
        return f"{channel} gate threshold (naive/unstained ref = well {ref})"

    infection_gate_widgets = mo.ui.dictionary({
        key: mo.ui.slider(
            start=0, stop=max(2000.0, round(val * 6)), step=max(50.0, round(val / 100)),
            value=round(val), show_value=True, full_width=True,
            label=_widget_label(key),
        )
        for key, val in _infection_gate_defaults.items()
    })
    infection_gate_widgets
    return (infection_gate_widgets,)


@app.cell(hide_code=True)
def _(
    debris_gate,
    go,
    infection_gate_widgets,
    mo,
    np,
    pd,
    replicate_labels,
    sample_sheet,
):
    # Arm 2's 12 wells (C1-C8, D1-D4) actually compare 3 lentiviral delivery
    # reagents for the AA239 effector-virus transduction (Polybrene, Protamine
    # Sulfate, Lentiboost), per the sample sheet's "Transduction Reagent" column --
    # 4 wells each (2 NT-GFP + 2 B2M-GFP), not surfaced in the original arm
    # summary. This chart lives at the top of Arm 2's tab, in the infection-gating
    # section (reagent choice affects transduction/infection efficiency, not
    # knockdown per se) -- separate from the knockdown-metric charts below.
    # Reactive: recomputes from infection_gate_widgets (the same Thy1.1-APC
    # effector-delivery FMO gate used everywhere else for this arm) and the
    # debris/doublet sliders, exactly like the other infection-gate charts.
    _ARM2_WELLS = ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "D1", "D2", "D3", "D4"]
    _arm2_thy11_gate = infection_gate_widgets.value["APC-A__D11_fmo"]

    _reagent_rows = []
    for _w in _ARM2_WELLS:
        _d = debris_gate(_w)
        _n = len(_d)
        _pct_pos = 100 * float(np.mean(_d["APC-A"].values > _arm2_thy11_gate)) if _n else np.nan
        _reagent_rows.append({
            "well": _w,
            "reagent": sample_sheet.loc[_w, "Transduction Reagent"],
            "guide": "NT-GFP" if sample_sheet.loc[_w, "Integrated Guide"] == "NT-GFP" else "B2M-GFP",
            "n": _n,
            "pct_infected": _pct_pos,
        })
    arm2_reagent_df = pd.DataFrame(_reagent_rows)

    _reagent_colors = {"Polybrene": "#4C78A8", "Protamine Sulfate": "#E45756", "Lentiboost": "#54A24B"}
    _rl = replicate_labels(_ARM2_WELLS)
    arm2_reagent_fig = go.Figure()
    for _reagent, _grp in arm2_reagent_df.groupby("reagent", sort=False):
        arm2_reagent_fig.add_trace(go.Bar(
            x=[_rl[w] for w in _grp["well"]], y=_grp["pct_infected"],
            name=_reagent, marker_color=_reagent_colors.get(_reagent, "#999"),
            text=[f"{v:.1f}%<br>(n={n:,})" for v, n in zip(_grp["pct_infected"], _grp["n"])],
            textposition="outside",
        ))
    arm2_reagent_fig.update_layout(
        title=f"Arm 2: % Thy1.1-APC+ (effector-infected) per replicate, by transduction reagent (gate = {_arm2_thy11_gate:,.0f})",
        yaxis_title="% infected (Thy1.1-APC+)", yaxis_range=[0, 100],
        height=380, margin=dict(t=60),
        legend_title="Transduction reagent",
    )

    arm2_reagent_section = mo.vstack([
        mo.md(
            "### Arm 2 infection-gating: transduction-reagent comparison\n"
            "Arm 2's 12 wells also compare 3 lentiviral delivery reagents for the "
            "AA239 effector virus (`Transduction Reagent` column in the sample "
            "sheet): **Polybrene** (C1-C4), **Protamine Sulfate** (C5-C8), "
            "**Lentiboost** (D1-D4) -- 4 wells each (2 NT-GFP + 2 B2M-GFP). This "
            "chart is reactive to the Thy1.1-APC effector-delivery gate slider "
            "above (`infection_gate_widgets`), same as the other infection-gate "
            "content in this notebook -- move that slider to see the infected "
            "fractions below update live."
        ),
        arm2_reagent_fig,
    ])
    arm2_reagent_section
    return (arm2_reagent_section,)


@app.cell(hide_code=True)
def _(
    biexp,
    biexp_ticks,
    calibration_trace_labels,
    debris_gate,
    gated_range,
    infection_gate_widgets,
    interactive_hist,
    mo,
    np,
    sample_sheet,
):
    _CALIBRATION_CFGS = [
        {"key": "FITC-A__A3", "channel": "FITC-A", "naive": "A3",
         "used_by": {"Arm 1 GFP guide marker (B1-B6)": ["B1", "B2", "B3", "B4", "B5", "B6"]}},
        {"key": "FITC-A__A7", "channel": "FITC-A", "naive": "A7",
         "used_by": {"Arm 2 GFP guide marker (C1-C8,D1-D4)": ["C1","C2","C3","C4","C5","C6","C7","C8","D1","D2","D3","D4"],
                     "Arm 3 GFP guide marker (D7-D10)": ["D7", "D8", "D9", "D10"]}},
        {"key": "APC-A__D11_fmo", "channel": "APC-A", "naive": "D11", "secondary_naive": "A7",
         "used_by": {"Arm 2 Thy1.1 effector marker (C1-C8,D1-D4)": ["C1","C2","C3","C4","C5","C6","C7","C8","D1","D2","D3","D4"],
                     "Arm 3 Thy1.1 effector marker (D7-D10)": ["D7", "D8", "D9", "D10"]}},
        {"key": "BV421-A__A3", "channel": "BV421-A", "naive": "A3",
         "used_by": {"Arm 1 ASGR1 marker (B1-B6)": ["B1", "B2", "B3", "B4", "B5", "B6"]}},
        {"key": "BV421-A__A7", "channel": "BV421-A", "naive": "A7",
         "used_by": {"Arm 2 ASGR1 marker (C1-C8,D1-D4)": ["C1","C2","C3","C4","C5","C6","C7","C8","D1","D2","D3","D4"],
                     "Arm 3 ASGR1 marker (D7-D10)": ["D7", "D8", "D9", "D10"]}},
    ]

    def build_calibration_panel(cfg):
        channel, naive, key = cfg["channel"], cfg["naive"], cfg["key"]
        slider = infection_gate_widgets.elements[key]
        gate = slider.value
        naive_vals = debris_gate(naive)[channel].values
        naive_cell_line = sample_sheet.loc[naive, "Cell Line"]
        is_fmo = key.endswith("_fmo")
        if is_fmo:
            naive_label = f"FMO ref ({naive_cell_line}, stained/uninfected, well {naive})"
        else:
            naive_label = f"Naive ({naive_cell_line}, unstained, well {naive})"

        traces = {naive_label: naive_vals}
        pooled_for_range = [naive_vals]

        secondary = cfg.get("secondary_naive")
        if secondary:
            sec_vals = debris_gate(secondary)[channel].values
            sec_cell_line = sample_sheet.loc[secondary, "Cell Line"]
            sec_label = f"(for comparison) Naive ({sec_cell_line}, unstained, well {secondary})"
            traces[sec_label] = sec_vals
            pooled_for_range.append(sec_vals)

        all_sample_wells = [w for wells in cfg["used_by"].values() for w in wells]
        well_labels = calibration_trace_labels(all_sample_wells)

        for w in all_sample_wells:
            vals = debris_gate(w)[channel].values
            traces[well_labels[w]] = vals
            pooled_for_range.append(vals)
        pooled = np.concatenate(pooled_for_range)
        xr = gated_range(naive_vals, pooled, hi_pct=99.0)
        cofactor = max(gate, 1.0) / 5.0

        def _tx(v):
            return biexp(v, cofactor)

        traces_t = {k: _tx(v) for k, v in traces.items()}
        xr_t = (float(_tx(xr[0])), float(_tx(xr[1])))
        gate_t = float(_tx(gate))
        fig = interactive_hist(
            traces_t, gate_t, f"gate = {gate:,.0f}", xr_t,
            f"{channel} gate calibration (reference = {naive_label})",
        )
        tickvals, ticktext = biexp_ticks(xr[0], xr[1], cofactor)
        fig.update_layout(xaxis=dict(tickvals=tickvals, ticktext=ticktext, title=f"{channel} (biexponential scale)"))

        blocks = [mo.md(f"**Used by:** {', '.join(cfg['used_by'].keys())}"), slider, fig]

        if is_fmo:
            blocks.append(mo.md(
                "*FMO-style reference: this gate is calibrated off a "
                "**stained-but-uninfected** control (full antibody panel applied, "
                "no guide/effector virus delivered, well D11) rather than a fully-"
                "unstained naive well -- per the FMO-gating convention. The fully-"
                "unstained trace is kept visible above for comparison.*"
            ))

        return mo.vstack(blocks)

    infection_calibration_accordion = mo.ui.tabs({
        (
            f"{cfg['channel']} FMO ref: {sample_sheet.loc[cfg['naive'], 'Cell Line']} (stained/uninfected)"
            if cfg["key"].endswith("_fmo")
            else f"{cfg['channel']} naive: {sample_sheet.loc[cfg['naive'], 'Cell Line']} (unstained)"
        ): build_calibration_panel(cfg)
        for cfg in _CALIBRATION_CFGS
    })

    infection_calibration_section = mo.vstack([
        mo.md(
            "## Infection/guide-delivery gate calibration\n"
            "FlowJo-style gating histograms (biexponential x-axis, same "
            "`arcsinh(x/cofactor)` transform used for the readout channels below) for "
            "every guide/effector-marker and ASGR1-stratification gate actually used "
            "in the pipeline. The Thy1.1-APC effector-delivery gate (arms 2/3) uses an "
            "FMO-style stained-but-uninfected reference well (D11) rather than a fully-"
            "unstained naive. **The sliders below are the live source of truth** for "
            "these gates -- moving one reactively re-runs every arm that uses that "
            "(channel, reference-well) pair and updates its knockdown numbers below."
        ),
        infection_calibration_accordion,
    ])
    infection_calibration_section
    return


@app.cell(hide_code=True)
def _(biexp, biexp_ticks, debris_gate, gated_range, interactive_hist, mo, np):
    _HEP_REFS = [
        {"label": "WTC11 (no virus)", "well": "H_A1"},
        {"label": "17_3 (integrated NT-GFP, no effector)", "well": "H_A2"},
        {"label": "17_3 (integrated B2M-GFP, no effector)", "well": "H_A3"},
        {"label": "17_3 (no guide, no effector)", "well": "H_A4"},
    ]

    def _hep_panel(channel, title):
        traces = {}
        pooled = []
        for ref in _HEP_REFS:
            vals = debris_gate(ref["well"])[channel].values
            traces[f"{ref['label']} ({ref['well']})"] = vals
            pooled.append(vals)
        pooled = np.concatenate(pooled)
        xr = gated_range(pooled, pooled, hi_pct=99.0)
        cofactor = max(float(np.median(np.abs(pooled))), 1.0) / 5.0

        def _tx(v):
            return biexp(v, cofactor)

        traces_t = {k: _tx(v) for k, v in traces.items()}
        xr_t = (float(_tx(xr[0])), float(_tx(xr[1])))
        fig = interactive_hist(traces_t, None, "", xr_t, title)
        tickvals, ticktext = biexp_ticks(xr[0], xr[1], cofactor)
        fig.update_layout(xaxis=dict(tickvals=tickvals, ticktext=ticktext, title=f"{channel} (biexponential scale)"))
        return fig

    hepatocyte_marker_tabs = mo.ui.tabs({
        "ASGR1-PE (informational)": _hep_panel("PE-A", "ASGR1-PE across cell-line reference wells (informational, not used for arm stratification)"),
        "Albumin-AF647 (informational)": _hep_panel("APC-A", "Albumin-AF647 across cell-line reference wells (informational, not used for arm stratification)"),
    })

    hepatocyte_marker_section = mo.vstack([
        mo.md(
            "## Hepatocyte marker gating (ASGR1, Albumin) -- informational only\n"
            "**No arm below is stratified by ASGR1-PE or Albumin-AF647** -- those "
            "antibodies were only stained on a separate set of reference wells "
            "(H_A1-H_A4), never co-stained together with the knockdown-arm readout "
            "antibodies. The tabs below exist purely to show what signal looks like "
            "on those reference wells, for context/QC. **ASGR1 stratification for "
            "the knockdown arms instead uses ASGR1-BV421**, which *is* co-stained "
            "in-well for all 3 arms -- see each arm's own \"ASGR1+\" sub-tab below. "
            "This dataset does not contain the previously-flagged non-specific "
            "ASGR1-FITC clone at all (ASGR1 here is BV421 or PE only)."
        ),
        hepatocyte_marker_tabs,
    ])
    hepatocyte_marker_section
    return


@app.cell(hide_code=True)
def _(debris_slider, go, mo, np, raw_wells):
    _pooled_fsc = np.concatenate([d["FSC-A"].values for d in raw_wells.values()])
    _pooled_fsc = _pooled_fsc[(_pooled_fsc > -50_000) & (_pooled_fsc < 3_000_000)]
    _fig = go.Figure()
    _fig.add_trace(go.Histogram(x=_pooled_fsc, nbinsx=150, marker_color="#4C78A8"))
    _fig.add_shape(
        type="line", x0=debris_slider.value, x1=debris_slider.value, y0=0, y1=1,
        yref="paper", line=dict(color="red", width=3, dash="dash"), layer="above",
    )
    _fig.add_annotation(
        x=debris_slider.value, y=1.03, yref="paper",
        text=f"debris gate = {debris_slider.value:,.0f}", showarrow=False,
        font=dict(color="red"),
    )
    _fig.update_layout(
        title="Pooled FSC-A (all loaded wells, fixed HLC Trial3 D28) -- debris/cell valley check",
        xaxis_title="FSC-A (linear)", yaxis_title="count", height=380, margin=dict(t=60),
    )
    debris_check = mo.vstack([
        mo.md(
            "**Mandatory debris-gate recheck.** SSC-A is only floored at >0 (no separate "
            "debris criterion). Confirm the dashed line below sits in the FSC-A valley "
            "between debris and real cells before trusting any downstream gate -- this "
            "dataset is a different day/run than D27 and the valley location may have shifted."
        ),
        _fig,
    ])
    debris_check
    return


@app.cell(hide_code=True)
def _(debris_slider, go, mo, np, raw_wells):
    _rng_d = np.random.default_rng(0)
    _xs_d, _ys_d = [], []
    for _d in raw_wells.values():
        _fsc = _d["FSC-A"].values
        _ssc = _d["SSC-A"].values
        # Same outlier clip as the debris_check histogram (FSC-A) plus an analogous
        # clip on SSC-A -- roughly 12% of events reach up to ~9e7 on SSC-A here
        # (vs. a real population sitting in 0-3e6), which without clipping stretches
        # Plotly's autorange so far that the real population is compressed into an
        # invisible sliver near the origin (looked like an "empty" plot, but the
        # points were always there -- just off-scale).
        _mask = (_fsc > -50_000) & (_fsc < 3_000_000) & (_ssc > -50_000) & (_ssc < 3_000_000)
        _fsc, _ssc = _fsc[_mask], _ssc[_mask]
        _n = len(_fsc)
        _cap = min(_n, 2500)
        if _cap == 0:
            continue
        _idx = _rng_d.choice(_n, _cap, replace=False)
        _xs_d.append(_fsc[_idx])
        _ys_d.append(_ssc[_idx])
    _xs_d = np.concatenate(_xs_d)
    _ys_d = np.concatenate(_ys_d)
    _fig_ds = go.Figure()
    _fig_ds.add_trace(go.Scattergl(
        x=_xs_d, y=_ys_d, mode="markers",
        marker=dict(size=2, opacity=0.25, color="#4C78A8"),
        name="events",
    ))
    _fig_ds.add_shape(
        type="line", x0=debris_slider.value, x1=debris_slider.value, y0=0, y1=1, yref="paper",
        line=dict(color="red", width=3, dash="dash"), layer="above",
    )
    _fig_ds.add_annotation(
        x=debris_slider.value, y=1.03, yref="paper",
        text=f"FSC-A debris gate = {debris_slider.value:,.0f}", showarrow=False, font=dict(color="red"),
    )
    _fig_ds.update_layout(
        title="Pooled FSC-A vs SSC-A -- debris/cell valley in 2D",
        xaxis_title="FSC-A", yaxis_title="SSC-A", height=420, margin=dict(t=60),
        xaxis_range=[-50_000, 3_000_000], yaxis_range=[-50_000, 3_000_000],
    )
    _fig_ds.add_annotation(
        text="downsampled to <=2,500 events/well for display, outliers >3e6 on either axis clipped", x=0, y=-0.14,
        xref="paper", yref="paper", showarrow=False, font=dict(size=10, color="gray"),
    )
    debris_scatter = mo.vstack([_fig_ds])
    debris_scatter
    return


@app.cell(hide_code=True)
def _(doublet_slider, go, mo, np, raw_wells):
    _pooled_w = np.concatenate([d["FSC-Width"].values for d in raw_wells.values()])
    _pooled_w = _pooled_w[(_pooled_w > 0) & (_pooled_w < 8000)]
    _fig_w = go.Figure()
    _fig_w.add_trace(go.Histogram(x=_pooled_w, nbinsx=150, marker_color="#54A24B"))
    _fig_w.add_shape(
        type="line", x0=doublet_slider.value, x1=doublet_slider.value, y0=0, y1=1,
        yref="paper", line=dict(color="red", width=3, dash="dash"), layer="above",
    )
    _fig_w.add_annotation(
        x=doublet_slider.value, y=1.03, yref="paper",
        text=f"doublet cutoff = {doublet_slider.value:,.0f}", showarrow=False,
        font=dict(color="red"),
    )
    _fig_w.update_layout(
        title="Pooled FSC-Width (all loaded wells) -- doublet/clump cutoff check",
        xaxis_title="FSC-Width (linear)", yaxis_title="count", height=380, margin=dict(t=60),
    )
    doublet_check = mo.vstack([
        mo.md(
            "**Doublet/clump discrimination.** Confirm the cutoff below sits past the "
            "main singlet peak, clipping the doublet tail, before trusting downstream gates."
        ),
        _fig_w,
    ])
    doublet_check
    return


@app.cell(hide_code=True)
def _(doublet_slider, go, mo, np, raw_wells):
    _rng_w = np.random.default_rng(0)
    _xs_w, _ys_w = [], []
    for _d in raw_wells.values():
        _fsc = _d["FSC-A"].values
        _wid = _d["FSC-Width"].values
        # Same outlier clip as the doublet_check histogram (FSC-Width in (0, 8000))
        # plus an analogous clip on FSC-A -- without this, extreme FSC-A outliers
        # stretch the autorange so far that the real singlet/doublet population is
        # compressed into an invisible sliver (looked "empty", points were off-scale).
        _mask = (_fsc > -50_000) & (_fsc < 3_000_000) & (_wid > 0) & (_wid < 8000)
        _fsc, _wid = _fsc[_mask], _wid[_mask]
        _n = len(_fsc)
        _cap = min(_n, 2500)
        if _cap == 0:
            continue
        _idx = _rng_w.choice(_n, _cap, replace=False)
        _xs_w.append(_fsc[_idx])
        _ys_w.append(_wid[_idx])
    _xs_w = np.concatenate(_xs_w)
    _ys_w = np.concatenate(_ys_w)
    _fig_ws = go.Figure()
    _fig_ws.add_trace(go.Scattergl(
        x=_xs_w, y=_ys_w, mode="markers",
        marker=dict(size=2, opacity=0.25, color="#54A24B"),
        name="events",
    ))
    _fig_ws.add_shape(
        type="line", x0=0, x1=1, xref="paper", y0=doublet_slider.value, y1=doublet_slider.value,
        line=dict(color="red", width=3, dash="dash"), layer="above",
    )
    _fig_ws.add_annotation(
        x=0.02, y=doublet_slider.value, xref="paper",
        text=f"FSC-Width doublet cutoff = {doublet_slider.value:,.0f}", showarrow=False,
        font=dict(color="red"), xanchor="left", yanchor="bottom",
    )
    _fig_ws.update_layout(
        title="Pooled FSC-A vs FSC-Width -- doublet band above the cutoff",
        xaxis_title="FSC-A", yaxis_title="FSC-Width", height=420, margin=dict(t=60),
        xaxis_range=[-50_000, 3_000_000], yaxis_range=[0, 8000],
    )
    _fig_ws.add_annotation(
        text="downsampled to <=2,500 events/well for display, outliers clipped (FSC-A>3e6, FSC-Width>8000)", x=0, y=-0.14,
        xref="paper", yref="paper", showarrow=False, font=dict(size=10, color="gray"),
    )
    doublet_scatter = mo.vstack([_fig_ws])
    doublet_scatter
    return


@app.cell(hide_code=True)
def _(mo):
    debris_slider = mo.ui.slider(
        start=0, stop=1_500_000, step=10_000, value=250_000,
        label="FSC-A debris/cell gate (events at/below this FSC-A value are excluded as debris)",
        full_width=True, show_value=True,
    )
    debris_slider
    return (debris_slider,)


@app.cell(hide_code=True)
def _(mo):
    doublet_slider = mo.ui.slider(
        start=0, stop=8000, step=100, value=3800,
        label="FSC-Width doublet/clump cutoff (events above this FSC-Width are excluded as doublets)",
        full_width=True, show_value=True,
    )
    doublet_slider
    return (doublet_slider,)


@app.cell(hide_code=True)
def _():
    import marimo as mo
    import numpy as np
    import pandas as pd
    import requests
    import urllib.parse
    import io
    import plotly.graph_objects as go

    return go, io, mo, np, pd, requests, urllib


@app.cell(hide_code=True)
def _(mo):
    token_input = mo.ui.text(
        kind="password",
        label="GCS access token (paste output of `gcloud auth print-access-token`; expires ~1hr, re-paste if you hit 401s)",
        full_width=True,
    )
    token_input
    return (token_input,)


@app.cell(hide_code=True)
def _(requests, token_input, urllib):
    BUCKET = "20240617-landerlab-crispri-project"
    PREFIX = "Flow_Data_CRISPRi/20260929_HLC_Trial3_D28_fixed/"

    def gcs_get(name: str) -> bytes:
        enc = urllib.parse.quote(name, safe="")
        url = f"https://storage.googleapis.com/storage/v1/b/{BUCKET}/o/{enc}?alt=media"
        headers = {"Authorization": f"Bearer {token_input.value}"}
        r = requests.get(url, headers=headers)
        r.raise_for_status()
        return r.content

    def gcs_list(prefix: str):
        url = f"https://storage.googleapis.com/storage/v1/b/{BUCKET}/o"
        headers = {"Authorization": f"Bearer {token_input.value}"}
        items = []
        params = {"prefix": prefix}
        while True:
            r = requests.get(url, headers=headers, params=params)
            r.raise_for_status()
            j = r.json()
            items.extend(j.get("items", []))
            if "nextPageToken" not in j:
                break
            params["pageToken"] = j["nextPageToken"]
        return items

    return PREFIX, gcs_get


@app.cell(hide_code=True)
def _(np, pd):
    def parse_fcs_bytes(data: bytes):
        text_start = int(data[10:18].decode().strip())
        text_end = int(data[18:26].decode().strip())
        data_start_hdr = int(data[26:34].decode().strip())
        data_end_hdr = int(data[34:42].decode().strip())
        text = data[text_start:text_end + 1].decode("latin-1")
        delim = text[0]
        parts = text[1:].split(delim)
        kv = {}
        it = iter(parts)
        for k, v in zip(it, it):
            kv[k.upper()] = v
        data_start = data_start_hdr or int(kv["$BEGINDATA"])
        data_end = data_end_hdr or int(kv["$ENDDATA"])
        npar = int(kv["$PAR"]); ntot = int(kv["$TOT"])
        endian = "<" if kv["$BYTEORD"].startswith("1,2") else ">"
        dt_map = {"F": "f4", "D": "f8", "I": "u4"}
        dtype = np.dtype(endian + dt_map[kv["$DATATYPE"]])
        raw = np.frombuffer(data[data_start:data_end + 1], dtype=dtype).reshape(ntot, npar)
        names = [kv.get(f"$P{i}S") or kv.get(f"$P{i}N", f"P{i}") for i in range(1, npar + 1)]
        return pd.DataFrame(raw, columns=names), kv

    return (parse_fcs_bytes,)


@app.cell(hide_code=True)
def _(PREFIX, gcs_get, io, pd):
    _xls_bytes = gcs_get(PREFIX + "20260929_SampleID copy.xls")
    sample_sheet = pd.read_excel(io.BytesIO(_xls_bytes), engine="xlrd", header=0)

    def _well_key(well_id: str) -> str:
        # Plate 2 ("02-Well-*") wells collide on bare code with plate 1 ("01-Well-*")
        # once the "0N-Well-" prefix is stripped (e.g. both have an "A1") -- prefix
        # plate-2 keys with "H_" (hepatocyte co-stain plate) to keep them unique and
        # matching the raw_wells dict keys used below.
        if well_id.startswith("02-Well-"):
            return "H_" + well_id.replace("02-Well-", "")
        return well_id.replace("01-Well-", "")

    sample_sheet["Well"] = sample_sheet["Well ID"].apply(_well_key)
    sample_sheet = sample_sheet.set_index("Well", drop=False)
    sample_sheet
    return (sample_sheet,)


@app.cell(hide_code=True)
def _(PREFIX, gcs_get, parse_fcs_bytes):
    # Wells used across the 3 knockdown arms + naive/single-stain/FMO references.
    # Excludes: K562 wells (02-Well-A5/A6/A7, ASGR1-antibody-specificity check only,
    # not part of the knockdown analysis per the user) and D5/D6 (guide-only,
    # no-effector singleton wells -- GFP compensation controls only, not meant for
    # replicate-based statistics, per the user).
    WELLS = [
        # WTC11 arm (transient B2M-GFP virus, no separate effector virus)
        "A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "B5", "B6", "E6", "E7",
        # 17_3 arms (integrated-guide+AA239 effector, and fully-transient two-virus)
        "A5", "A6", "A7", "A8",
        "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8",
        "D1", "D2", "D3", "D4", "D7", "D8", "D9", "D10", "D11",
        "E1", "E2", "E3", "E4", "E5", "E8", "E9", "E10",
    ]
    # Hepatocyte co-stain wells (ASGR1-PE + Albumin-AF647), informational only --
    # not used for in-arm stratification since these antibodies were not co-stained
    # in the arm wells themselves (only ASGR1-BV421 was, within the arm wells).
    HEP_WELLS = ["02-Well-A1", "02-Well-A2", "02-Well-A3", "02-Well-A4"]

    def _load_well(w):
        _data = gcs_get(PREFIX + f"01-Well-{w}.fcs")
        _df, _kv = parse_fcs_bytes(_data)
        return _df

    def _load_hep_well(name):
        _data = gcs_get(PREFIX + f"{name}.fcs")
        _df, _kv = parse_fcs_bytes(_data)
        return _df

    raw_wells = {w: _load_well(w) for w in WELLS}
    raw_wells.update({name.replace("02-Well-", "H_"): _load_hep_well(name) for name in HEP_WELLS})
    f"Loaded {len(raw_wells)} wells, {sum(len(d) for d in raw_wells.values()):,} total events"
    return (raw_wells,)


@app.cell(hide_code=True)
def _(debris_slider, doublet_slider, go, mo, np, pd, raw_wells, sample_sheet):
    MIN_CELLS_PER_GATE = 300
    _CONTROL_TOKENS = ["nt-", "non-targeting", "unstained", "no guide", "control"]

    def is_control_label(text: str) -> bool:
        t = (text or "").lower()
        return any(tok in t for tok in _CONTROL_TOKENS)

    def debris_gate(well: str) -> pd.DataFrame:
        d = raw_wells[well]
        return d[
            (d["FSC-A"] > debris_slider.value)
            & (d["SSC-A"] > 0)
            & (d["FSC-Width"] <= doublet_slider.value)
        ]

    def pctile_gate(well: str, chan: str, pct: float = 99.0) -> float:
        d = debris_gate(well)
        return float(np.percentile(d[chan].values, pct))

    def apply_markers(well: str, marker_specs: list) -> pd.DataFrame:
        d = debris_gate(well)
        mask = np.ones(len(d), dtype=bool)
        for chan, thr in marker_specs:
            mask &= d[chan].values > thr
        return d[mask]

    def gated_range(naive_vals, pooled_vals, pad_frac: float = 0.10, hi_pct: float = 99.0):
        naive_vals = np.asarray(naive_vals)
        pooled_vals = np.asarray(pooled_vals)
        lo_anchor = np.percentile(naive_vals, 0.5) if len(naive_vals) else np.percentile(pooled_vals, 0.5)
        hi_anchor = np.percentile(pooled_vals, hi_pct)
        span = max(hi_anchor - lo_anchor, 1.0)
        return lo_anchor - pad_frac * span, hi_anchor + pad_frac * span

    def add_threshold(fig, x, label, color="red", y=1.05):
        fig.add_shape(
            type="line", x0=x, x1=x, y0=0, y1=1, yref="paper",
            line=dict(color=color, width=3, dash="dash"), layer="above",
        )
        fig.add_annotation(
            x=x, y=y, yref="paper", text=label, showarrow=False, font=dict(color=color),
        )

    def interactive_hist(
        traces: dict, gate, gate_label: str, xrange, title: str, nbins: int = 60,
        gate2=None, gate2_label: str = "", gate2_color: str = "purple",
    ):
        fig = go.Figure()
        colors = ["#4C78A8", "#E45756", "#54A24B", "#F58518", "#B279A2", "#72B7B2"]
        for i, (label, vals) in enumerate(traces.items()):
            vals = np.asarray(vals)
            vals = vals[(vals >= xrange[0]) & (vals <= xrange[1])]
            counts, edges = np.histogram(vals, bins=nbins, range=xrange)
            counts = counts / counts.max() * 100 if counts.max() > 0 else counts
            centers = (edges[:-1] + edges[1:]) / 2
            fig.add_trace(go.Scatter(
                x=centers, y=counts, mode="lines", fill="tozeroy",
                name=f"{label} (n={len(vals):,})", line=dict(color=colors[i % len(colors)]),
                opacity=0.55,
            ))
        if gate is not None:
            add_threshold(fig, gate, gate_label, color="red", y=1.05)
        if gate2 is not None:
            add_threshold(fig, gate2, gate2_label, color=gate2_color, y=1.12)
        fig.update_layout(
            title=title, xaxis_title="signal", yaxis_title="% of peak (normalized)",
            height=380, margin=dict(t=70), xaxis_range=list(xrange),
            legend=dict(itemclick="toggle", itemdoubleclick="toggleothers"),
        )
        return fig

    def scatter_gate(
        traces: dict, xrange, title: str, y_chan: str = "SSC-A", n_show: int = 2000, seed: int = 0,
        gate=None, gate_label: str = "",
        gate2=None, gate2_label: str = "", gate2_color: str = "purple",
    ):
        fig = go.Figure()
        colors = ["#4C78A8", "#E45756", "#54A24B", "#F58518", "#B279A2", "#72B7B2"]
        rng = np.random.default_rng(seed)
        for i, (label, (xv, yv)) in enumerate(traces.items()):
            xv = np.asarray(xv); yv = np.asarray(yv)
            if len(xv) > n_show:
                idx = rng.choice(len(xv), n_show, replace=False)
                xv, yv = xv[idx], yv[idx]
            fig.add_trace(go.Scattergl(
                x=xv, y=yv, mode="markers", name=label,
                marker=dict(size=3, opacity=0.4, color=colors[i % len(colors)]),
            ))
        if gate is not None:
            add_threshold(fig, gate, gate_label, color="red", y=1.05)
        if gate2 is not None:
            add_threshold(fig, gate2, gate2_label, color=gate2_color, y=1.12)
        fig.update_layout(
            title=title, xaxis_title="signal", yaxis_title=y_chan, height=380,
            margin=dict(t=70), xaxis_range=list(xrange),
            legend=dict(itemclick="toggle", itemdoubleclick="toggleothers"),
        )
        fig.add_annotation(
            text="downsampled to <=2,000 events/sample for display", x=0, y=-0.18,
            xref="paper", yref="paper", showarrow=False, font=dict(size=10, color="gray"),
        )
        return fig

    def tidy_summary_table(df: pd.DataFrame):
        _df = df.copy()
        for c in _df.columns:
            if "pct" in c.lower():
                _df[c] = _df[c].round(1)
        rename = {
            "well": "Well", "label": "Condition", "n": "n (gated cells)",
            "pct_below_gate": "% below 1st-pctile-mock gate",
            "knockdown_pp": "Knockdown, 1st pctile (pp vs mock)",
            "pct_below_gate50": "% below 50th-pctile-mock gate",
            "knockdown_pp_50": "Knockdown, 50th pctile (pp vs mock)",
        }
        _df = _df.rename(columns={k: v for k, v in rename.items() if k in _df.columns})
        return mo.ui.table(_df, selection=None)

    def replicate_labels(wells: list) -> dict:
        labels = {}
        for i, w in enumerate(wells, start=1):
            cell_line = sample_sheet.loc[w, "Cell Line"]
            labels[w] = f"{cell_line} rep{i} ({w})"
        return labels

    def biexp(x, cofactor: float):
        """arcsinh-based biexponential display transform: linear-like near zero,
        log-like at high magnitude, handles negative/near-zero compensated values
        gracefully (unlike a true log axis)."""
        return np.arcsinh(np.asarray(x, dtype=float) / cofactor)

    def biexp_ticks(lo_raw: float, hi_raw: float, cofactor: float):
        vals = [0.0]
        mult = 1.0
        while True:
            v = cofactor * mult
            if v > hi_raw * 1.02:
                break
            vals.append(v)
            mult *= 3.0
        if lo_raw < -cofactor * 0.5:
            vals = [lo_raw] + vals
        vals = sorted(set(round(v, 0) for v in vals))
        tickvals = [biexp(v, cofactor) for v in vals]
        ticktext = [f"{v:,.0f}" for v in vals]
        return tickvals, ticktext

    def _condition_label(well: str) -> str:
        row = sample_sheet.loc[well]
        cell_line = str(row["Cell Line"])

        def _clean(x):
            s = str(x).strip() if x is not None else ""
            return s if s and s.lower() not in ("nan", "no", "none") else None

        parts = [cell_line]
        guide_c = _clean(row.get("Integrated Guide"))
        if guide_c:
            parts.append(f"integrated-{guide_c}")
        v1c, v2c = _clean(row.get("Virus 1 Type")), _clean(row.get("Virus 2 Type"))
        if v1c:
            parts.append(v1c.replace(" ", "-"))
        if v2c:
            parts.append(v2c.replace(" ", "-"))
        label = " + ".join(parts)
        return label

    def calibration_trace_labels(wells: list) -> dict:
        base = {w: _condition_label(w) for w in wells}
        totals = {}
        for w in wells:
            totals[base[w]] = totals.get(base[w], 0) + 1
        seen = {}
        out = {}
        for w in wells:
            b = base[w]
            if totals[b] > 1:
                seen[b] = seen.get(b, 0) + 1
                out[w] = f"{b} (rep{seen[b]})"
            else:
                out[w] = b
        return out

    return (
        MIN_CELLS_PER_GATE,
        apply_markers,
        biexp,
        biexp_ticks,
        calibration_trace_labels,
        debris_gate,
        gated_range,
        interactive_hist,
        is_control_label,
        pctile_gate,
        replicate_labels,
        scatter_gate,
        tidy_summary_table,
    )


@app.cell(hide_code=True)
def _(
    MIN_CELLS_PER_GATE,
    apply_markers,
    biexp,
    biexp_ticks,
    debris_gate,
    gated_range,
    go,
    interactive_hist,
    is_control_label,
    mo,
    np,
    pd,
    replicate_labels,
    sample_sheet,
    scatter_gate,
    tidy_summary_table,
):
    def build_arm(
        *,
        title: str,
        description: str,
        groups: list,
        readout_channel: str,
        naive_well: str,
        control_group_index: int = 0,
        caveat: str = None,
        summary_flag: str = "",
    ):
        """groups: list of (group_label, wells, marker_specs) -- marker_specs is a
        list of (channel, threshold) applied as an AND gate (empty list = no
        infection/guide gate)."""
        naive_vals = debris_gate(naive_well)[readout_channel].values

        per_well_rows = []
        gated_by_well = {}
        for label, wells, marker_specs in groups:
            for w in wells:
                d = apply_markers(w, marker_specs) if marker_specs else debris_gate(w)
                gated_by_well[w] = (label, d)

        control_label, control_wells, control_marker_specs = groups[control_group_index]
        control_vals = np.concatenate([
            gated_by_well[w][1][readout_channel].values for w in control_wells
        ]) if all(len(gated_by_well[w][1]) > 0 for w in control_wells) else np.array([])
        gate = float(np.percentile(control_vals, 1)) if len(control_vals) else None
        gate50 = float(np.percentile(control_vals, 50)) if len(control_vals) else None

        pooled_vals = np.concatenate([d[readout_channel].values for _, d in gated_by_well.values() if len(d)])
        xr = gated_range(naive_vals, pooled_vals, hi_pct=99.0)

        cofactor = (max(gate, 1.0) / 5.0) if gate else (max(float(np.median(np.abs(pooled_vals))), 1.0) / 5.0)

        def _tx(v):
            return biexp(v, cofactor)

        naive_vals_t = _tx(naive_vals)
        xr_t = (float(_tx(xr[0])), float(_tx(xr[1])))
        gate_t = float(_tx(gate)) if gate is not None else None

        # One clickable histogram trace PER REPLICATE WELL (never pooled across
        # replicates within a condition), matching the scatter's existing
        # per-well convention.
        hist_traces = {}
        scatter_traces = {}
        for w, (label, d) in gated_by_well.items():
            if len(d):
                vals_t = _tx(d[readout_channel].values)
                hist_traces[f"{label}: {w}"] = vals_t
                scatter_traces[f"{label}: {w}"] = (vals_t, d["SSC-A"].values)
        hist_traces = {"unstained (naive)": naive_vals_t, **hist_traces}

        gate50_t = float(_tx(gate50)) if gate50 is not None else None
        _gate1_label = f"1st %ile mock ({control_label}) = {gate:,.0f}" if gate else ""
        _gate50_label = f"50th %ile mock ({control_label}) = {gate50:,.0f}" if gate50 else ""

        fig_hist = interactive_hist(
            hist_traces, gate_t, _gate1_label,
            xr_t, f"{title}: {readout_channel} distribution (biexponential x-axis, cofactor={cofactor:,.0f})",
            gate2=gate50_t, gate2_label=_gate50_label, gate2_color="purple",
        )
        fig_scatter = scatter_gate(
            scatter_traces, xr_t, f"{title}: {readout_channel} vs SSC-A (biexponential x-axis)", y_chan="SSC-A",
            gate=gate_t, gate_label=_gate1_label, gate2=gate50_t, gate2_label=_gate50_label, gate2_color="purple",
        )

        _tickvals, _ticktext = biexp_ticks(xr[0], xr[1], cofactor)
        fig_hist.update_layout(xaxis=dict(
            tickvals=_tickvals, ticktext=_ticktext, title=f"{readout_channel} (biexponential scale, original units)",
        ))
        fig_scatter.update_layout(xaxis=dict(
            tickvals=_tickvals, ticktext=_ticktext, title=f"{readout_channel} (biexponential scale, original units)",
        ))

        excluded = []
        for label, wells, _ in groups:
            for w in wells:
                d = gated_by_well[w][1]
                n = len(d)
                if n < MIN_CELLS_PER_GATE and not is_control_label(label):
                    excluded.append((w, label, n))
                    continue
                pct = 100 * np.mean(d[readout_channel].values < gate) if gate is not None and n else np.nan
                pct50 = 100 * np.mean(d[readout_channel].values < gate50) if gate50 is not None and n else np.nan
                per_well_rows.append({
                    "well": w, "label": label, "n": n,
                    "pct_below_gate": pct, "pct_below_gate50": pct50,
                })

        df = pd.DataFrame(per_well_rows)
        if len(df) and gate is not None:
            ctrl_mean = df.loc[df["label"] == control_label, "pct_below_gate"].mean()
            df["knockdown_pp"] = df["pct_below_gate"] - ctrl_mean
        if len(df) and gate50 is not None:
            ctrl_mean50 = df.loc[df["label"] == control_label, "pct_below_gate50"].mean()
            df["knockdown_pp_50"] = df["pct_below_gate50"] - ctrl_mean50

        bar_fig = go.Figure()
        bar_fig_50 = go.Figure()
        if len(df):
            rl = {}
            for label, wells, _ in groups:
                rl.update(replicate_labels([w for w in wells if w in df["well"].values]))
            x_labels = [rl.get(w, w) for w in df["well"]]
            colors_map = {g[0]: c for g, c in zip(groups, ["#4C78A8", "#E45756", "#54A24B", "#F58518"])}
            bar_fig.add_trace(go.Bar(
                x=x_labels, y=df["pct_below_gate"],
                text=[f"{v:.1f}%<br>(n={n:,})" for v, n in zip(df["pct_below_gate"], df["n"])],
                textposition="outside",
                marker_color=[colors_map.get(l, "#999") for l in df["label"]],
            ))
            ymax = max(50, float(df["pct_below_gate"].max()) * 1.25) if df["pct_below_gate"].notna().any() else 100
            bar_fig.update_layout(
                title=f"{title}: % below 1st-percentile-of-mock gate (per replicate; tail effect)",
                yaxis_title="% below gate", yaxis_range=[0, min(100, ymax)],
                height=380, margin=dict(t=60),
            )
            bar_fig_50.add_trace(go.Bar(
                x=x_labels, y=df["pct_below_gate50"],
                text=[f"{v:.1f}%<br>(n={n:,})" for v, n in zip(df["pct_below_gate50"], df["n"])],
                textposition="outside",
                marker_color=[colors_map.get(l, "#999") for l in df["label"]],
            ))
            ymax50 = max(50, float(df["pct_below_gate50"].max()) * 1.25) if df["pct_below_gate50"].notna().any() else 100
            bar_fig_50.update_layout(
                title=f"{title}: % below 50th-percentile-of-mock gate (per replicate; population-level shift)",
                yaxis_title="% below gate", yaxis_range=[0, min(100, ymax50)],
                height=380, margin=dict(t=60),
            )

        excl_note = ""
        if excluded:
            items = "; ".join(f"{lbl} ({w}), n={n}" for w, lbl, n in excluded)
            excl_note = f"*Excluded (population < {MIN_CELLS_PER_GATE} cells):* {items}"

        meta_cols = ["Well ID", "Sample ID #", "Cell Line", "Integrated Guide",
                     "Doxycycline Induction", "Virus 1 Type", "Virus 2 Type"]
        all_wells = [w for _, wells, _ in groups for w in wells]
        meta_table = mo.ui.table(sample_sheet.loc[all_wells, meta_cols], selection=None)

        _metric_note = mo.md(
            "*Two complementary metrics below: the red **1st-percentile-of-mock** "
            "gate captures a tail/low-signal effect (cells shifted furthest down); "
            "the purple **50th-percentile-of-mock** gate captures a broader "
            "population-level shift (whether the bulk of the distribution moved, "
            "not just the tail).*"
        )
        blocks = [mo.md(description)]
        if caveat:
            blocks.append(mo.callout(mo.md(caveat), kind="warn"))
        blocks.append(mo.ui.tabs({"Well metadata": meta_table}))
        blocks.append(_metric_note)
        blocks.append(fig_hist)
        blocks.append(fig_scatter)
        blocks.append(bar_fig)
        blocks.append(bar_fig_50)
        if excl_note:
            blocks.append(mo.md(excl_note))
        if len(df):
            blocks.append(tidy_summary_table(df.sort_values(["label", "well"])))

        primary_label, primary_wells, _ = groups[-1]
        _cell_line = sample_sheet.loc[primary_wells[0], "Cell Line"] if primary_wells else ""
        summary = {
            "arm": title,
            "cell_line": _cell_line,
            "mock_group": control_label,
            "guide_group": primary_label,
            "knockdown_pp_1st": float(df.loc[df["label"] == primary_label, "knockdown_pp"].mean()) if len(df) and "knockdown_pp" in df else float("nan"),
            "knockdown_pp_50th": float(df.loc[df["label"] == primary_label, "knockdown_pp_50"].mean()) if len(df) and "knockdown_pp_50" in df else float("nan"),
            "flag": summary_flag,
        }
        return mo.vstack(blocks), summary

    return (build_arm,)


@app.cell(hide_code=True)
def _(build_arm, infection_gate_widgets):
    arm1_content_base, arm1_summary = build_arm(
        title="WTC11 transient B2M-GFP virus",
        description=(
            "**Construct:** `NT-GFP` (non-targeting guide, GFP-marked) vs `B2M-GFP` "
            "(B2M-targeting guide, GFP-marked) delivered as a single transient virus "
            "into WTC11 cells (no separate effector virus in this arm -- WTC11 "
            "presumably carries the CRISPRi effector machinery already). "
            "GFP+ = confirmed transduction. Readout: B2M-APC."
        ),
        groups=[
            ("NT-GFP (control)", ["B1", "B2"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
            ("B2M-GFP (guide)", ["B3", "B4", "B5", "B6"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
        ],
        readout_channel="APC-A",
        naive_well="A3",
    )
    arm1_content_base
    return arm1_content_base, arm1_summary


@app.cell(hide_code=True)
def _(build_arm, infection_gate_widgets):
    arm2_content_base, arm2_summary = build_arm(
        title="17_3 stably-integrated guide + AA239 effector (two-component)",
        description=(
            "**Construct:** the guide (`NT-GFP` control or `B2M-GFP` B2M-targeting) "
            "is stably integrated into 17_3 cells; `AA239 dCas9-KRAB Thy1.1` is "
            "delivered as a separate transient effector virus (Thy1.1-APC+ = "
            "confirmed effector delivery, detected via antibody stain since no "
            "fluorescent Thy1.1 marker is otherwise used in these wells). Functional "
            "knockdown requires both the integrated guide AND effector delivery "
            "(GFP+ AND Thy1.1-APC+). Readout: B2M-PE-Dazzle-594 (read on the "
            "B610-ECD-A channel -- verified empirically against the B2M-PE-Dazzle "
            "single-stain control well, E9)."
        ),
        groups=[
            ("NT-GFP (control)", ["C1", "C2", "C5", "C6", "D1", "D2"],
             [("FITC-A", infection_gate_widgets.value["FITC-A__A7"]), ("APC-A", infection_gate_widgets.value["APC-A__D11_fmo"])]),
            ("B2M-GFP (guide)", ["C3", "C4", "C7", "C8", "D3", "D4"],
             [("FITC-A", infection_gate_widgets.value["FITC-A__A7"]), ("APC-A", infection_gate_widgets.value["APC-A__D11_fmo"])]),
        ],
        readout_channel="B610-ECD-A",
        naive_well="A7",
    )
    arm2_content_base
    return arm2_content_base, arm2_summary


@app.cell(hide_code=True)
def _(build_arm, infection_gate_widgets):
    arm3_content_base, arm3_summary = build_arm(
        title="17_3 fully-transient two-virus (guide + AA239 effector)",
        description=(
            "**Construct:** same target (B2M) and same `AA239 dCas9-KRAB Thy1.1` "
            "effector virus as the integrated-guide arm above, but here the guide "
            "itself (`NT-GFP`/`B2M-GFP`) is ALSO delivered as a transient virus "
            "rather than being stably integrated -- a fully-transient two-virus "
            "delivery system, kept as its own arm rather than pooled with the "
            "integrated-guide arm since guide-delivery mode is a distinct biological "
            "variable. GFP+ AND Thy1.1-APC+ double-positive gate, same idiom as the "
            "arm above. Readout: B2M-PE-Dazzle-594 (B610-ECD-A)."
        ),
        groups=[
            ("NT-GFP (control)", ["D7", "D8"],
             [("FITC-A", infection_gate_widgets.value["FITC-A__A7"]), ("APC-A", infection_gate_widgets.value["APC-A__D11_fmo"])]),
            ("B2M-GFP (guide)", ["D9", "D10"],
             [("FITC-A", infection_gate_widgets.value["FITC-A__A7"]), ("APC-A", infection_gate_widgets.value["APC-A__D11_fmo"])]),
        ],
        readout_channel="B610-ECD-A",
        naive_well="A7",
        summary_flag="Only n=2 wells per group -- smaller replicate count than arm 2's n=6.",
    )
    arm3_content_base
    return arm3_content_base, arm3_summary


@app.cell(hide_code=True)
def _(
    MIN_CELLS_PER_GATE,
    apply_markers,
    infection_gate_widgets,
    is_control_label,
    np,
    pd,
):
    # Marker-stratified knockdown: same 1st/50th-percentile-of-mock metric as every
    # arm above, but computed on top of (AND'd with) each arm's existing
    # guide/effector-delivery gate -- never in place of it. Only ASGR1 (BV421-A,
    # co-stained in-well) is available for stratification in this dataset -- no
    # Albumin co-stain exists in any arm's own wells (see the informational-only
    # hepatocyte-marker section above).
    def stratum_filter(d, asgr1_spec):
        if asgr1_spec is not None:
            d = d[d[asgr1_spec[0]].values > asgr1_spec[1]]
        return d

    def stratified_knockdown_rows(arm_label, groups, readout_channel, control_group_index, asgr1_gate=None):
        strata = [("All cells (unstratified)", None)]
        if asgr1_gate is not None:
            strata.append(("ASGR1+", ("BV421-A", asgr1_gate)))

        control_label, control_wells, _ = groups[control_group_index]
        rows = []
        for stratum_label, asgr1_spec in strata:
            infection_gated = {}
            for label, wells, marker_specs in groups:
                for w in wells:
                    d = apply_markers(w, marker_specs)
                    d = stratum_filter(d, asgr1_spec)
                    infection_gated[w] = (label, d)

            control_vals = np.concatenate([infection_gated[w][1][readout_channel].values for w in control_wells]) \
                if all(len(infection_gated[w][1]) > 0 for w in control_wells) else np.array([])
            gate1 = float(np.percentile(control_vals, 1)) if len(control_vals) else None
            gate50 = float(np.percentile(control_vals, 50)) if len(control_vals) else None

            stratum_rows = []
            for label, wells, _ in groups:
                for w in wells:
                    d = infection_gated[w][1]
                    n = len(d)
                    if n < MIN_CELLS_PER_GATE and not is_control_label(label):
                        stratum_rows.append({"arm": arm_label, "stratum": stratum_label, "well": w, "label": label, "n": n,
                                              "pct_below_1st": np.nan, "pct_below_50th": np.nan, "excluded_low_n": True})
                        continue
                    p1 = 100 * np.mean(d[readout_channel].values < gate1) if gate1 is not None and n else np.nan
                    p50 = 100 * np.mean(d[readout_channel].values < gate50) if gate50 is not None and n else np.nan
                    stratum_rows.append({"arm": arm_label, "stratum": stratum_label, "well": w, "label": label, "n": n,
                                          "pct_below_1st": p1, "pct_below_50th": p50, "excluded_low_n": False})

            _sdf = pd.DataFrame(stratum_rows)
            ctrl_mean_1 = _sdf.loc[_sdf["label"] == control_label, "pct_below_1st"].mean()
            ctrl_mean_50 = _sdf.loc[_sdf["label"] == control_label, "pct_below_50th"].mean()
            _sdf["knockdown_pp_1st"] = _sdf["pct_below_1st"] - ctrl_mean_1
            _sdf["knockdown_pp_50th"] = _sdf["pct_below_50th"] - ctrl_mean_50
            rows.extend(_sdf.to_dict("records"))
        return rows

    STRAT_ARM_CFGS = [
        {
            "arm_label": "WTC11 transient B2M-GFP virus", "naive_well": "A3",
            "groups": [
                ("NT-GFP (control)", ["B1", "B2"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
                ("B2M-GFP (guide)", ["B3", "B4", "B5", "B6"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
            ],
            "readout_channel": "APC-A", "control_group_index": 0,
            "asgr1_gate": infection_gate_widgets.value["BV421-A__A3"],
        },
        {
            "arm_label": "17_3 stably-integrated guide + AA239 effector (two-component)", "naive_well": "A7",
            "groups": [
                ("NT-GFP (control)", ["C1", "C2", "C5", "C6", "D1", "D2"],
                 [("FITC-A", infection_gate_widgets.value["FITC-A__A7"]), ("APC-A", infection_gate_widgets.value["APC-A__D11_fmo"])]),
                ("B2M-GFP (guide)", ["C3", "C4", "C7", "C8", "D3", "D4"],
                 [("FITC-A", infection_gate_widgets.value["FITC-A__A7"]), ("APC-A", infection_gate_widgets.value["APC-A__D11_fmo"])]),
            ],
            "readout_channel": "B610-ECD-A", "control_group_index": 0,
            "asgr1_gate": infection_gate_widgets.value["BV421-A__A7"],
        },
        {
            "arm_label": "17_3 fully-transient two-virus (guide + AA239 effector)", "naive_well": "A7",
            "groups": [
                ("NT-GFP (control)", ["D7", "D8"],
                 [("FITC-A", infection_gate_widgets.value["FITC-A__A7"]), ("APC-A", infection_gate_widgets.value["APC-A__D11_fmo"])]),
                ("B2M-GFP (guide)", ["D9", "D10"],
                 [("FITC-A", infection_gate_widgets.value["FITC-A__A7"]), ("APC-A", infection_gate_widgets.value["APC-A__D11_fmo"])]),
            ],
            "readout_channel": "B610-ECD-A", "control_group_index": 0,
            "asgr1_gate": infection_gate_widgets.value["BV421-A__A7"],
        },
    ]

    _strat_all_rows = []
    for _cfg in STRAT_ARM_CFGS:
        _strat_all_rows.extend(stratified_knockdown_rows(
            _cfg["arm_label"], _cfg["groups"], _cfg["readout_channel"], _cfg["control_group_index"],
            asgr1_gate=_cfg["asgr1_gate"],
        ))

    marker_strat_df = pd.DataFrame(_strat_all_rows)
    marker_strat_df
    return STRAT_ARM_CFGS, marker_strat_df, stratum_filter


@app.cell(hide_code=True)
def _(
    STRAT_ARM_CFGS,
    apply_markers,
    biexp,
    debris_gate,
    gated_range,
    go,
    interactive_hist,
    marker_strat_df,
    mo,
    np,
    pd,
    replicate_labels,
    scatter_gate,
    stratum_filter,
):
    def _tidy_strat_table(df):
        _d = df.copy()
        for c in ["pct_below_1st", "pct_below_50th", "knockdown_pp_1st", "knockdown_pp_50th"]:
            _d[c] = _d[c].round(1)
        _d = _d.drop(columns=["excluded_low_n"])
        _d = _d.rename(columns={
            "stratum": "Stratum", "well": "Well", "label": "Condition", "n": "n (gated cells)",
            "pct_below_1st": "% below 1st-pctile-mock", "pct_below_50th": "% below 50th-pctile-mock",
            "knockdown_pp_1st": "Knockdown, 1st pctile (pp)", "knockdown_pp_50th": "Knockdown, 50th pctile (pp)",
        })
        return mo.ui.table(_d, selection=None)

    def _strat_bar_chart(arm_label, stratum_label, metric_col, title_suffix):
        _sub = marker_strat_df[(marker_strat_df["arm"] == arm_label) & (marker_strat_df["stratum"] == stratum_label)]
        if not len(_sub):
            return None
        groups = next(c for c in STRAT_ARM_CFGS if c["arm_label"] == arm_label)["groups"]
        rl = {}
        for label, wells, _ in groups:
            rl.update(replicate_labels([w for w in wells if w in _sub["well"].values]))
        colors_map = {g[0]: c for g, c in zip(groups, ["#4C78A8", "#E45756", "#54A24B", "#F58518"])}
        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=[rl.get(w, w) for w in _sub["well"]], y=_sub[metric_col],
            text=[f"{v:.1f}%<br>(n={n:,})" for v, n in zip(_sub[metric_col], _sub["n"])],
            textposition="outside",
            marker_color=[colors_map.get(l, "#999") for l in _sub["label"]],
        ))
        ymax = max(50, float(_sub[metric_col].max()) * 1.25) if _sub[metric_col].notna().any() else 100
        fig.update_layout(
            title=f"{arm_label} -- {stratum_label}: {title_suffix} (per replicate)",
            yaxis_title="% below gate", yaxis_range=[0, min(100, ymax)], height=360, margin=dict(t=60),
        )
        return fig

    def build_marker_stratum_block(cfg, stratum_label, asgr1_spec):
        groups = cfg["groups"]
        readout_channel = cfg["readout_channel"]
        control_group_index = cfg["control_group_index"]
        naive_well = cfg["naive_well"]
        control_label, control_wells, _ = groups[control_group_index]
        naive_vals_raw = debris_gate(naive_well)[readout_channel].values

        gated_by_well = {}
        for label, wells, marker_specs in groups:
            for w in wells:
                d = apply_markers(w, marker_specs)
                d = stratum_filter(d, asgr1_spec)
                gated_by_well[w] = (label, d)

        control_vals = np.concatenate([gated_by_well[w][1][readout_channel].values for w in control_wells]) \
            if all(len(gated_by_well[w][1]) > 0 for w in control_wells) else np.array([])
        if not len(control_vals):
            return mo.md(f"**{stratum_label}:** too few control cells to draw a gate -- skipped.")
        gate1 = float(np.percentile(control_vals, 1))
        gate50 = float(np.percentile(control_vals, 50))

        pooled_vals = np.concatenate([d[readout_channel].values for _, d in gated_by_well.values() if len(d)])
        xr = gated_range(naive_vals_raw, pooled_vals, hi_pct=99.0)
        cofactor = max(gate1, 1.0) / 5.0

        def _tx(v):
            return biexp(v, cofactor)

        xr_t = (float(_tx(xr[0])), float(_tx(xr[1])))
        gate1_t, gate50_t = float(_tx(gate1)), float(_tx(gate50))

        hist_traces = {"unstained (naive)": _tx(naive_vals_raw)}
        scatter_traces = {}
        for w, (label, d) in gated_by_well.items():
            if len(d):
                vals_t = _tx(d[readout_channel].values)
                hist_traces[f"{label}: {w}"] = vals_t
                scatter_traces[f"{label}: {w}"] = (vals_t, d["SSC-A"].values)

        fig_hist = interactive_hist(
            hist_traces, gate1_t, f"1st %ile mock ({control_label}) = {gate1:,.0f}",
            xr_t, f"{cfg['arm_label']} -- {stratum_label}: {readout_channel} distribution",
            gate2=gate50_t, gate2_label=f"50th %ile mock ({control_label}) = {gate50:,.0f}", gate2_color="purple",
        )
        fig_scatter = scatter_gate(
            scatter_traces, xr_t, f"{cfg['arm_label']} -- {stratum_label}: {readout_channel} vs SSC-A",
            y_chan="SSC-A", gate=gate1_t, gate_label="1st %ile mock", gate2=gate50_t, gate2_label="50th %ile mock", gate2_color="purple",
        )
        bar1 = _strat_bar_chart(cfg["arm_label"], stratum_label, "pct_below_1st", "% below 1st-pctile-mock gate")
        bar50 = _strat_bar_chart(cfg["arm_label"], stratum_label, "pct_below_50th", "% below 50th-pctile-mock gate")

        blocks = [fig_hist, fig_scatter]
        if bar1 is not None:
            blocks.append(bar1)
        if bar50 is not None:
            blocks.append(bar50)
        blocks.append(_tidy_strat_table(marker_strat_df[
            (marker_strat_df["arm"] == cfg["arm_label"]) & (marker_strat_df["stratum"] == stratum_label)
        ]))
        return mo.vstack(blocks)

    def _arm_strat_comparison(cfg):
        primary_label = cfg["groups"][-1][0]
        _sub = marker_strat_df[(marker_strat_df["arm"] == cfg["arm_label"]) & (marker_strat_df["label"] == primary_label)]
        agg = _sub.groupby("stratum", as_index=False)[["knockdown_pp_1st", "knockdown_pp_50th"]].mean()
        _order = ["All cells (unstratified)", "ASGR1+"]
        agg["stratum"] = pd.Categorical(agg["stratum"], categories=_order, ordered=True)
        agg = agg.sort_values("stratum")

        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=agg["stratum"].astype(str), y=agg["knockdown_pp_1st"], name="1st pctile (tail effect)",
            marker_color="red", text=[f"{v:.1f}pp" for v in agg["knockdown_pp_1st"]], textposition="outside",
        ))
        fig.add_trace(go.Bar(
            x=agg["stratum"].astype(str), y=agg["knockdown_pp_50th"], name="50th pctile (population shift)",
            marker_color="purple", text=[f"{v:.1f}pp" for v in agg["knockdown_pp_50th"]], textposition="outside",
        ))
        _all_vals = list(agg["knockdown_pp_1st"]) + list(agg["knockdown_pp_50th"])
        _all_vals = [v for v in _all_vals if v == v]
        ymax = max(50, max(_all_vals) * 1.25) if _all_vals else 50
        ymin = min(0, min(_all_vals) * 1.25) if _all_vals else 0
        fig.update_layout(
            title=f"{cfg['arm_label']}: knockdown across strata ({primary_label} vs. mock)",
            barmode="group", yaxis_title="Knockdown (pp)", yaxis_range=[ymin, ymax],
            height=380, margin=dict(t=60),
            legend=dict(itemclick="toggle", itemdoubleclick="toggleothers"),
        )
        return fig, agg

    def _arm_strat_tldr(cfg, agg):
        _base = agg[agg["stratum"] == "All cells (unstratified)"]
        if not len(_base):
            return mo.md("**TL;DR:** insufficient data to compare strata.")
        base1 = float(_base["knockdown_pp_1st"].iloc[0])
        base50 = float(_base["knockdown_pp_50th"].iloc[0])
        row = agg[agg["stratum"] == "ASGR1+"]
        if not len(row):
            return mo.md(f"**TL;DR:** unstratified knockdown is {base1:.1f}pp (1st pctile) / {base50:.1f}pp (50th pctile); no ASGR1+ stratum available.")
        v1 = float(row["knockdown_pp_1st"].iloc[0])
        v50 = float(row["knockdown_pp_50th"].iloc[0])
        d1, d50 = v1 - base1, v50 - base50
        dir1 = "lowers" if d1 < -0.5 else "raises" if d1 > 0.5 else "leaves essentially unchanged"
        dir50 = "lowers" if d50 < -0.5 else "raises" if d50 > 0.5 else "leaves essentially unchanged"
        return mo.md(
            f"**TL;DR:** **ASGR1+** selection {dir1} 1st-percentile knockdown ({base1:.1f}pp -> {v1:.1f}pp, {d1:+.1f}pp) "
            f"and {dir50} the 50th-percentile metric ({base50:.1f}pp -> {v50:.1f}pp, {d50:+.1f}pp)."
        )

    def build_arm_with_strat(cfg, base_content):
        asgr1_gate = cfg["asgr1_gate"]
        tabs = {"All cells": base_content}
        if asgr1_gate is not None:
            tabs["ASGR1+"] = build_marker_stratum_block(cfg, "ASGR1+", ("BV421-A", asgr1_gate))
        if len(tabs) == 1:
            return base_content

        comparison_fig, agg = _arm_strat_comparison(cfg)
        tldr = _arm_strat_tldr(cfg, agg)

        return mo.vstack([
            comparison_fig,
            tldr,
            mo.md(
                "*ASGR1-BV421 (co-stained in-well) is the only hepatocyte marker "
                "available for stratification in this dataset -- Albumin was only "
                "stained on separate, non-arm reference wells (see the informational "
                "hepatocyte-marker section above). This dataset does not contain the "
                "previously-flagged non-specific ASGR1-FITC clone at all.*"
            ),
            mo.ui.tabs(tabs),
        ])

    return (build_arm_with_strat,)


@app.cell(hide_code=True)
def _(
    STRAT_ARM_CFGS,
    arm1_content_base,
    arm2_content_base,
    arm2_reagent_section,
    arm3_content_base,
    build_arm_with_strat,
    mo,
):
    _cfg_by_label = {c["arm_label"]: c for c in STRAT_ARM_CFGS}
    arm1_content = build_arm_with_strat(_cfg_by_label["WTC11 transient B2M-GFP virus"], arm1_content_base)
    arm2_content = mo.vstack([
        arm2_reagent_section,
        build_arm_with_strat(_cfg_by_label["17_3 stably-integrated guide + AA239 effector (two-component)"], arm2_content_base),
    ])
    arm3_content = build_arm_with_strat(_cfg_by_label["17_3 fully-transient two-virus (guide + AA239 effector)"], arm3_content_base)
    mo.md("Marker-stratified content folded into all 3 arms; Arm 2 also gets a transduction-reagent infection-gating chart at the top of its tab.")
    return arm1_content, arm2_content, arm3_content


if __name__ == "__main__":
    app.run()
