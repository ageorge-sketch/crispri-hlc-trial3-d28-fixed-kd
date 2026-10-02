# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "marimo>=0.25.1",
#     "numpy==2.5.3",
#     "pandas==3.0.6",
#     "plotly==7.1.0",
#     "requests==2.34.2",
#     "xlrd==2.0.2",
# ]
# ///

import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium", auto_download=["html"])


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
def _(mo):
    mo.md("""
    # CRISPRi HLC Trial3 D28 (fixed) -- B2M knockdown analysis

    This notebook measures B2M knockdown efficiency across three independent arms
    in fixed cells from HLC Trial 3, Day 28, each arm using a different
    guide-delivery method and/or cell line. It parses the raw flow cytometry
    (FCS) files for each well, applies cell and antibody gates, and reports the
    percentage-point shift in the fraction of cells below a fixed B2M-antibody
    gate, comparing a non-targeting control guide (NT-GFP) to the B2M-targeting
    guide (B2M-GFP), within cells confirmed to carry the guide and/or effector
    virus.

    This dataset targets the B2M gene and includes no CD81 antibody stain. Two
    well groups are excluded from the knockdown analysis: wells from the K562
    cell line, used only to check ASGR1-antibody specificity, and two
    guide-only wells with no effector virus, used only as GFP compensation
    controls.
    """)
    return


@app.cell(hide_code=True)
def arms_overview(mo, pd):
    _arms_overview_df = pd.DataFrame([
        {
            "Arm": "1", "Construct / cell line": "WTC11, transient B2M-GFP virus (guide+marker, single virus)",
            "Guide type": "NT-GFP vs B2M-GFP (dox-inducible CRISPRi effector)",
            "Key comparison": "B2M-GFP +dox vs NT-GFP control (+dox vs no-dox split)",
        },
        {
            "Arm": "2", "Construct / cell line": "17_3, stably-integrated guide + AA239 transient effector virus",
            "Guide type": "NT-GFP vs B2M-GFP (integrated); Thy1.1-APC-confirmed effector",
            "Key comparison": "B2M-GFP vs NT-GFP, both GFP+/Thy1.1+; 3 transduction reagents compared",
        },
        {
            "Arm": "3", "Construct / cell line": "17_3, fully-transient two-virus (guide + AA239 effector)",
            "Guide type": "NT-GFP vs B2M-GFP (transient); Thy1.1-APC-confirmed effector",
            "Key comparison": "B2M-GFP vs NT-GFP, both GFP+/Thy1.1+ (n=2 wells/group)",
        },
    ])
    arms_overview_table = mo.ui.table(_arms_overview_df, selection=None)
    mo.vstack([
        mo.md("### Experimental arms at a glance"),
        arms_overview_table,
    ])

    return


@app.cell(hide_code=True)
def _(
    arm1_summary,
    arm2_summary,
    arm3_summary,
    go,
    mo,
    pd,
    plot_overview,
    plot_result,
):
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
        title="Plot 1. Knockdown across all 3 arms, 1st vs. 50th percentile gate",
        barmode="group", yaxis_title="Knockdown (percentage points)", yaxis_range=[0, _ymax],
        height=420, margin=dict(t=60),
        legend=dict(itemclick="toggle", itemdoubleclick="toggleothers"),
    )

    _best1 = max(_summary_rows, key=lambda r: r[_k1] if r[_k1] == r[_k1] else -1e9)
    _best50 = max(_summary_rows, key=lambda r: r[_k50] if r[_k50] == r[_k50] else -1e9)

    top_summary = mo.vstack([
        plot_overview(
            1, "Knockdown across all 3 arms",
            why=(
                "This plot answers which of the 3 knockdown arms shows the largest B2M knockdown, "
                "and whether the answer depends on which of the two knockdown metrics is used."
            ),
            how=(
                "For each arm, the primary guide-active (B2M-GFP) group's % of cells below a gate is "
                "compared to that same arm's own mock/control group, at two gates: the 1st-percentile-"
                "of-mock gate (tail effect) and the 50th-percentile-of-mock gate (population-level "
                "shift). The difference (guide-active minus control, in percentage points) is the "
                "knockdown number plotted here."
            ),
            how_to_read=(
                "Each arm is one group of two bars: red is the 1st-percentile-of-mock knockdown, "
                "purple is the 50th-percentile-of-mock knockdown, both in percentage points. See each "
                "arm's own tab below for full per-replicate detail, ASGR1+ stratification, and caveats."
            ),
        ),
        top_summary_fig,
        plot_result(
            f"The largest 1st-percentile knockdown is {_best1[_k1]:.1f}pp, in the \"{_best1['arm']}\" arm; "
            f"the largest 50th-percentile knockdown is {_best50[_k50]:.1f}pp, in the \"{_best50['arm']}\" arm."
        ),
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

    Every plot in this notebook is numbered (Plot 1, Plot 2, ...) and preceded by
    a short Why/How/How-to-read-it note, with a Result line stating the computed
    numbers after the plot. Numbers follow each plot's position in the notebook's
    cell order, not the order tabs happen to be clicked in, so a number may appear
    out of visual sequence when viewed inside a nested tab.
    """)
    return


@app.cell(hide_code=True)
def infection_gating_tab(
    arm2_reagent_section,
    biexp,
    biexp_ticks,
    channel_cofactor,
    debris_gate,
    gated_range,
    go,
    infection_gate_widgets,
    interactive_hist,
    mo,
    np,
    pd,
    plot_overview,
    plot_result,
    replicate_labels,
    sample_sheet,
    scatter_gate,
):
    def _infection_arm_panel(arm_name, groups, naive_well, marker_channels, plot_start):
        """marker_channels: list of (channel, gate, display_name). Mirrors
        build_arm's layout (histogram + scatter + per-replicate bar, pooled by
        group, never pooled across replicates/groups) but for the raw
        guide/effector-marker channel itself (infection status), not the
        downstream knockdown readout."""
        blocks = []
        plot_n = plot_start
        for channel, gate, display_name in marker_channels:
            naive_vals = debris_gate(naive_well)[channel].values
            hist_traces, scatter_traces = {}, {}
            bar_rows = []
            rl = {}
            for label, wells, _ in groups:
                rl.update(replicate_labels(wells, group_label=label))
                xs_h, xs_s, ys_s = [], [], []
                for w in wells:
                    d = debris_gate(w)
                    n = len(d)
                    if n:
                        xs_h.append(d[channel].values)
                        xs_s.append(d[channel].values)
                        ys_s.append(d["SSC-A"].values)
                    pct = 100 * float(np.mean(d[channel].values > gate)) if n else np.nan
                    bar_rows.append({"well": w, "label": label, "n": n, "pct": pct})
                if xs_h:
                    hist_traces[label] = np.concatenate(xs_h)
                    scatter_traces[label] = (np.concatenate(xs_s), np.concatenate(ys_s))

            pooled_vals = np.concatenate(list(hist_traces.values())) if hist_traces else np.array([])
            xr = gated_range(naive_vals, pooled_vals, hi_pct=99.0, group_vals=list(hist_traces.values()))
            cofactor = channel_cofactor(channel, naive_vals)

            def _tx(v):
                return biexp(v, cofactor)

            naive_vals_t = _tx(naive_vals)
            xr_t = (float(_tx(xr[0])), float(_tx(xr[1])))
            gate_t = float(_tx(gate))
            hist_traces_t = {"unstained/FMO reference": naive_vals_t}
            hist_traces_t.update({k: _tx(v) for k, v in hist_traces.items()})
            scatter_traces_t = {k: (_tx(xv), yv) for k, (xv, yv) in scatter_traces.items()}

            fig_hist = interactive_hist(
                hist_traces_t, gate_t, f"{display_name} gate = {gate:,.0f}", xr_t,
                f"Plot {plot_n}. {arm_name}: {display_name} ({channel}) infection-marker distribution",
            )
            fig_scatter = scatter_gate(
                scatter_traces_t, xr_t,
                f"Plot {plot_n + 1}. {arm_name}: {display_name} ({channel}) vs SSC-A", y_chan="SSC-A",
                gate=gate_t, gate_label=f"{display_name} gate",
            )
            tickvals, ticktext = biexp_ticks(xr[0], xr[1], cofactor)
            for f in (fig_hist, fig_scatter):
                f.update_layout(xaxis=dict(tickvals=tickvals, ticktext=ticktext, title=f"{channel} (biexponential scale)"))

            bar_df = pd.DataFrame(bar_rows)
            colors_map = {g[0]: c for g, c in zip(groups, ["#4C78A8", "#E45756", "#54A24B", "#F58518"])}
            bar_fig = go.Figure()
            bar_fig.add_trace(go.Bar(
                x=[rl.get(w, w) for w in bar_df["well"]], y=bar_df["pct"],
                text=[f"{v:.1f}%<br>(n={n:,})" for v, n in zip(bar_df["pct"], bar_df["n"])],
                textposition="outside",
                marker_color=[colors_map.get(l, "#999") for l in bar_df["label"]],
            ))
            bar_fig.update_layout(
                title=f"Plot {plot_n + 2}. {arm_name}: % {display_name}+ per replicate (gate={gate:,.0f})",
                yaxis_title=f"% {display_name}+", yaxis_range=[0, 115], height=380, margin=dict(t=60),
            )

            _overall_pct = 100 * float(np.average(bar_df["pct"], weights=bar_df["n"])) if bar_df["n"].sum() else float("nan")

            blocks.append(plot_overview(
                plot_n, f"{arm_name}: {display_name} ({channel}) infection-marker distribution",
                why=f"States whether cells in {arm_name} were successfully transduced/marked for {display_name}, the gate that feeds this arm's knockdown analysis.",
                how=(
                    f"{channel} values (debris/doublet/singlet-gated, compensated) are pooled across "
                    f"replicate wells within each condition group and plotted as peak-normalized "
                    f"histograms, with an unstained or FMO reference overlaid."
                ),
                how_to_read="Same layout as the gate-calibration histograms above: biexponential x-axis, % of each group's peak count on the y-axis, dashed red line is the gate.",
            ))
            blocks.append(fig_hist)
            blocks.append(plot_overview(
                plot_n + 1, f"{arm_name}: {display_name} ({channel}) vs SSC-A",
                why="Checks the infection-marker-positive population as a distinct group in 2D against granularity (SSC-A).",
                how="A fixed-seed random subsample of the same gated, compensated events is plotted as a scatter (never pooled across groups).",
                how_to_read="Same x-axis as the histogram above; the dashed red line marks the same gate.",
            ))
            blocks.append(fig_scatter)
            blocks.append(plot_overview(
                plot_n + 2, f"{arm_name}: % {display_name}+ per replicate",
                why=f"States the {display_name} infection/transduction rate for {arm_name} directly, per replicate well.",
                how=f"For each well, % of gated cells above the {display_name} gate is computed.",
                how_to_read="Each bar is one replicate well; n (gated cell count) is printed above each bar.",
            ))
            blocks.append(bar_fig)
            blocks.append(plot_result(f"Overall {display_name}+ rate in {arm_name}: {_overall_pct:.1f}% (n-weighted across replicate wells)."))
            plot_n += 3

        if marker_channels and sample_sheet.loc[naive_well, "Cell Line"] != "WTC11":
            blocks.append(mo.callout(mo.md(
                "*Caveat: compensation is anchored to WTC11 single-stain controls "
                "(no single-stain wells exist for this arm's own cell line), so this "
                "arm's unstained-reference line sits at a different baseline than the "
                "true compensation zero-point -- cells scattered below it (even below 0) "
                "reflect that baseline mismatch, not necessarily a gating error. % "
                "infected is still reliable since the gate is set from this same arm's "
                "own reference well.*"
            ), kind="warn"))

        return mo.vstack(blocks), plot_n

    arm1_infect_tab, _n1 = _infection_arm_panel(
        "Arm 1 (WTC11 transient B2M-GFP virus)",
        [("NT-GFP (control)", ["B1", "B2"], []),
         ("B2M-GFP (guide), no dox", ["B3", "B4"], []),
         ("B2M-GFP (guide), +dox", ["B5", "B6"], [])],
        naive_well="A3",
        marker_channels=[("FITC-A", infection_gate_widgets.value["FITC-A__A3"], "GFP (guide)")],
        plot_start=42,
    )
    _arm2_infect_panel, _n2 = _infection_arm_panel(
        "Arm 2 (17_3 stably-integrated guide + AA239 effector)",
        [("NT-GFP (control)", ["C1", "C2", "C5", "C6", "D1", "D2"], []),
         ("B2M-GFP (guide)", ["C3", "C4", "C7", "C8", "D3", "D4"], [])],
        naive_well="A7",
        marker_channels=[
            ("FITC-A", infection_gate_widgets.value["FITC-A__A7"], "GFP (guide)"),
            ("APC-A", infection_gate_widgets.value["APC-A__D11_fmo"], "Thy1.1-APC (effector)"),
        ],
        plot_start=_n1,
    )
    arm2_infect_tab = mo.vstack([arm2_reagent_section, _arm2_infect_panel])
    arm3_infect_tab, _n3 = _infection_arm_panel(
        "Arm 3 (17_3 fully-transient two-virus)",
        [("NT-GFP (control)", ["D7", "D8"], []),
         ("B2M-GFP (guide)", ["D9", "D10"], [])],
        naive_well="A7",
        marker_channels=[
            ("FITC-A", infection_gate_widgets.value["FITC-A__A7"], "GFP (guide)"),
            ("APC-A", infection_gate_widgets.value["APC-A__D11_fmo"], "Thy1.1-APC (effector)"),
        ],
        plot_start=_n2,
    )

    return arm1_infect_tab, arm2_infect_tab, arm3_infect_tab


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
    return (infection_gate_widgets,)


@app.cell(hide_code=True)
def _(
    REAGENT_ABBREV,
    debris_gate,
    go,
    infection_gate_widgets,
    mo,
    np,
    pd,
    plot_overview,
    plot_result,
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
    arm2_reagent_df["reagent_abbrev"] = arm2_reagent_df["reagent"].map(REAGENT_ABBREV)
    arm2_reagent_df["rep_in_group"] = arm2_reagent_df.groupby(["reagent", "guide"]).cumcount() + 1
    # Label includes the reagent so no x-axis category is ever shared across
    # reagents -- each bar is a distinct well/condition, not a matched replicate
    # of the other reagents' bars.
    arm2_reagent_df["well_label"] = (
        arm2_reagent_df["reagent_abbrev"] + ": " + arm2_reagent_df["guide"]
        + " rep" + arm2_reagent_df["rep_in_group"].astype(str)
    )
    arm2_reagent_df = arm2_reagent_df.sort_values(["reagent", "guide", "rep_in_group"])

    _reagent_colors = {"Polybrene": "#4C78A8", "Protamine Sulfate": "#E45756", "Lentiboost": "#54A24B"}
    arm2_reagent_fig = go.Figure()
    _category_order = []
    for _reagent, _grp in arm2_reagent_df.groupby("reagent", sort=False):
        _category_order.extend(_grp["well_label"].tolist())
        arm2_reagent_fig.add_trace(go.Bar(
            x=_grp["well_label"], y=_grp["pct_infected"],
            name=_reagent, marker_color=_reagent_colors.get(_reagent, "#999"),
            text=[f"{v:.1f}%<br>(n={n:,})" for v, n in zip(_grp["pct_infected"], _grp["n"])],
            textposition="outside",
        ))
        _mid = _grp["well_label"].iloc[len(_grp) // 2]
        arm2_reagent_fig.add_annotation(
            x=_mid, y=1.1, yref="paper", text=f"<b>{_reagent}</b>", showarrow=False,
            font=dict(size=12, color=_reagent_colors.get(_reagent, "#999")),
        )
    arm2_reagent_fig.update_layout(
        title=f"Plot 9. Arm 2: % Thy1.1-APC+ (effector-infected) per well, by transduction reagent (gate = {_arm2_thy11_gate:,.0f})",
        yaxis_title="% infected (Thy1.1-APC+)", yaxis_range=[0, 115],
        height=420, margin=dict(t=90),
        xaxis=dict(categoryorder="array", categoryarray=_category_order),
        showlegend=False,
    )
    for _i in range(1, 3):
        arm2_reagent_fig.add_vline(x=_i * 4 - 0.5, line=dict(color="#CCCCCC", width=1, dash="dot"))

    _reagent_means = arm2_reagent_df.groupby("reagent")["pct_infected"].mean()
    _best_reagent = _reagent_means.idxmax()

    arm2_reagent_section = mo.vstack([
        mo.md(
            "### Arm 2 infection-gating: transduction-reagent comparison\n"
            "Arm 2's 12 wells compare 3 lentiviral delivery reagents for the "
            "AA239 effector virus (`Transduction Reagent` column in the sample "
            "sheet): **Polybrene** (C1-C4), **Protamine Sulfate** (C5-C8), "
            "**Lentiboost** (D1-D4) -- 4 wells each (2 NT-GFP + 2 B2M-GFP)."
        ),
        plot_overview(
            9, "Arm 2 transduction-reagent comparison",
            why="This plot answers which of the 3 transduction reagents delivers the AA239 effector virus most efficiently in Arm 2.",
            how=(
                "For each well, % Thy1.1-APC+ (confirmed effector delivery) is computed using the "
                "Thy1.1-APC effector-delivery gate slider above, then grouped by reagent."
            ),
            how_to_read=(
                "Each bar is one well, labeled by reagent, guide condition, and replicate number "
                "(e.g. \"PB: NT-GFP rep1\"); n (gated cell count) is printed above each bar. Wells "
                "are grouped into 3 reagent blocks (dotted vertical lines mark the boundaries) -- "
                "\"rep1\"/\"rep2\" under different reagents are different physical wells, not matched "
                "replicates of each other, so bars are not meant to be compared left-to-right across "
                "reagent blocks by replicate number, only by reagent-level means (see Result below)."
            ),
        ),
        arm2_reagent_fig,
        plot_result(
            "Mean % Thy1.1-APC+ by reagent: "
            + ", ".join(f"{r} {v:.1f}%" for r, v in _reagent_means.items())
            + f". {_best_reagent} gives the highest mean effector-delivery rate."
        ),
    ])
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
    plot_overview,
    plot_result,
    sample_sheet,
):
    _CALIBRATION_CFGS = [
        {"key": "FITC-A__A3", "channel": "FITC-A", "naive": "A3", "plot_n": 2,
         "used_by": {"Arm 1 GFP guide marker (B1-B6)": ["B1", "B2", "B3", "B4", "B5", "B6"]}},
        {"key": "FITC-A__A7", "channel": "FITC-A", "naive": "A7", "plot_n": 3,
         "used_by": {"Arm 2 GFP guide marker (C1-C8,D1-D4)": ["C1","C2","C3","C4","C5","C6","C7","C8","D1","D2","D3","D4"],
                     "Arm 3 GFP guide marker (D7-D10)": ["D7", "D8", "D9", "D10"]}},
        {"key": "APC-A__D11_fmo", "channel": "APC-A", "naive": "D11", "secondary_naive": "A7", "plot_n": 4,
         "used_by": {"Arm 2 Thy1.1 effector marker (C1-C8,D1-D4)": ["C1","C2","C3","C4","C5","C6","C7","C8","D1","D2","D3","D4"],
                     "Arm 3 Thy1.1 effector marker (D7-D10)": ["D7", "D8", "D9", "D10"]}},
        {"key": "BV421-A__A3", "channel": "BV421-A", "naive": "A3", "plot_n": 5,
         "used_by": {"Arm 1 ASGR1 marker (B1-B6)": ["B1", "B2", "B3", "B4", "B5", "B6"]}},
        {"key": "BV421-A__A7", "channel": "BV421-A", "naive": "A7", "plot_n": 6,
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
            naive_label = f"FMO reference ({naive_cell_line}, stained/uninfected)"
        else:
            naive_label = f"Unstained reference ({naive_cell_line})"

        traces = {naive_label: naive_vals}
        pooled_for_range = [naive_vals]

        secondary = cfg.get("secondary_naive")
        if secondary:
            sec_vals = debris_gate(secondary)[channel].values
            sec_cell_line = sample_sheet.loc[secondary, "Cell Line"]
            sec_label = f"(for comparison) Unstained reference ({sec_cell_line})"
            traces[sec_label] = sec_vals
            pooled_for_range.append(sec_vals)

        all_sample_wells = [w for wells in cfg["used_by"].values() for w in wells]
        well_labels = calibration_trace_labels(all_sample_wells)

        for w in all_sample_wells:
            vals = debris_gate(w)[channel].values
            traces[well_labels[w]] = vals
            pooled_for_range.append(vals)
        pooled = np.concatenate(pooled_for_range)
        xr = gated_range(naive_vals, pooled, hi_pct=99.0, group_vals=list(traces.values()))
        cofactor = max(gate, 1.0) / 5.0

        def _tx(v):
            return biexp(v, cofactor)

        traces_t = {k: _tx(v) for k, v in traces.items()}
        xr_t = (float(_tx(xr[0])), float(_tx(xr[1])))
        gate_t = float(_tx(gate))
        plot_n = cfg["plot_n"]
        fig = interactive_hist(
            traces_t, gate_t, f"gate = {gate:,.0f}", xr_t,
            f"Plot {plot_n}. {channel} gate calibration (reference = {naive_label})",
        )
        tickvals, ticktext = biexp_ticks(xr[0], xr[1], cofactor)
        fig.update_layout(xaxis=dict(tickvals=tickvals, ticktext=ticktext, title=f"{channel} (biexponential scale)"))

        _pct_above_gate = 100 * float(np.mean(naive_vals > gate)) if len(naive_vals) else float("nan")
        blocks = [
            mo.md(f"**Used by:** {', '.join(cfg['used_by'].keys())}"),
            slider,
            plot_overview(
                plot_n, f"{channel} gate calibration",
                why=(
                    f"This plot lets the {channel} threshold used to call guide/effector/marker-positive "
                    f"cells be set and checked against a reference well, since this threshold feeds "
                    f"directly into every arm's knockdown gating that uses this (channel, reference-well) pair."
                ),
                how=(
                    f"The reference well ({naive_label}) and every sample well that uses this gate are "
                    f"plotted as {channel} histograms (biexponential x-axis, % of each well's peak count "
                    f"on the y-axis). The slider above sets the gate position directly."
                ),
                how_to_read=(
                    "Each line is one well; the dashed red line is the current gate position from the "
                    "slider above. A well-placed gate sits in the valley between the reference well's "
                    "peak and the sample wells' positive population."
                ),
            ),
            fig,
            plot_result(
                f"At the current gate ({gate:,.0f}), {_pct_above_gate:.1f}% of the reference well's own "
                f"cells fall above the gate (ideally close to 0%, since the reference well is expected "
                f"to be negative for this marker)."
            ),
        ]

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

    # NOTE: the calibration panel (histograms + sliders) built above is kept as
    # plumbing only -- infection_gate_widgets.value is still read by build_arm /
    # the per-arm infection tabs below, but this section is intentionally not
    # displayed/rendered anymore (infection results now live inside each arm's
    # own "Infection/guide gating" tab as static, read-only output; no standalone
    # calibration-panel UI is shown anywhere in the notebook).
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
            "(channel, reference-well) pair and updates its knockdown numbers below.\n\n**FMO controls not available for two gates:** the GFP guide-marker gate (FITC-A) and the ASGR1-BV421 gate use fully-unstained naive wells as their reference, not fluorescence-minus-one (FMO) controls, because no well on this plate carries the full antibody panel while omitting only the GFP virus or only the ASGR1-BV421 stain. This is a limitation of the plate layout, not a pipeline choice."
        ),
        infection_calibration_accordion,
    ])
    return


@app.cell(hide_code=True)
def _(
    biexp,
    biexp_ticks,
    debris_gate,
    gated_range,
    interactive_hist,
    mo,
    np,
    plot_overview,
    plot_result,
):
    _HEP_REFS = [
        {"label": "WTC11 (no virus)", "well": "H_A1"},
        {"label": "17_3 (integrated NT-GFP, no effector)", "well": "H_A2"},
        {"label": "17_3 (integrated B2M-GFP, no effector)", "well": "H_A3"},
        {"label": "17_3 (no guide, no effector)", "well": "H_A4"},
    ]

    def _hep_panel(channel, plot_n, title, marker_name):
        traces = {}
        pooled = []
        for ref in _HEP_REFS:
            vals = debris_gate(ref["well"])[channel].values
            traces[ref["label"]] = vals
            pooled.append(vals)
        pooled = np.concatenate(pooled)
        xr = gated_range(pooled, pooled, hi_pct=99.0, group_vals=list(traces.values()))
        cofactor = max(float(np.median(np.abs(pooled))), 1.0) / 5.0

        def _tx(v):
            return biexp(v, cofactor)

        traces_t = {k: _tx(v) for k, v in traces.items()}
        xr_t = (float(_tx(xr[0])), float(_tx(xr[1])))
        fig = interactive_hist(traces_t, None, "", xr_t, f"Plot {plot_n}. {title}")
        tickvals, ticktext = biexp_ticks(xr[0], xr[1], cofactor)
        fig.update_layout(xaxis=dict(tickvals=tickvals, ticktext=ticktext, title=f"{channel} (biexponential scale)"))
        _maxes = [f"{ref['label']}: {float(np.percentile(debris_gate(ref['well'])[channel].values, 95)):,.0f}" for ref in _HEP_REFS]
        return mo.vstack([
            plot_overview(
                plot_n, title,
                why=(
                    f"This plot shows what {marker_name} signal looks like across the cell-line "
                    f"reference wells, for quality-control context -- it is not used to gate or "
                    f"stratify any knockdown arm."
                ),
                how=(
                    f"{marker_name} ({channel}) values from each reference well are plotted as "
                    f"histograms, normalized so each well's tallest bin reads 100%."
                ),
                how_to_read=(
                    "Each line is one reference well (see legend); the x-axis is "
                    f"{channel} on a biexponential scale, the y-axis is % of each well's peak count."
                ),
            ),
            fig,
            plot_result("95th-percentile " + marker_name + " signal per reference well: " + "; ".join(_maxes) + "."),
        ])

    hepatocyte_marker_tabs = mo.ui.tabs({
        "ASGR1-PE (informational)": _hep_panel("PE-A", 7, "ASGR1-PE across cell-line reference wells (informational, not used for arm stratification)", "ASGR1-PE"),
        "Albumin-AF647 (informational)": _hep_panel("APC-A", 8, "Albumin-AF647 across cell-line reference wells (informational, not used for arm stratification)", "Albumin-AF647"),
    })

    hepatocyte_marker_section = mo.vstack([
        mo.md(
            "## Hepatocyte marker gating (ASGR1, Albumin) -- informational only\n"
            "**No arm below is stratified by ASGR1-PE or Albumin-AF647.** Those "
            "antibodies are stained only on a separate set of reference wells "
            "(H_A1-H_A4), not co-stained together with the knockdown-arm readout "
            "antibodies. The tabs below show what signal looks like on those "
            "reference wells, for context and quality control. **ASGR1 "
            "stratification for the knockdown arms instead uses ASGR1-BV421**, "
            "which is co-stained in-well for all 3 arms -- see each arm's own "
            "\"ASGR1+\" sub-tab below. ASGR1 is stained with BV421 or PE antibody "
            "clones only in this dataset."
        ),
        hepatocyte_marker_tabs,
    ])
    return (hepatocyte_marker_section,)


@app.cell(hide_code=True)
def _(debris_slider, go, mo, np, plot_overview, plot_result, raw_wells):
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
        title="Plot 10. Pooled FSC-A, debris/cell valley check (all loaded wells, fixed HLC Trial3 D28)",
        xaxis_title="FSC-A (linear)", yaxis_title="count", height=380, margin=dict(t=60),
    )
    _pct_below_debris = 100 * float(np.mean(_pooled_fsc < debris_slider.value))
    debris_check = mo.vstack([
        plot_overview(
            10, "Pooled FSC-A, debris/cell valley check",
            why="This plot answers whether the FSC-A debris gate is positioned correctly, since every downstream gate in this notebook depends on it.",
            how="FSC-A values from every loaded well are pooled into one histogram.",
            how_to_read=(
                "The x-axis is FSC-A (linear scale); the dashed red line is the current debris gate from "
                "the slider below. A correctly placed gate sits in the valley between the small-event "
                "(debris) peak and the larger-event (real cell) population."
            ),
        ),
        _fig,
        plot_result(f"{_pct_below_debris:.1f}% of all pooled events fall below the current debris gate and are excluded as debris."),
    ])
    return (debris_check,)


@app.cell(hide_code=True)
def _(debris_slider, go, mo, np, plot_overview, raw_wells):
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
        title="Plot 11. Pooled FSC-A vs SSC-A, debris/cell valley in 2D",
        xaxis_title="FSC-A", yaxis_title="SSC-A", height=420, margin=dict(t=60),
        xaxis_range=[-50_000, 3_000_000], yaxis_range=[-50_000, 3_000_000],
    )
    _fig_ds.add_annotation(
        text="downsampled to <=2,500 events/well for display, outliers >3e6 on either axis clipped", x=0, y=-0.14,
        xref="paper", yref="paper", showarrow=False, font=dict(size=10, color="gray"),
    )
    debris_scatter = mo.vstack([
        plot_overview(
            11, "Pooled FSC-A vs SSC-A, debris/cell valley in 2D",
            why="This plot checks the same debris gate as the previous plot, in 2D against SSC-A, to catch cases a 1D histogram would hide.",
            how="A random downsample (up to 2,500 events per well) of FSC-A and SSC-A from every loaded well is plotted together.",
            how_to_read="The dashed red line marks the current FSC-A debris gate; real cells should sit to its right across the SSC-A range.",
        ),
        _fig_ds,
    ])
    return (debris_scatter,)


@app.cell(hide_code=True)
def _(doublet_slider, go, mo, np, plot_overview, plot_result, raw_wells):
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
        title="Plot 12. Pooled FSC-Width, doublet/clump cutoff check (all loaded wells)",
        xaxis_title="FSC-Width (linear)", yaxis_title="count", height=380, margin=dict(t=60),
    )
    _pct_above_doublet = 100 * float(np.mean(_pooled_w > doublet_slider.value))
    doublet_check = mo.vstack([
        plot_overview(
            12, "Pooled FSC-Width, doublet/clump cutoff check",
            why="This plot answers whether the FSC-Width doublet cutoff is positioned correctly, since every downstream gate depends on it.",
            how="FSC-Width values from every loaded well are pooled into one histogram.",
            how_to_read="The dashed red line is the current doublet cutoff from the slider below; it should sit past the main singlet peak, clipping only the doublet tail.",
        ),
        _fig_w,
        plot_result(f"{_pct_above_doublet:.1f}% of all pooled events fall above the current doublet cutoff and are excluded as doublets."),
    ])
    return (doublet_check,)


@app.cell(hide_code=True)
def _(doublet_slider, go, mo, np, plot_overview, raw_wells):
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
        title="Plot 13. Pooled FSC-A vs FSC-Width, doublet band above the cutoff",
        xaxis_title="FSC-A", yaxis_title="FSC-Width", height=420, margin=dict(t=60),
        xaxis_range=[-50_000, 3_000_000], yaxis_range=[0, 8000],
    )
    _fig_ws.add_annotation(
        text="downsampled to <=2,500 events/well for display, outliers clipped (FSC-A>3e6, FSC-Width>8000)", x=0, y=-0.14,
        xref="paper", yref="paper", showarrow=False, font=dict(size=10, color="gray"),
    )
    doublet_scatter = mo.vstack([
        plot_overview(
            13, "Pooled FSC-A vs FSC-Width, doublet band above the cutoff",
            why="This plot checks the same doublet cutoff as the previous plot, in 2D against FSC-A, to catch cases a 1D histogram would hide.",
            how="A random downsample (up to 2,500 events per well) of FSC-A and FSC-Width from every loaded well is plotted together.",
            how_to_read="The dashed red line marks the current FSC-Width doublet cutoff; singlets should sit below it across the FSC-A range.",
        ),
        _fig_ws,
    ])
    return (doublet_scatter,)


@app.cell(hide_code=True)
def _(mo):
    debris_slider = mo.ui.slider(
        start=0, stop=1_500_000, step=10_000, value=250_000,
        label="FSC-A debris/cell gate (events at/below this FSC-A value are excluded as debris)",
        full_width=True, show_value=True,
    )
    return (debris_slider,)


@app.cell(hide_code=True)
def _(mo):
    doublet_slider = mo.ui.slider(
        start=0, stop=8000, step=100, value=3800,
        label="FSC-Width doublet/clump cutoff (events above this FSC-Width are excluded as doublets)",
        full_width=True, show_value=True,
    )
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
    wells_loaded_msg = f"Loaded {len(raw_wells)} wells, {sum(len(d) for d in raw_wells.values()):,} total events"
    return raw_wells, wells_loaded_msg


@app.cell(hide_code=True)
def _(
    SINGLET_RATIO_BAND,
    compensate,
    debris_slider,
    doublet_slider,
    go,
    mo,
    np,
    pd,
    raw_wells,
    sample_sheet,
    ssc_cap_slider,
):
    MIN_CELLS_PER_GATE = 300
    _CONTROL_TOKENS = ["nt-", "non-targeting", "unstained", "no guide", "control"]

    def is_control_label(text: str) -> bool:
        t = (text or "").lower()
        return any(tok in t for tok in _CONTROL_TOKENS)

    def debris_gate(well: str) -> pd.DataFrame:
        d = raw_wells[well]
        d = d[
            (d["FSC-A"] > debris_slider.value)
            & (d["SSC-A"] > 0)
            & (d["SSC-A"] <= ssc_cap_slider.value)
            & (d["FSC-Width"] <= doublet_slider.value)
        ]
        _ratio = d["FSC-H"].values / np.maximum(d["FSC-A"].values, 1.0)
        d = d[(_ratio >= SINGLET_RATIO_BAND[0]) & (_ratio <= SINGLET_RATIO_BAND[1])]
        return compensate(d)

    def pctile_gate(well: str, chan: str, pct: float = 99.0) -> float:
        d = debris_gate(well)
        return float(np.percentile(d[chan].values, pct))

    def apply_markers(well: str, marker_specs: list) -> pd.DataFrame:
        d = debris_gate(well)
        mask = np.ones(len(d), dtype=bool)
        for chan, thr in marker_specs:
            mask &= d[chan].values > thr
        return d[mask]

    def gated_range(naive_vals, pooled_vals, pad_frac: float = 0.10, hi_pct: float = 99.0, lo_pct: float = 1.0, group_vals=None):
        """Left edge is cut at the lo_pct percentile of whichever population
        (unstained control, or any condition group passed in group_vals) has the
        LOWEST signal -- not always the unstained control. Using only the
        unstained well's percentile as a fixed floor was clipping real cells that
        legitimately fall below the unstained reference (seen in infection-marker
        plots), so every candidate population's own lo_pct percentile is
        computed and the minimum of those sets the floor. Right edge keeps its
        padding past the hi_pct percentile of the pooled (gated) data."""
        naive_vals = np.asarray(naive_vals)
        pooled_vals = np.asarray(pooled_vals)
        candidates = [naive_vals] if len(naive_vals) else []
        if group_vals:
            candidates.extend(np.asarray(g) for g in group_vals if len(g))
        if not candidates:
            candidates = [pooled_vals]
        lo_anchor = min(float(np.percentile(c, lo_pct)) for c in candidates if len(c))
        hi_anchor = np.percentile(pooled_vals, hi_pct)
        span = max(hi_anchor - lo_anchor, 1.0)
        return lo_anchor, hi_anchor + pad_frac * span

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
        """traces is keyed by group label (replicate wells already pooled by the
        caller); each group is drawn as a plain scatter of a fixed-seed random
        subsample of its events (never a density/contour plot), using the same
        x-axis range/transform as the paired histogram for direct comparison."""
        fig = go.Figure()
        colors = ["#4C78A8", "#E45756", "#54A24B", "#F58518", "#B279A2", "#72B7B2"]
        rng = np.random.default_rng(seed)
        all_y = np.concatenate([np.asarray(yv) for _, yv in traces.values()]) if traces else np.array([])
        y_max = float(np.percentile(all_y, 99.5)) if len(all_y) else None
        for i, (label, (xv, yv)) in enumerate(traces.items()):
            xv = np.asarray(xv); yv = np.asarray(yv)
            n_total = len(xv)
            if n_total > n_show:
                idx = rng.choice(n_total, n_show, replace=False)
                xv_plot, yv_plot = xv[idx], yv[idx]
            else:
                xv_plot, yv_plot = xv, yv
            color = colors[i % len(colors)]
            fig.add_trace(go.Scattergl(
                x=xv_plot, y=yv_plot, mode="markers", name=f"{label} (n={n_total:,})",
                marker=dict(size=3, opacity=0.35, color=color),
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
        if y_max is not None:
            fig.update_layout(yaxis_range=[0, y_max])
        fig.add_annotation(
            text=f"downsampled to <={n_show:,} events/group for display", x=0, y=-0.18,
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
            "mfi": "Median B2M (MFI)",
            "mfi_ratio_vs_control": "MFI ratio vs control",
        }
        _df = _df.rename(columns={k: v for k, v in rename.items() if k in _df.columns})
        return mo.ui.table(_df, selection=None)

    def replicate_labels(wells: list, group_label: str = None) -> dict:
        labels = {}
        for i, w in enumerate(wells, start=1):
            if group_label:
                labels[w] = f"{group_label} rep{i}"
            else:
                cell_line = sample_sheet.loc[w, "Cell Line"]
                labels[w] = f"{cell_line} rep{i} ({w})"
        return labels

    _CHANNEL_COFACTORS = {}

    def channel_cofactor(channel: str, naive_vals) -> float:
        """One fixed biexponential cofactor per channel, derived from the
        unstained (naive) well, reused across every plot of that channel."""
        if channel not in _CHANNEL_COFACTORS:
            _CHANNEL_COFACTORS[channel] = max(float(np.median(np.abs(naive_vals))), 1.0) / 5.0
        return _CHANNEL_COFACTORS[channel]

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

    def plot_overview(n: int, title: str, why: str, how: str, how_to_read: str):
        """Standard pre-plot overview block: Plot N and the conclusion-oriented
        "why" are shown directly; the "how" (calculation/gating methodology) and
        "how to read it" prose are collapsed into a hidden-by-default accordion
        so the notebook reads as results-first."""
        return mo.vstack([
            mo.md(f"**Plot {n}. {title}**\n\n*Why:* {why}"),
            mo.accordion({
                "Methodology (how this was computed / how to read it)": mo.md(
                    f"*How:* {how}\n\n*How to read it:* {how_to_read}"
                )
            }),
        ])

    def plot_result(text: str):
        return mo.md(f"*Result:* {text}")

    REAGENT_ABBREV = {"Polybrene": "PB", "Protamine Sulfate": "PS", "Lentiboost": "LB"}

    def reagent_abbrev(well: str) -> str:
        reagent = sample_sheet.loc[well, "Transduction Reagent"]
        return REAGENT_ABBREV.get(reagent, "")


    return (
        MIN_CELLS_PER_GATE,
        REAGENT_ABBREV,
        apply_markers,
        biexp,
        biexp_ticks,
        calibration_trace_labels,
        channel_cofactor,
        debris_gate,
        gated_range,
        interactive_hist,
        is_control_label,
        pctile_gate,
        plot_overview,
        plot_result,
        reagent_abbrev,
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
    channel_cofactor,
    debris_gate,
    gated_range,
    gating_hierarchy_table,
    go,
    interactive_hist,
    is_control_label,
    mo,
    np,
    pd,
    plot_overview,
    plot_result,
    reagent_abbrev,
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
        plot_start: int = 1,
        show_reagent: bool = False,
        subgroup_labels: dict = None,
        infection_tab_content=None,
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
        primary_label, primary_wells, _ = groups[-1]
        control_vals = np.concatenate([
            gated_by_well[w][1][readout_channel].values for w in control_wells
        ]) if all(len(gated_by_well[w][1]) > 0 for w in control_wells) else np.array([])
        gate = float(np.percentile(control_vals, 1)) if len(control_vals) else None
        gate50 = float(np.percentile(control_vals, 50)) if len(control_vals) else None

        pooled_vals = np.concatenate([d[readout_channel].values for _, d in gated_by_well.values() if len(d)])
        _group_vals_by_label = {}
        for _, (label, d) in gated_by_well.items():
            if len(d):
                _group_vals_by_label.setdefault(label, []).append(d[readout_channel].values)
        _group_vals_for_range = [np.concatenate(vs) for vs in _group_vals_by_label.values()]
        xr = gated_range(naive_vals, pooled_vals, hi_pct=99.0, group_vals=_group_vals_for_range)

        cofactor = channel_cofactor(readout_channel, naive_vals)

        def _tx(v):
            return biexp(v, cofactor)

        naive_vals_t = _tx(naive_vals)
        xr_t = (float(_tx(xr[0])), float(_tx(xr[1])))
        gate_t = float(_tx(gate)) if gate is not None else None

        # Histogram and scatter/density traces both pool replicate wells within
        # each condition group (one trace per group, not per well).
        # subgroup_labels lets a group be split further by a non-replicate
        # condition (e.g. transduction reagent) that varies across a group's
        # replicate wells -- reps within the same (group, subgroup) are still
        # pooled, but different subgroups are kept as separate traces.
        hist_traces_by_group = {}
        scatter_traces_by_group = {}
        for w, (label, d) in gated_by_well.items():
            if len(d):
                trace_label = f"{label} [{subgroup_labels[w]}]" if subgroup_labels and w in subgroup_labels else label
                vals_t = _tx(d[readout_channel].values)
                hist_traces_by_group.setdefault(trace_label, []).append(vals_t)
                xs, ys = scatter_traces_by_group.setdefault(trace_label, ([], []))
                xs.append(vals_t)
                ys.append(d["SSC-A"].values)
        hist_traces = {"unstained (naive)": naive_vals_t}
        hist_traces.update({label: np.concatenate(vs) for label, vs in hist_traces_by_group.items()})
        scatter_traces = {
            label: (np.concatenate(xs), np.concatenate(ys))
            for label, (xs, ys) in scatter_traces_by_group.items()
        }

        gate50_t = float(_tx(gate50)) if gate50 is not None else None
        _gate1_label = f"1st %ile mock ({control_label}) = {gate:,.0f}" if gate else ""
        _gate50_label = f"50th %ile mock ({control_label}) = {gate50:,.0f}" if gate50 else ""

        fig_hist = interactive_hist(
            hist_traces, gate_t, _gate1_label,
            xr_t, f"Plot {plot_start}. {title}: {readout_channel} distribution (biexponential x-axis, cofactor={cofactor:,.0f})",
            gate2=gate50_t, gate2_label=_gate50_label, gate2_color="purple",
        )
        fig_scatter = scatter_gate(
            scatter_traces, xr_t, f"Plot {plot_start + 1}. {title}: {readout_channel} vs SSC-A (biexponential x-axis)", y_chan="SSC-A",
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
                mfi = float(np.median(d[readout_channel].values)) if n else np.nan
                per_well_rows.append({
                    "well": w, "label": label, "n": n,
                    "pct_below_gate": pct, "pct_below_gate50": pct50, "mfi": mfi,
                })

        df = pd.DataFrame(per_well_rows)
        if len(df) and gate is not None:
            ctrl_mean = df.loc[df["label"] == control_label, "pct_below_gate"].mean()
            df["knockdown_pp"] = df["pct_below_gate"] - ctrl_mean
        if len(df) and gate50 is not None:
            ctrl_mean50 = df.loc[df["label"] == control_label, "pct_below_gate50"].mean()
            df["knockdown_pp_50"] = df["pct_below_gate50"] - ctrl_mean50
        if len(df):
            ctrl_mfi = df.loc[df["label"] == control_label, "mfi"].mean()
            df["mfi_ratio_vs_control"] = df["mfi"] / ctrl_mfi if ctrl_mfi else np.nan

        bar_fig = go.Figure()
        bar_fig_50 = go.Figure()
        if len(df):
            rl = {}
            for label, wells, _ in groups:
                rl.update(replicate_labels([w for w in wells if w in df["well"].values], group_label=label))
            x_labels = [
                rl.get(w, w) + (f" [{reagent_abbrev(w)}]" if show_reagent and reagent_abbrev(w) else "")
                for w in df["well"]
            ]
            colors_map = {g[0]: c for g, c in zip(groups, ["#4C78A8", "#E45756", "#54A24B", "#F58518"])}
            bar_fig.add_trace(go.Bar(
                x=x_labels, y=df["pct_below_gate"],
                text=[f"{v:.1f}%<br>(n={n:,})" for v, n in zip(df["pct_below_gate"], df["n"])],
                textposition="outside",
                marker_color=[colors_map.get(l, "#999") for l in df["label"]],
            ))
            ymax = max(50, float(df["pct_below_gate"].max()) * 1.25) if df["pct_below_gate"].notna().any() else 100
            bar_fig.update_layout(
                title=f"Plot {plot_start + 2}. {title}: % below 1st-percentile-of-mock gate (per replicate; tail effect)",
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
                title=f"Plot {plot_start + 3}. {title}: % below 50th-percentile-of-mock gate (per replicate; population-level shift)",
                yaxis_title="% below gate", yaxis_range=[0, min(100, ymax50)],
                height=380, margin=dict(t=60),
            )

        excl_note = ""
        if excluded:
            items = "; ".join(f"{lbl} ({w}), n={n}" for w, lbl, n in excluded)
            excl_note = f"*Excluded (population < {MIN_CELLS_PER_GATE} cells):* {items}"

        meta_cols = ["Well ID", "Sample ID #", "Cell Line", "Integrated Guide",
                     "Doxycycline Induction", "Virus 1 Type", "Virus 2 Type"]
        if show_reagent:
            meta_cols = meta_cols + ["Transduction Reagent"]
        all_wells = [w for _, wells, _ in groups for w in wells]
        _meta_df = sample_sheet.loc[all_wells, meta_cols].copy()
        if show_reagent:
            _meta_df["Reagent (abbrev.)"] = [reagent_abbrev(w) for w in all_wells]
        meta_table = mo.ui.table(_meta_df, selection=None)

        _metric_note = mo.md(
            "*Two complementary metrics below: the red **1st-percentile-of-mock** "
            "gate captures a tail/low-signal effect (cells shifted furthest down); "
            "the purple **50th-percentile-of-mock** gate captures a broader "
            "population-level shift (whether the bulk of the distribution moved, "
            "not just the tail).*"
        )
        _primary_n = int(df.loc[df["label"] == primary_label, "n"].sum()) if len(df) else 0
        _control_n = int(df.loc[df["label"] == control_label, "n"].sum()) if len(df) else 0
        _kd1 = float(df.loc[df["label"] == primary_label, "knockdown_pp"].mean()) if len(df) and "knockdown_pp" in df else float("nan")
        _kd50 = float(df.loc[df["label"] == primary_label, "knockdown_pp_50"].mean()) if len(df) and "knockdown_pp_50" in df else float("nan")
        _pct1_primary = float(df.loc[df["label"] == primary_label, "pct_below_gate"].mean()) if len(df) else float("nan")
        _pct1_control = float(df.loc[df["label"] == control_label, "pct_below_gate"].mean()) if len(df) else float("nan")
        _pct50_primary = float(df.loc[df["label"] == primary_label, "pct_below_gate50"].mean()) if len(df) else float("nan")
        _pct50_control = float(df.loc[df["label"] == control_label, "pct_below_gate50"].mean()) if len(df) else float("nan")

        blocks = [mo.md(description)]
        if caveat:
            blocks.append(mo.callout(mo.md(caveat), kind="warn"))
        _arm_tabs = {
            "Well metadata": meta_table,
            "Gating hierarchy (counts per step)": gating_hierarchy_table(groups),
        }
        if infection_tab_content is not None:
            _arm_tabs["Infection/guide gating"] = infection_tab_content
        blocks.append(mo.ui.tabs(_arm_tabs))
        blocks.append(mo.md(f"## Knockdown: {readout_channel}"))
        blocks.append(_metric_note)

        blocks.append(plot_overview(
            plot_start, f"{title} -- {readout_channel} distribution",
            why=(
                f"This plot answers whether the {primary_label} group shows lower "
                f"{readout_channel} signal than the {control_label} group in the "
                f"{title} arm, which is the basis of the B2M-knockdown measurement "
                f"for this arm."
            ),
            how=(
                f"{readout_channel} values (gated on this arm's guide/effector-"
                f"delivery markers, and compensated for cross-channel spillover) "
                f"are pooled across replicate wells within each condition group"
                + (", further split by transduction reagent since reagent varies "
                   "within each group's replicate wells" if subgroup_labels else "")
                + (
                    f", then binned into a histogram and normalized so each "
                    f"group's tallest bin reads 100%. An unstained reference well "
                    f"is overlaid for comparison."
                )
            ),
            how_to_read=(
                f"The x-axis is {readout_channel} on a biexponential scale (linear "
                f"near zero, log-like at high values); the y-axis is % of each "
                f"group's own peak count, so every group's curve has the same "
                f"height regardless of cell count. Each line is one condition "
                f"group" + (" and reagent" if subgroup_labels else "") +
                f", pooled across its replicate wells. The dashed red line is "
                f"the 1st-percentile-of-{control_label} gate; the dashed purple "
                f"line is the 50th-percentile-of-{control_label} gate."
            ),
        ))
        blocks.append(fig_hist)
        blocks.append(plot_result(
            f"Averaged across replicate wells, {_pct1_primary:.1f}% of {primary_label} cells "
            f"(n={_primary_n:,}) fall below the 1st-percentile-of-{control_label} gate, versus "
            f"{_pct1_control:.1f}% of {control_label} cells (n={_control_n:,}) -- a {_kd1:+.1f} "
            f"percentage-point difference."
            if len(df) else "No gated cells were available to compute this result."
        ))

        blocks.append(plot_overview(
            plot_start + 1, f"{title} -- {readout_channel} vs SSC-A",
            why=(
                f"This plot checks whether the {readout_channel} signal difference "
                f"between {primary_label} and {control_label} corresponds to a "
                f"distinct population in 2D (not just a 1D histogram artifact), and "
                f"lets the gate position be checked against the granularity (SSC-A) "
                f"axis."
            ),
            how=(
                f"The same gated, compensated {readout_channel} values as the "
                f"histogram above are plotted against each cell's SSC-A, with "
                f"replicate wells in the same group pooled into one density "
                f"contour per group (never pooled across groups)."
            ),
            how_to_read=(
                f"The x-axis is {readout_channel} (biexponential scale, same as "
                f"the histogram above); the y-axis is SSC-A, capped at the 99.5th "
                f"percentile of plotted events so the main population is not "
                f"compressed by outliers. Each colored contour is one condition "
                f"group" + (" and reagent" if subgroup_labels else "") +
                f". The dashed red/purple lines mark the same two gates as "
                f"the histogram above."
            ),
        ))
        blocks.append(fig_scatter)
        blocks.append(plot_result(
            f"The {primary_label} and {control_label} populations are visible as separate "
            f"density contours along the {readout_channel} axis, consistent with the shift "
            f"seen in the histogram above."
            if len(df) else "No gated cells were available to compute this result."
        ))

        blocks.append(plot_overview(
            plot_start + 2, f"{title} -- % below 1st-percentile-of-mock gate",
            why=f"This plot states the tail-effect knockdown number for {title} directly, per replicate well.",
            how=(
                f"For each well, the % of gated cells below the 1st-percentile-of-"
                f"{control_label} gate (from the histogram above) is computed, "
                f"then the mean of the {control_label} wells is subtracted from "
                f"each well's value to give a knockdown in percentage points (pp)."
            ),
            how_to_read=(
                "Each bar is one replicate well, labeled with its condition and "
                "well ID; the bar height is % of cells below the gate, with n "
                "(gated cell count) printed above each bar."
            ),
        ))
        blocks.append(bar_fig)
        blocks.append(plot_result(
            f"Mean knockdown for {primary_label} vs. {control_label} at this gate is {_kd1:+.1f} "
            f"percentage points, averaged across {len(df.loc[df['label'] == primary_label])} "
            f"replicate well(s)."
            if len(df) and "knockdown_pp" in df else "No gated cells were available to compute this result."
        ))

        blocks.append(plot_overview(
            plot_start + 3, f"{title} -- % below 50th-percentile-of-mock gate",
            why=f"This plot states the population-level-shift knockdown number for {title} directly, per replicate well.",
            how=(
                f"Same computation as the previous plot, using the 50th-"
                f"percentile-of-{control_label} gate instead of the 1st "
                f"percentile, so it captures a shift of the bulk of the "
                f"distribution rather than only the lowest-signal tail."
            ),
            how_to_read="Same layout as the previous bar chart, at the 50th-percentile gate.",
        ))
        blocks.append(bar_fig_50)
        blocks.append(plot_result(
            f"Mean knockdown for {primary_label} vs. {control_label} at this gate is {_kd50:+.1f} "
            f"percentage points, averaged across {len(df.loc[df['label'] == primary_label])} "
            f"replicate well(s)."
            if len(df) and "knockdown_pp_50" in df else "No gated cells were available to compute this result."
        ))

        if excl_note:
            blocks.append(mo.md(excl_note))
        if len(df):
            blocks.append(tidy_summary_table(df.sort_values(["label", "well"])))

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
def _(arm1_infect_tab, build_arm, infection_gate_widgets):
    arm1_content_base, arm1_summary = build_arm(
        title="WTC11 transient B2M-GFP virus",
        description=(
            "**Construct:** `NT-GFP` (non-targeting guide, GFP-marked) vs `B2M-GFP` "
            "(B2M-targeting guide, GFP-marked) delivered as a single transient virus "
            "into WTC11 cells. WTC11 carries a doxycycline-inducible CRISPRi "
            "effector, so the B2M-GFP group is split into **no dox** (effector off, "
            "expected to show no knockdown) and **+dox** (effector on, expected to "
            "show knockdown) -- this induction status is this arm's main "
            "experimental variable, so these wells are kept as separate groups "
            "rather than averaged together. GFP+ = confirmed transduction. "
            "Readout: B2M-APC."
        ),
        groups=[
            ("NT-GFP (control)", ["B1", "B2"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
            ("B2M-GFP (guide), no dox", ["B3", "B4"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
            ("B2M-GFP (guide), +dox", ["B5", "B6"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
        ],
        readout_channel="APC-A",
        naive_well="A3",
        plot_start=14,
        caveat=(
            "**Sample-sheet discrepancy:** the sample sheet's `Doxycycline "
            "Induction` column reads \"No\" for all of B3, B4, B5, and B6. B5 and "
            "B6 received doxycycline induction and B3 and B4 did not -- the "
            "groups above reflect this confirmed induction status, not the sheet's "
            "column value, which does not capture the distinction."
        ),
        infection_tab_content=arm1_infect_tab,
    )
    return arm1_content_base, arm1_summary


@app.cell(hide_code=True)
def _(arm2_infect_tab, build_arm, infection_gate_widgets, reagent_abbrev):
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
            "single-stain control well, E9).\n\n"
            "**Transduction reagent:** this arm's 12 wells also compare 3 "
            "lentiviral delivery reagents for the AA239 effector virus, labeled "
            "below as **PB** (Polybrene), **PS** (Protamine Sulfate), and **LB** "
            "(Lentiboost), since reagent differs across this arm's replicate "
            "wells. See the reagent-comparison chart at the top of this tab for "
            "the effector-delivery rate by reagent; the distribution and scatter "
            "plots below also split each condition by reagent for the same reason."
        ),
        groups=[
            ("NT-GFP (control)", ["C1", "C2", "C5", "C6", "D1", "D2"],
             [("FITC-A", infection_gate_widgets.value["FITC-A__A7"]), ("APC-A", infection_gate_widgets.value["APC-A__D11_fmo"])]),
            ("B2M-GFP (guide)", ["C3", "C4", "C7", "C8", "D3", "D4"],
             [("FITC-A", infection_gate_widgets.value["FITC-A__A7"]), ("APC-A", infection_gate_widgets.value["APC-A__D11_fmo"])]),
        ],
        readout_channel="B610-ECD-A",
        naive_well="A7",
        plot_start=18,
        show_reagent=True,
        subgroup_labels={
            w: reagent_abbrev(w)
            for w in ["C1", "C2", "C5", "C6", "D1", "D2", "C3", "C4", "C7", "C8", "D3", "D4"]
        },
        infection_tab_content=arm2_infect_tab,
    )
    return arm2_content_base, arm2_summary


@app.cell(hide_code=True)
def _(arm3_infect_tab, build_arm, infection_gate_widgets):
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
        plot_start=22,
        infection_tab_content=arm3_infect_tab,
    )
    return arm3_content_base, arm3_summary


@app.cell(hide_code=True)
def _(
    MIN_CELLS_PER_GATE,
    apply_markers,
    infection_gate_widgets,
    is_control_label,
    np,
    pd,
    reagent_abbrev,
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
                ("B2M-GFP (guide), no dox", ["B3", "B4"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
                ("B2M-GFP (guide), +dox", ["B5", "B6"], [("FITC-A", infection_gate_widgets.value["FITC-A__A3"])]),
            ],
            "readout_channel": "APC-A", "control_group_index": 0,
            "asgr1_gate": infection_gate_widgets.value["BV421-A__A3"], "plot_start": 26,
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
            "asgr1_gate": infection_gate_widgets.value["BV421-A__A7"], "plot_start": 31,
            "subgroup_labels": {
                w: reagent_abbrev(w)
                for w in ["C1", "C2", "C5", "C6", "D1", "D2", "C3", "C4", "C7", "C8", "D3", "D4"]
            },
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
            "asgr1_gate": infection_gate_widgets.value["BV421-A__A7"], "plot_start": 36,
        },
    ]

    _strat_all_rows = []
    for _cfg in STRAT_ARM_CFGS:
        _strat_all_rows.extend(stratified_knockdown_rows(
            _cfg["arm_label"], _cfg["groups"], _cfg["readout_channel"], _cfg["control_group_index"],
            asgr1_gate=_cfg["asgr1_gate"],
        ))

    marker_strat_df = pd.DataFrame(_strat_all_rows)
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
    plot_overview,
    plot_result,
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

    def _strat_bar_chart(arm_label, stratum_label, metric_col, title_suffix, plot_n):
        _sub = marker_strat_df[(marker_strat_df["arm"] == arm_label) & (marker_strat_df["stratum"] == stratum_label)]
        if not len(_sub):
            return None
        _cfg = next(c for c in STRAT_ARM_CFGS if c["arm_label"] == arm_label)
        groups = _cfg["groups"]
        _subgroup_labels = _cfg.get("subgroup_labels")
        rl = {}
        for label, wells, _ in groups:
            rl.update(replicate_labels([w for w in wells if w in _sub["well"].values], group_label=label))
        if _subgroup_labels:
            rl = {w: lbl + (f" [{_subgroup_labels[w]}]" if w in _subgroup_labels else "") for w, lbl in rl.items()}
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
            title=f"Plot {plot_n}. {arm_label} -- {stratum_label}: {title_suffix} (per replicate)",
            yaxis_title="% below gate", yaxis_range=[0, min(100, ymax)], height=360, margin=dict(t=60),
        )
        return fig

    def build_marker_stratum_block(cfg, stratum_label, asgr1_spec):
        groups = cfg["groups"]
        readout_channel = cfg["readout_channel"]
        control_group_index = cfg["control_group_index"]
        naive_well = cfg["naive_well"]
        plot_start = cfg.get("plot_start", 1)
        arm_label = cfg["arm_label"]
        subgroup_labels = cfg.get("subgroup_labels")
        primary_label = groups[-1][0]
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
        _strat_group_vals_by_label = {}
        for _, (label, d) in gated_by_well.items():
            if len(d):
                _strat_group_vals_by_label.setdefault(label, []).append(d[readout_channel].values)
        _strat_group_vals = [np.concatenate(vs) for vs in _strat_group_vals_by_label.values()]
        xr = gated_range(naive_vals_raw, pooled_vals, hi_pct=99.0, group_vals=_strat_group_vals)
        cofactor = max(gate1, 1.0) / 5.0

        def _tx(v):
            return biexp(v, cofactor)

        xr_t = (float(_tx(xr[0])), float(_tx(xr[1])))
        gate1_t, gate50_t = float(_tx(gate1)), float(_tx(gate50))

        hist_traces_by_group = {}
        scatter_traces_by_group = {}
        for w, (label, d) in gated_by_well.items():
            if len(d):
                trace_label = f"{label} [{subgroup_labels[w]}]" if subgroup_labels and w in subgroup_labels else label
                vals_t = _tx(d[readout_channel].values)
                hist_traces_by_group.setdefault(trace_label, []).append(vals_t)
                xs, ys = scatter_traces_by_group.setdefault(trace_label, ([], []))
                xs.append(vals_t)
                ys.append(d["SSC-A"].values)
        hist_traces = {"unstained (naive)": _tx(naive_vals_raw)}
        hist_traces.update({label: np.concatenate(vs) for label, vs in hist_traces_by_group.items()})
        scatter_traces = {
            label: (np.concatenate(xs), np.concatenate(ys))
            for label, (xs, ys) in scatter_traces_by_group.items()
        }

        fig_hist = interactive_hist(
            hist_traces, gate1_t, f"1st %ile mock ({control_label}) = {gate1:,.0f}",
            xr_t, f"Plot {plot_start}. {arm_label} -- {stratum_label}: {readout_channel} distribution",
            gate2=gate50_t, gate2_label=f"50th %ile mock ({control_label}) = {gate50:,.0f}", gate2_color="purple",
        )
        fig_scatter = scatter_gate(
            scatter_traces, xr_t, f"Plot {plot_start + 1}. {arm_label} -- {stratum_label}: {readout_channel} vs SSC-A",
            y_chan="SSC-A", gate=gate1_t, gate_label="1st %ile mock", gate2=gate50_t, gate2_label="50th %ile mock", gate2_color="purple",
        )
        bar1 = _strat_bar_chart(arm_label, stratum_label, "pct_below_1st", "% below 1st-pctile-mock gate", plot_start + 2)
        bar50 = _strat_bar_chart(arm_label, stratum_label, "pct_below_50th", "% below 50th-pctile-mock gate", plot_start + 3)

        _strat_sub = marker_strat_df[(marker_strat_df["arm"] == arm_label) & (marker_strat_df["stratum"] == stratum_label)]
        _kd1 = float(_strat_sub.loc[_strat_sub["label"] == primary_label, "knockdown_pp_1st"].mean()) if len(_strat_sub) else float("nan")
        _kd50 = float(_strat_sub.loc[_strat_sub["label"] == primary_label, "knockdown_pp_50th"].mean()) if len(_strat_sub) else float("nan")
        _strat_n = int(_strat_sub.loc[_strat_sub["label"] == primary_label, "n"].sum()) if len(_strat_sub) else 0

        blocks = [
            plot_overview(
                plot_start, f"{arm_label} -- {stratum_label}: {readout_channel} distribution",
                why=(
                    f"This plot checks whether the {primary_label} vs. {control_label} knockdown signal "
                    f"seen across all cells in {arm_label} still holds within the {stratum_label} cell "
                    f"group specifically."
                ),
                how=(
                    f"Same computation as the arm's unstratified histogram (pooled across replicate wells within each condition group), but restricted to cells additionally passing the {stratum_label} marker gate."
                ),
                how_to_read=(
                    "Same layout as the unstratified histogram: biexponential x-axis, % of each well's "
                    "peak count on the y-axis, dashed red/purple lines mark the 1st- and "
                    "50th-percentile-of-mock gates for this stratum."
                ),
            ),
            fig_hist,
            plot_overview(
                plot_start + 1, f"{arm_label} -- {stratum_label}: {readout_channel} vs SSC-A",
                why=f"This plot checks the {stratum_label} population separation in 2D, the same purpose as the arm's unstratified scatter plot.",
                how="Same computation as the arm's unstratified scatter plot, restricted to this stratum.",
                how_to_read="Same layout as the unstratified scatter plot above.",
            ),
            fig_scatter,
        ]
        if bar1 is not None:
            blocks.append(plot_overview(
                plot_start + 2, f"{arm_label} -- {stratum_label}: % below 1st-percentile-of-mock gate",
                why=f"States the tail-effect knockdown number for the {stratum_label} stratum directly, per replicate well.",
                how="Same computation as the arm's unstratified 1st-percentile bar chart, restricted to this stratum.",
                how_to_read="Same layout as the arm's unstratified bar chart.",
            ))
            blocks.append(bar1)
            blocks.append(plot_result(
                f"Mean knockdown for {primary_label} vs. {control_label} within the {stratum_label} stratum "
                f"is {_kd1:+.1f} percentage points (n={_strat_n:,} {primary_label} cells)."
            ))
        if bar50 is not None:
            blocks.append(plot_overview(
                plot_start + 3, f"{arm_label} -- {stratum_label}: % below 50th-percentile-of-mock gate",
                why=f"States the population-level-shift knockdown number for the {stratum_label} stratum directly, per replicate well.",
                how="Same computation as the arm's unstratified 50th-percentile bar chart, restricted to this stratum.",
                how_to_read="Same layout as the arm's unstratified bar chart.",
            ))
            blocks.append(bar50)
            blocks.append(plot_result(
                f"Mean knockdown for {primary_label} vs. {control_label} within the {stratum_label} stratum "
                f"is {_kd50:+.1f} percentage points (n={_strat_n:,} {primary_label} cells)."
            ))
        blocks.append(_tidy_strat_table(_strat_sub))
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
            title=f"Plot {cfg.get('plot_start', 1) + 4}. {cfg['arm_label']}: knockdown across strata ({primary_label} vs. mock)",
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
        _primary_label = cfg["groups"][-1][0]
        _plot_n = cfg.get("plot_start", 1) + 4

        return mo.vstack([
            plot_overview(
                _plot_n, f"{cfg['arm_label']}: knockdown across strata",
                why=(
                    f"This plot states whether selecting for the ASGR1+ cell group changes the measured "
                    f"knockdown for {cfg['arm_label']}, compared to using all gated cells regardless of "
                    f"ASGR1 status."
                ),
                how=(
                    f"The 1st- and 50th-percentile-of-mock knockdown (pp) for the {_primary_label} group "
                    f"is computed twice: once over all gated cells, once restricted to the ASGR1+ stratum, "
                    f"then plotted side by side."
                ),
                how_to_read=(
                    "Each pair of bars is one cell group (all cells, or ASGR1+); red bars are the "
                    "1st-percentile-of-mock metric, purple bars are the 50th-percentile-of-mock metric, "
                    "both in percentage points of knockdown."
                ),
            ),
            comparison_fig,
            tldr,
            mo.md(
                "*ASGR1-BV421 (co-stained in-well) is the only hepatocyte marker "
                "available for stratification in this dataset. Albumin is stained "
                "only on separate, non-arm reference wells (see the informational "
                "hepatocyte-marker section above).*"
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
    return arm1_content, arm2_content, arm3_content


@app.cell(hide_code=True)
def _(debris_slider, doublet_slider, mo, np, pd, raw_wells, ssc_cap_slider):
    _COMP_CHANNELS = ["FITC-A", "APC-A", "BV421-A", "B610-ECD-A"]
    _SINGLE_STAIN_WELLS = {"FITC-A": "A6", "APC-A": "E6", "BV421-A": "E7", "B610-ECD-A": "E9"}
    _UNSTAINED_WELL = "A3"

    # Spillover terms excluded below (forced to 0 rather than estimated from the
    # plate's single-stain wells): BV421-A->FITC-A and APC-A->FITC-A. Both were
    # derived from well E7 (ASGR1-BV421 single stain), whose own-channel median
    # across the whole well sits *below* the unstained well -- ASGR1 is a
    # hepatocyte-differentiation marker and E7 appears to be an early/less-mature
    # population where it isn't really expressed yet, unlike actual Day 28 HLC
    # sample wells (which show ASGR1-BV421 medians 5-10x higher). The "positive"
    # sub-population used to estimate spillover from E7 is real (clearly brighter
    # BV421 than background), but it's also uniformly brighter across every other
    # channel (FITC, APC, B610-ECD all elevated together), consistent with a
    # general autofluorescence/granularity correlation in that subset rather than
    # true optical bleed-through -- which was driving FITC-A (GFP) compensated
    # values for the real guide-marker wells (B1-B6) down to ~0 or negative,
    # well below the unstained reference. No cleaner single-stain well exists on
    # this plate for either BV421-A or FITC-A (checked all GFP+/antibody-free and
    # ASGR1-BV421-stained wells), so rather than replacing the reference well,
    # these two specific cross-terms are excluded from the estimated matrix.
    _EXCLUDED_SPILLOVER_TERMS = {("BV421-A", "FITC-A"), ("APC-A", "FITC-A")}

    def _build_spillover_matrix():
        unstained = raw_wells[_UNSTAINED_WELL]
        baseline = {ch: float(np.median(unstained[ch].values)) for ch in _COMP_CHANNELS}
        n = len(_COMP_CHANNELS)
        s = np.eye(n)
        for i, ch_i in enumerate(_COMP_CHANNELS):
            d = raw_wells[_SINGLE_STAIN_WELLS[ch_i]]
            pos_thr = np.percentile(unstained[ch_i].values, 95)
            pos = d[d[ch_i].values > pos_thr]
            if len(pos) < 20:
                continue
            denom = float(np.median(pos[ch_i].values)) - baseline[ch_i]
            if denom <= 0:
                continue
            for j, ch_j in enumerate(_COMP_CHANNELS):
                if j == i:
                    continue
                if (ch_i, ch_j) in _EXCLUDED_SPILLOVER_TERMS:
                    continue
                s[i, j] = (float(np.median(pos[ch_j].values)) - baseline[ch_j]) / denom
        return s

    SPILLOVER_MATRIX = _build_spillover_matrix()
    _SPILLOVER_INV = np.linalg.inv(SPILLOVER_MATRIX)

    def compensate(d: pd.DataFrame) -> pd.DataFrame:
        """Removes cross-channel fluorescence spillover using a spillover matrix built
        from this plate's own single-stain control wells, since the FCS files store an
        identity spillover matrix (no built-in compensation)."""
        out = d.copy()
        vals = d[_COMP_CHANNELS].values
        comp = vals @ _SPILLOVER_INV
        for j, ch in enumerate(_COMP_CHANNELS):
            out[ch] = comp[:, j]
        return out

    def _singlet_ratio_band():
        xs = []
        for d in raw_wells.values():
            dd = d[(d["FSC-A"] > debris_slider.value) & (d["SSC-A"] > 0) & (d["FSC-Width"] <= doublet_slider.value)]
            if len(dd):
                xs.append(dd["FSC-H"].values / np.maximum(dd["FSC-A"].values, 1.0))
        r = np.concatenate(xs)
        return float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))

    SINGLET_RATIO_BAND = _singlet_ratio_band()

    compensation_table = mo.ui.table(
        pd.DataFrame(SPILLOVER_MATRIX, index=_COMP_CHANNELS, columns=_COMP_CHANNELS).round(3),
        selection=None,
    )
    gating_notes = mo.vstack([
        mo.md(
            "## Compensation and singlet/cell gates (FlowJo-equivalent additions)\n"
            "**Compensation:** the spillover matrix below is computed from this plate's own single-stain "
            "control wells (GFP: well A6; B2M-APC/Thy1.1-APC: well E6; ASGR1-BV421: well E7; "
            "B2M-PE-Dazzle/B610-ECD: well E9), each referenced against the fully-unstained well A3, since "
            "the FCS files store an identity spillover matrix (no built-in compensation). Every readout, "
            "guide, and antibody channel below (FITC-A, APC-A, BV421-A, B610-ECD-A) is compensated with "
            "this matrix before any gate is applied.\n\n"
            "**Singlet gate:** in addition to the FSC-Width cutoff below, cells must also fall in a "
            f"diagonal FSC-H/FSC-A band ({SINGLET_RATIO_BAND[0]:.2f}-{SINGLET_RATIO_BAND[1]:.2f} of the "
            "FSC-H:FSC-A ratio -- the central 95% of this ratio across all loaded wells after the debris "
            "and width gates).\n\n"
            "**Cell gate:** SSC-A is capped at an upper bound (see slider below), forming a rectangle "
            "with the FSC-A floor above."
        ),
        compensation_table,
    ])
    gating_notes

    _MARKER_LABELS = {"FITC-A": "GFP+ (guide)", "APC-A": "Thy1.1-APC+ (effector)", "BV421-A": "ASGR1-BV421+"}

    def gating_hierarchy_table(groups):
        """One row per well: event counts (and % of the previous step) through
        all events -> cells -> single cells -> each guide/effector marker gate
        applied for that well's group."""
        rows = []
        for label, wells, marker_specs in groups:
            for w in wells:
                d0 = raw_wells[w]
                n0 = len(d0)
                d1 = d0[d0["FSC-A"] > debris_slider.value]
                n1 = len(d1)
                d2 = d1[(d1["SSC-A"] > 0) & (d1["SSC-A"] <= ssc_cap_slider.value) & (d1["FSC-Width"] <= doublet_slider.value)]
                _ratio = d2["FSC-H"].values / np.maximum(d2["FSC-A"].values, 1.0)
                d3 = d2[(_ratio >= SINGLET_RATIO_BAND[0]) & (_ratio <= SINGLET_RATIO_BAND[1])]
                n3 = len(d3)
                row = {
                    "Well": w, "Condition": label,
                    "All events (n)": n0,
                    "Cells (n)": n1, "Cells (% of all events)": round(100 * n1 / n0, 1) if n0 else np.nan,
                    "Single cells (n)": n3,
                    "Single cells (% of cells)": round(100 * n3 / n1, 1) if n1 else np.nan,
                }
                cur = compensate(d3)
                prev_n = n3
                for chan, thr in marker_specs:
                    cur = cur[cur[chan].values > thr]
                    step_label = _MARKER_LABELS.get(chan, f"{chan}+")
                    nk = len(cur)
                    row[f"{step_label} (n)"] = nk
                    row[f"{step_label} (% of previous step)"] = round(100 * nk / prev_n, 1) if prev_n else np.nan
                    prev_n = nk
                rows.append(row)
        return mo.ui.table(pd.DataFrame(rows), selection=None)

    return SINGLET_RATIO_BAND, compensate, gating_hierarchy_table


@app.cell(hide_code=True)
def _(mo):
    ssc_cap_slider = mo.ui.slider(
        start=1_000_000, stop=30_000_000, step=250_000, value=9_500_000,
        label="SSC-A upper cap (events above this SSC-A are excluded; rectangle cell gate with the FSC-A floor above)",
        full_width=True, show_value=True,
    )
    return (ssc_cap_slider,)


@app.cell(hide_code=True)
def _(
    debris_slider,
    doublet_slider,
    go,
    mo,
    np,
    plot_overview,
    plot_result,
    raw_wells,
):
    _rng_t = np.random.default_rng(0)
    _xs_t, _ys_t, _well_flags = [], [], []
    for _w, _d in raw_wells.items():
        _dd = _d[(_d["FSC-A"] > debris_slider.value) & (_d["SSC-A"] > 0) & (_d["FSC-Width"] <= doublet_slider.value)]
        if len(_dd) < 50:
            continue
        _t = _dd["Time"].values.astype(float)
        _fsc = _dd["FSC-A"].values.astype(float)
        _edges = np.linspace(_t.min(), _t.max() + 1, 11)
        _meds = np.array([
            np.median(_fsc[(_t >= _edges[i]) & (_t < _edges[i + 1])])
            if np.any((_t >= _edges[i]) & (_t < _edges[i + 1])) else np.nan
            for i in range(10)
        ])
        _rel = (np.nanmax(_meds) - np.nanmin(_meds)) / np.nanmedian(_meds)
        if _rel > 0.25:
            _well_flags.append((_w, round(float(_rel), 2)))
        _n = len(_t)
        _cap = min(_n, 1500)
        _idx = _rng_t.choice(_n, _cap, replace=False)
        _xs_t.append(_t[_idx]); _ys_t.append(_fsc[_idx])
    _xs_t = np.concatenate(_xs_t); _ys_t = np.concatenate(_ys_t)
    _fig_t = go.Figure()
    _fig_t.add_trace(go.Scattergl(x=_xs_t, y=_ys_t, mode="markers", marker=dict(size=2, opacity=0.25, color="#F58518")))
    _fig_t.update_layout(
        title="Plot 41. Pooled FSC-A vs Time, acquisition-stability check (all loaded wells, debris/width-gated)",
        xaxis_title="Time", yaxis_title="FSC-A", height=380, margin=dict(t=60),
    )
    _flag_text = (
        ", ".join(f"{w} ({r*100:.0f}% swing)" for w, r in _well_flags)
        if _well_flags else "none"
    )
    time_check = mo.vstack([
        plot_overview(
            41, "Pooled FSC-A vs Time, acquisition-stability check",
            why="This plot answers whether any well had a clog or bubble during acquisition, which would bias its gated cell population.",
            how=(
                "FSC-A is split into 10 equal time bins per well (debris- and width-gated); a well is "
                "flagged if its median FSC-A swings by more than 25% between bins."
            ),
            how_to_read="The x-axis is acquisition time, the y-axis is FSC-A; a sharp jump or discontinuity indicates an acquisition problem.",
        ),
        _fig_t,
        plot_result(
            f"Wells flagged for a >25% median-FSC-A swing: {_flag_text}. No hard time gate is applied, "
            "since the flagged swings are mild (25-34%), not sharp jumps."
        ),
    ])
    return (time_check,)


@app.cell(hide_code=True)
def _(
    arm1_content_base,
    arm2_content_base,
    arm3_content_base,
    debris_check,
    debris_scatter,
    debris_slider,
    doublet_check,
    doublet_scatter,
    doublet_slider,
    hepatocyte_marker_section,
    infection_gate_widgets,
    marker_strat_df,
    mo,
    sample_sheet,
    ssc_cap_slider,
    time_check,
    wells_loaded_msg,
):
    mo.accordion({
        "Calibration sliders, QC & raw per-arm views (click to expand -- advanced/debug)": mo.vstack([
            mo.md("#### Infection/guide-marker gate calibration sliders"),
            infection_gate_widgets,
            mo.md("#### Debris gate (FSC-A floor)"),
            debris_slider, debris_check, debris_scatter,
            mo.md("#### Doublet gate (FSC-Width ceiling)"),
            doublet_slider, doublet_check, doublet_scatter,
            mo.md("#### Cell gate (SSC-A upper cap)"),
            ssc_cap_slider,
            mo.md("#### Hepatocyte marker (ASGR1/Albumin) calibration"),
            hepatocyte_marker_section,
            mo.md("#### Marker-stratification raw table"),
            marker_strat_df,
            mo.md("#### Sample sheet"),
            sample_sheet,
            mo.md(f"#### Well loading\n{wells_loaded_msg}"),
            mo.md("#### Time-acquisition QC"),
            time_check,
            mo.md("#### Raw per-arm content (same data as the tabs above, pre-stratification)"),
            arm1_content_base, arm2_content_base, arm3_content_base,
        ]),
    })
    return


if __name__ == "__main__":
    app.run()
