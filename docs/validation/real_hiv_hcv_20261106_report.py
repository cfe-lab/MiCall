"""Generate pseudonymous scientific reports and charts for real HIV/HCV A/B."""
from datetime import datetime, timezone
from html import escape
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs/validation"
WORK = ROOT / ".validation-work"
METRICS = json.loads((WORK / "analysis-private.json").read_text())
DETAILS = {x["alias"]: x for x in
           json.loads((WORK / "per_contig_private.json").read_text())}
DONE = [m for m in METRICS if m["complete"]]
FAILED = [m for m in METRICS if not m["complete"]]
FIG = OUT / "figures"
FIG.mkdir(exist_ok=True)
N = lambda value: f"{value:,}"


def totals():
    return {
        "complete": len(DONE), "target": len(METRICS),
        "hcv": sum(m["group"] == "HCV" for m in DONE),
        "hiv": sum(m["group"] == "HIV" for m in DONE),
        "gained": sum(m["gained_callable"] for m in DONE),
        "lost": sum(m["lost_callable"] for m in DONE),
        "discordant": sum(m["discordant_bases"] for m in DONE),
        "loss_sites": sum(m["depth_positions_decreased"] for m in DONE),
        "threshold_gained": sum(m["positions_cross_100_up"] for m in DONE),
        "threshold_lost": sum(m["positions_cross_100_down"] for m in DONE),
        "low_fraction": sum(DETAILS.get(m["alias"], {}).get(
            "gained_max100_avg_below100", 0) for m in DONE),
        "high_fraction": sum(DETAILS.get(m["alias"], {}).get(
            "gained_avg_at_least100", 0) for m in DONE),
        "sidecar_bugs": sum(m["sidecar_violations"] for m in DONE),
        "pair_bugs": sum(not m["reported_unique_matches"] for m in DONE),
        "placement_gain": sum(m["new_placements"] -
                              m["baseline_placements"] for m in DONE),
    }


def plots():
    if not DONE: return
    labels = [m["alias"] for m in DONE]
    yy = list(range(len(DONE)))
    size = (10.5, max(5.0, .45 * len(DONE) + 1.8))
    fig, ax = plt.subplots(figsize=size)
    ax.barh(yy, [m["gained_callable"] for m in DONE],
            color="#376da6", label="Callable gained")
    ax.barh(yy, [-m["lost_callable"] for m in DONE],
            color="#c47e43", label="Callable lost")
    ax.set_yticks(yy, labels); ax.invert_yaxis()
    ax.axvline(0, color="#777777", lw=.8)
    ax.set_xlabel("Q15 MAX callable positions; fixed input contigs")
    ax.set_title("Coverage changes relative to upstream MiCall")
    ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(FIG / "callable_delta.svg"); plt.close(fig)

    fig, ax = plt.subplots(figsize=size)
    cross = [m["cross_contig_pairs"] for m in DONE]
    intra = [m["within_contig_repeat_pairs"] for m in DONE]
    ax.barh(yy, cross, color="#376da6", label="Across contigs")
    ax.barh(yy, intra, left=cross, color="#ddac69",
            label="Different locations, same contig")
    ax.set_yticks(yy, labels); ax.invert_yaxis()
    ax.set_xlabel("Original read pairs with multiple placements")
    ax.set_title("Where reported ambiguity comes from")
    ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(FIG / "ambiguity.svg"); plt.close(fig)

    fig, ax = plt.subplots(figsize=size)
    for i, m in enumerate(DONE):
        b, c = m["upstream_total_seconds"], m["nonstealing_total_seconds"]
        ax.plot([b, c], [i, i], color="#bac7d3", lw=1.5)
    ax.scatter([m["upstream_total_seconds"] for m in DONE], yy,
               s=45, color="#6b7783", label="Upstream")
    ax.scatter([m["nonstealing_total_seconds"] for m in DONE], yy,
               s=45, color="#376da6", label="Non-stealing")
    ax.set_yticks(yy, labels); ax.invert_yaxis()
    ax.set_xlabel("Wall seconds: remap + pairing + consensus")
    ax.set_title("Exploratory timings; shared concurrent host")
    ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(FIG / "runtime.svg"); plt.close(fig)

    # Matplotlib indents SVG path commands with trailing spaces.
    # Normalize text-only formatting so generated charts pass git diff --check.
    for chart in FIG.glob("*.svg"):
        content = chart.read_text()
        chart.write_text("\n".join(line.rstrip()
                           for line in content.splitlines()) + "\n")


def table_rows():
    rows = []
    for m in DONE:
        net = m["callable_nonstealing"] - m["callable_upstream"]
        ratio = m["nonstealing_total_seconds"] / max(
            1e-9, m["upstream_total_seconds"])
        rows.append((m["alias"], "IVA" if m["panel"] == "archived_IVA"
                     else "Ivanka", str(m["reference_count"]),
                     N(m["original_read_pairs"]),
                     N(m["cross_contig_pairs"]), f"{net:+,}",
                     N(m["discordant_bases"]), f"{ratio:.2f}"))
    return rows


def markdown_table():
    header = ("Case", "Contigs from", "Mapped references", "Read pairs",
              "Cross-contig alternatives", "Callable delta", "Base flips",
              "Wall-time ratio")
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join(["---"] + ["---:"] * (len(header) - 1)) + "|"]
    out += ["| " + " | ".join(row) + " |" for row in table_rows()]
    return "\n".join(out)


def run():
    plots()
    t = totals()
    generation = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    qa = "PASS" if not (t["sidecar_bugs"] or t["pair_bugs"]) else "CHECK FAIL"
    exceptions = []
    for m in FAILED:
        if m["alias"] == "HCV-01":
            exceptions.append(
                "- HCV-01: candidate completed, but upstream MiCall failed "
                "inside sam2aln on one malformed CSV row whose reference "
                "value was missing. This is not a valid paired A/B success.")
        else:
            exceptions.append("- " + m["alias"] + ": " +
                              m.get("reason", "incomplete A/B"))
    if not exceptions: exceptions.append("- None.")
    failures = "\n".join(exceptions)
    notable = sorted(DONE, key=lambda m: m["gained_callable"], reverse=True)[:4]
    highlights = "\n".join(
        "- " + m["alias"] + ": +" + N(m["gained_callable"]) +
        " callable, -" + N(m["lost_callable"]) + " lost, " +
        N(m["discordant_bases"]) + " discordant; " +
        N(DETAILS.get(m["alias"], {}).get("gained_max100_avg_below100", 0)) +
        " newly callable sites with fractional depth <100."
        for m in notable)
    executive = f"""# Real HIV/HCV non-stealing remapper — executive validation

**Decision: promising experimental behavior, but NOT a clinical or production approval.**

**Dataset:** {t['complete']} of {t['target']} completed controlled A/B conditions,
including {t['hcv']} HCV-focused and {t['hiv']} HIV-focused conditions.
The 16 selected conditions represent 14 unique archived specimens; two HCV
specimens were evaluated with both IVA and Ivanka contig panels.

## Findings at a glance

- **{N(t['gained'])} newly callable Q15 MAX positions; {N(t['lost'])} lost; {N(t['discordant'])} callable-base disagreements.** These totals are case-condition coordinates, not statistically independent biological sites.
- Of the newly callable positions, **{N(t['low_fraction'])}** have candidate fractional depth below 100 despite maximum depth meeting that threshold; **{N(t['high_fraction'])}** also have fractional depth at least 100.
- **{N(t['loss_sites'])}** positions show lower maximum depth than upstream, even if some do not affect consensus calls.
- **{N(t['placement_gain'])}** additional mapping placements were reported; these are *alternative hypotheses*, not newly observed molecules.
- **Internal checks: {qa}.** {t['sidecar_bugs']} depth-report disagreements and {t['pair_bugs']} unique-pair accounting mismatches in completed conditions.

![Callable position gains and losses](figures/callable_delta.svg)

## Selected outcomes

{highlights}

**Scientific caution — HCV-05:** For the preserved 16-contig Ivanka S48
panel, 2,594 Q15 MAX positions became callable, three were lost, and one
allele changed from T to C. Upstream Q15 observations were T=223/C=219;
non-stealing observations T=817/C=926. The supplied fixed contig carries
C, but is **not independent ground truth**. This illustrates the risk that
a new placement hypothesis can change a clinically relevant allele call.

## Scope and exceptions

{failures}

Archived real paired-end FASTQ and fixed assembled contigs were used for
both arms. This is **not** full de novo reassembly, G2P, primer-trimming or
clinical report verification; no independently adjudicated genotype/variant
truth labels are available. The 16 conditions were purposively selected,
with some repeated specimens and many ambiguous contigs.

**Recommendation:** keep this PR experimental pending orthogonal variant
truth, haplotype/mixed-infection controls and controlled performance tests.

See the scientific report and anonymized machine-readable metrics.
Generated {generation}.
"""
    tech = f"""# Real HIV/HCV remapping — scientific and engineering report

## Research question

With genuine paired FASTQ inputs and **identical fixed contigs**, how does
a bounded equal-best multi-placement remapper alter MiCall coverage, consensus,
molecular accounting, and performance compared with the unmodified mapper?

## Versions, data, and experimental controls

- Upstream MiCall base commit: **cba1898a**.
- Non-stealing candidate commit: **2bc94b59366b0ff22a8a673bdb961c28ba75487a**.
- Python 3.12 and Bowtie2 2.2.8; same machine and installed dependencies.
- Archived private 198-sample dataset; **16 selected case-conditions from
  14 distinct samples**: 9 HCV-focused and 7 HIV-focused conditions.
  Two HCV specimens were run twice with *different* preserved contig panels.
  Some HIV V3-loop libraries include HLA-B contigs.
- Each arm receives byte-identical decompressed original paired FASTQ and the
  same archived original IVA or Ivanka contig CSV. SHA-256 of each compressed
  FASTQ and contig CSV was checked between arms.
- Pipeline boundary: MiCall map_to_contigs, sam2aln, aln2counts. Not an entire
  de novo assembly, primer trimming, G2P, or downstream clinical workflow.
- Candidate uses Bowtie2 -k 8 and retains complete, equal-best paired
  alignments. The maximum-depth interpretation supplies consensus; average
  (fractional) depth is diagnostic.
- Q15 MAX consensus is read from conseq.csv (q-cutoff 15, MAX cutoff).
  Callable means the base is A/C/G/T. All comparisons are coordinate-aware
  using contig reference ID and consensus offset. MiCall's established
  coverage threshold is 100.
- Upstream per-base coverage is reconstructed from Q15 aligned.csv using
  exact contiguous A/C/G/T runs. Candidate coverage is independently
  reconstructed the same way and cross-checked against max_depth from its
  sidecar. R1/R2 QNAME flags separately verify the physical sequenced
  input-pair count; PCR duplicates are **not** collapsed.

## Study results — one row per controlled input condition

{markdown_table()}

![Callable changes](figures/callable_delta.svg)

![Alternative placements by kind](figures/ambiguity.svg)

![Comparative wall-clock times](figures/runtime.svg)

## Aggregated results and integrity checks

- Valid paired conditions: **{t['complete']}/{t['target']}**.
- Newly callable Q15 MAX case-condition positions: **{N(t['gained'])}**.
  Lost callable positions: **{N(t['lost'])}**.
  Directly conflicting base calls: **{N(t['discordant'])}**.
- Coverage decreased on **{N(t['loss_sites'])}** contig positions. These
  are not automatically consensus regressions, but require investigation.
- Site-level depth crossings at 100: {N(t['threshold_gained'])} newly
  at or above 100; {N(t['threshold_lost'])} fell below 100.
- Newly callable sites with fractional depth below 100: **{N(t['low_fraction'])}**;
  newly callable with fractional depth at least 100: **{N(t['high_fraction'])}**.
  These are *conditional* weights, not molecule-level probabilities.
- Additional alignments across completed conditions: **{N(t['placement_gain'])}**.
  Correlated placements do not increase the number of sequenced pairs.
- Independently recomputed maximum depth vs sidecar mismatches:
  **{t['sidecar_bugs']}**; unique paired-read accounting mismatches:
  **{t['pair_bugs']}**.

## Detailed HCV allele-discordance audit

For HCV-05 (S48, **16-contig archived Ivanka input**), the candidate gained
2,594 Q15 MAX callable positions and lost three, including one allele change
in contig 2-HCV-1a at the internal coordinate 1407:

| Evidence | T | C |
|---|---:|---:|
| Upstream, Q15 | 223 | 219 |
| Non-stealing, Q15 | 817 | 926 |

The candidate maximum depth at this position was **1,745**, while fractional
depth was approximately **639.47**. The input contig itself carried C.
That agreement is an internal consistency indication, **not independent
evidence of a true C allele**. The shift is highly relevant to clinical
minority-variant claims and cannot be signed off from coverage alone.

## Negative controls and unusual results

- Single-contig original IVA controls HCV-03 and HCV-04 retained identical
  Q15 MAX consensus. The candidate still added alternative alignments
  within repeated regions of the same contig, demonstrating that not all
  additional records represent *cross-contig* stealing recovery.
- The bound of eight equally best placements can truncate genuine ties.
  A contig with fewer than eight alternatives does not prove the candidate
  found every possible compatible alignment.
- **Upstream parser exception:** the HCV-01 upstream arm failed in sam2aln
  because a malformed SAM-derived CSV row had an absent reference field;
  the candidate arm completed. This is an upstream baseline failure, **not
  evidence of a comparative consensus benefit**, and is excluded from
  the paired totals. The malformed record is retained privately for audit.
- All successful A/B cases used identical input hashes. The candidate
  per-base max-depth output was independently reconstructed for validation.

## Performance and memory cautions

Time observations shown are elapsed seconds with up to seven A/B workers
sharing the same host. Contention, operating-system cache and workload
heterogeneity can dominate wall times; no meaningful p-value or
repeatability claim is possible. Additional alignments scale downstream
processing and CSV volume roughly with retained placements. The compact
sparse/packed depth accumulator limits memory overhead but does not remove
the compute/I/O cost of repeatedly considering a read.

A publication-quality runtime comparison requires serial randomized
replicates, controlled CPU affinity and memory profiling.

## Important limitations

1. **No independent biological truth.** Archived contigs and upstream
   MiCall calls are not validated patient genotypes; reference agreement
   or increased coverage cannot prove correctness.
2. **Fixed contig experiments.** These tests control remapping and consensus
   after assembly. They do not test upstream IVA/Ivanka assembly correctness.
3. **Nonrandom purposeful cohort.** Case-conditions are not independent
   specimens; no statistical extrapolation to 198 samples.
4. **Upper-bound semantics.** One read may count for multiple contig hypotheses.
   Maximum depth is not independent molecule depth; average depth is
   conditional on the bounded accepted placement set.
5. **No PCR/optical deduplication or phase truth.** QNAME-pair accounting
   corrects repeated alignment records, not true molecular duplication.
6. **Concurrent runtime tests.** Benchmark timings are descriptive only.
7. **MiCall version and workflow limitations.** A fixed Bowtie2 2.2.8 was
   used for A/B. The complete current production pipeline was not rerun.

## Reproducibility and privacy

Scripts beside this report are the repeatable runner, statistical analyzer,
position-level private diagnostics and report generator. The source data
come from the private archived real MiCall dataset; private alias-to-specimen
mapping, FASTQ/contigs, full per-base variant details, logs and SAM/CSV
intermediates are kept only in Git-ignored .validation-work/.
Tracked reports, charts and machine-readable cohort_metrics.csv use
**anonymous case labels** and aggregate numeric statistics only.

## Decision and follow-up

**Scientific verdict: experimental, not ready for clinical production.**
Observed improvements in remapping coverage are real as conditional
per-contig changes and have passed internal accounting checks, but the
same changes can alter MAX alleles and promote sites with fractional
depth below 100. Before merge for production use, obtain independently
genotyped controls, minority-haplotype spike-ins and known mixtures,
end-to-end full-MiCall validation, and repeated profiling.

Generated {generation}.
"""
    (OUT / "real_hiv_hcv_20261106_executive.md").write_text(executive)
    (OUT / "real_hiv_hcv_20261106_scientific.md").write_text(tech)
    cards = [("Valid A/B conditions", f"{t['complete']}/{t['target']}"),
             ("New callable positions", N(t["gained"])),
             ("Lost callable positions", N(t["lost"])),
             ("Changed nucleotides", N(t["discordant"]))]
    tiles = "".join("<article><small>" + escape(name) + "</small><strong>" +
                    escape(value) + "</strong></article>" for name, value in cards)
    table = "".join("<tr>" + "".join("<td>" + escape(value) +
                    "</td>" for value in row) + "</tr>" for row in table_rows())
    def svg_for_html(name):
        source = (FIG / name).read_text()
        return source[source.index("<svg"):source.index("</svg>") + 6]

    pics = "".join("<figure>" + svg_for_html(name) + "</figure>"
                   for name in ("callable_delta.svg", "ambiguity.svg", "runtime.svg"))
    html = f"""<!doctype html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Real HIV/HCV remapping A/B validation</title>
<style>
body{{font-family:system-ui,sans-serif;background:#f5f8fb;color:#152d43;
max-width:1150px;margin:0 auto;padding:25px;line-height:1.55}}
h1{{font-size:37px;line-height:1.2}}.hero,article,figure,.scroll{{
background:white;border:1px solid #dce5ee;border-radius:12px;padding:22px}}
.hero{{padding:34px}}.kicker{{color:#376da6;letter-spacing:.14em;
text-transform:uppercase;font-weight:bold;font-size:12px}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:14px;margin:24px 0}}
article small{{display:block;color:#637487}}article strong{{font-size:29px}}
.notice{{border-left:4px solid #ca833d;padding:14px 20px;background:#fff7ec}}
figure{{margin:24px 0}}figure svg{{width:100%;height:auto;max-height:750px}}
.scroll{{overflow:auto}}table{{border-collapse:collapse;width:100%}}
th,td{{text-align:left;padding:10px 13px;border-bottom:1px solid #e1e6ec}}
th{{background:#eaf0f5}}.footer{{font-size:13px;color:#607486}}
</style></head><body><section class="hero">
<p class="kicker">Scientific validation · MiCall</p>
<h1>Non-stealing remapping on real HIV/HCV data</h1>
<p>Archived paired sequencing reads · fixed contig panels · upstream A/B</p>
<div class="notice"><b>Experimental only.</b> Maximum depth counts
compatible hypotheses, not independent molecules. There is no
independent genotype ground truth in this cohort.</div></section>
<section class="cards">{tiles}</section>{pics}
<h2>Case-level evidence</h2><div class="scroll"><table><thead><tr>
<th>Case</th><th>Contig panel</th><th>References</th><th>Input pairs</th>
<th>Cross-contig alternatives</th><th>Callable delta</th><th>Base flips</th>
<th>Time ratio</th></tr></thead><tbody>{table}</tbody></table></div>
<h2>Interpretation</h2>
<p>Some multi-contig HCV cases recover substantial callable coverage.
However, fractional depth remains below the usual threshold at many newly
callable sites and at least one consensus allele changes. This is not
sufficient evidence for clinical suitability.</p>
<p class="footer">See real_hiv_hcv_20261106_scientific.md for complete
methods, limitations and reproducibility. Generated {generation}.</p>
</body></html>"""
    (OUT / "real_hiv_hcv_20261106_dashboard.html").write_text(html)
    print("REPORT", json.dumps(t, sort_keys=True))
    for p in sorted(OUT.glob("real_hiv_hcv_20261106_*")):
        if p.suffix in (".md", ".csv", ".html"):
            print(p.name, p.stat().st_size)


if __name__ == "__main__":
    run()
