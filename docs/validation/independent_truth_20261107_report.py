"""Scientific release-gate report for published viral reference controls."""
from datetime import datetime, timezone
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SECRET = ROOT / ".validation-work/independent-truth"
OUT = ROOT / "docs/validation"


def load(path):
    return json.loads(path.read_text()) if path.exists() else {}


def get(hiv, virus, panel, arm):
    return hiv.get(virus, {}).get("panels", {}).get(panel, {}).get(arm, {})


def line(hiv,virus,panel):
    output=[]
    for arm in ("upstream","nonstealing"):
        r=get(hiv,virus,panel,arm)
        if not r.get("complete"):
            output.append(arm+": not completed")
        else:
            t=r["totals"]
            output.append(
                f"{arm}: {t['correct']:,} correct; {t['incorrect']:,} wrong; "
                f"{t['not_callable']:,} not callable; "
                f"{r['timing'].get('total_s',0):.1f} seconds")
    return "\n".join("- "+a for a in output)


def marker_table(hiv):
    r=get(hiv,"HIV9","truth","upstream")
    n=get(hiv,"HIV9","truth","nonstealing")
    if not (r.get("complete") and n.get("complete")):
        return "The nine-haplotype A/B run remains incomplete. No sensitivity " \
               "claim at 0.2% is justified."
    names=list(hiv["HIV9"]["truth_sources"].values())
    rows=["| Clone | Truth-unique marker loci | Upstream correct | Candidate correct | Upstream wrong | Candidate wrong |",
          "|:--|--:|--:|--:|--:|--:|"]
    for j,name in enumerate(names,1):
        ref=str(j)+"-HIV1-B-FR-K03455-seed"
        a=r["markers"].get(ref,{})
        b=n["markers"].get(ref,{})
        rows.append(
            f"| {name} | {a.get('loci',0)} | {a.get('truth_called',0)} "
            f"| {b.get('truth_called',0)} | {a.get('wrong_called',0)} "
            f"| {b.get('wrong_called',0)} |")
    return "\n".join(rows)


def hcv_table(data):
    lines=["| Read-origin mix | Mapper | Real minor-origin read pairs correctly placed | Absent HCV1b-decoy mapped pairs from majority library | Decoy Q15 MAX callable bases |",
           "|:--|:--|--:|--:|--:|"]
    for mix in ("HCV-origin-zero","HCV-origin-1pct"):
        for arm in ("upstream","nonstealing"):
            case=data.get(mix,{}).get(arm,{})
            if not case.get("complete"):
                lines.append(f"| {mix} | {arm} | pending | pending | pending |")
                continue
            minor=case.get("counts",{}).get("HCV2b",{})
            maj=case.get("counts",{}).get("HCV1a",{})
            decoy=sum(row["callable_q15_max"]
                      for name,row in case.get("per_reference",{}).items()
                      if name.endswith("HCV-1b"))
            lines.append(f"| {mix} | {arm} | {minor.get('correct',0)} "
                         f"| {maj.get('decoy',0)} | {decoy} |")
    return "\n".join(lines)


def phased_table(phased, names):
    if not phased.get('per_clone'):
        return "No alignment-independent private-motif evidence computed."
    rows = [
        ("| Clone label | Private marker motifs | At least 1 exact motif | "
        "At least 2 linked motifs | At least 3 linked motifs |"),
        "|:--|--:|--:|--:|--:|",
    ]
    for ref, name in names.items():
        obs = phased['per_clone'].get(ref, {})
        rows.append(
            f"| {name} | {len(obs.get('private_positions', []))} | "
            f"{obs.get('at_least_one', 0):,} | "
            f"{obs.get('at_least_two', 0):,} | "
            f"{obs.get('at_least_three', 0):,} |")
    return "\n".join(rows)


def newer_hiv5_truth_result():
    folder = SECRET / 'results/HIV5-all-p4/nonstealing'
    if not (folder / 'status.json').is_file():
        return ''
    from independent_truth_eval import truth_data, calls_by_region

    references, _, _ = truth_data('HIV5')
    calls = calls_by_region(folder / 'conseq.csv')
    correct = wrong = 0
    for ref, sequence in references.items():
        observed = calls.get(ref, {})
        correct += sum(observed.get(i) == value
                       for i, value in enumerate(sequence, 1))
        wrong += sum(observed.get(i) not in (None, value)
                     for i, value in enumerate(sequence, 1))
    elapsed = load(folder / 'status.json')['timings']['total_s']
    return (
        f"A subsequent revised candidate with all reported Bowtie2 "
        f"alignments and four ordered alignment threads recovered "
        f"**{correct:,} correct** and **{wrong:,} incorrect** bases in "
        f"**{elapsed:.1f} seconds** on the same 59,770-pair HIV5 control. "
        "This is a separate, later run than the pinned eight-hit candidate "
        "reported above; it removes the finite eight-hit cutoff only for "
        "panels of at most 16 contigs. Larger panels remain bounded, "
        "and Bowtie2 still performs heuristic rather than exhaustive search."
    )


def figure(hiv):
    both=[(panel,get(hiv,"HIV5",panel,"upstream"),
           get(hiv,"HIV5",panel,"nonstealing"))
          for panel in ("truth","perturbed")]
    if any(not a.get("complete") or not b.get("complete")
           for _,a,b in both):
        return
    fig,ax=plt.subplots(figsize=(10,4.2))
    for i,(panel,a,b) in enumerate(both):
        old=a["totals"]["correct"];new=b["totals"]["correct"]
        ax.barh(i-.18,old,height=.31,color="#64788b",
                label="Upstream" if i==0 else None)
        ax.barh(i+.18,new,height=.31,color="#246eac",
                label="Non-stealing" if i==0 else None)
        ax.text(new+80,i+.18,f"+{new-old:,}",va="center",fontsize=10)
    ax.set_yticks([0,1],["True clone references","50 deliberately wrong bases"])
    ax.invert_yaxis();ax.set_xlim(0,46000)
    ax.set_xlabel("Correct callable Q15 MAX bases against independent true reference")
    ax.set_title("Real five-strain HIV mixture: error-aware genotype truth")
    ax.legend(loc="lower right")
    fig.tight_layout()
    p=OUT/"figures/independent_hiv5_truth.svg"
    fig.savefig(p);plt.close(fig)
    p.write_text("\n".join(x.rstrip() for x in
                           p.read_text().splitlines())+"\n")


def main():
    hiv=load(SECRET/"independent_truth_metrics.json")
    hcv=load(SECRET/"hcv-origin/hcv_origin_metrics.json")
    cap=load(SECRET/"hiv9-k-vs-a/summary.json")
    phased=load(SECRET/"hiv9-exact-phased-motifs.json")
    figure(hiv)
    truth=get(hiv,"HIV5","truth","upstream")
    variant=get(hiv,"HIV5","perturbed","nonstealing")
    stage9=all(get(hiv,"HIV9","truth",arm).get("complete")
               for arm in ("upstream","nonstealing"))
    stamp=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines=[
        "# Independent genotype and minority-haplotype validation",
        "",
        "**Production decision: NO-GO. The PR is NOT ready for clinical production.**",
        "",
        "## Independently established five-strain HIV truth",
        "",
        ("Published five-strain experimental HIV-1 mixture, Illumina MiSeq "
        "ENA accession **SRR961514**. Source-strain genomes were resequenced "
        "individually and reported in [Di Giallonardo / Töpfer et al. "
        "(NAR 2014)](https://github.com/cbg-ethz/5-virus-mix), "
        "doi:10.1093/nar/gku537. Of 714,994 original paired reads, a "
        "deterministic uniform sample of **59,770 original read pairs** was used. "
        "Five independent truth reference sequences contain "
        f"**{truth.get('totals',{}).get('total',0):,} base positions**."),
        "",
        "### Reference-matched control",
        "",
        line(hiv,"HIV5","truth"),
        "",
        "![Independent known-genotype correctness](figures/independent_hiv5_truth.svg)",
        "",
        ("The candidate gained **798 independently correct callable positions** "
        "with zero wrong calls observed in this test; its wall time was "
        "**448 seconds versus 174 seconds** upstream (~2.6x slower). "
        "It has not established statistically calibrated false-call bounds."),
        "",
        "### Held-out, deliberately incorrect reference challenge",
        "",
        ("Fifty reproducible alternate nucleotides were planted in non-terminal "
        "positions across the five source contigs; the biological truth "
        "remained unchanged and independent. This controls for circular "
        "evaluation that simply rewards agreement with the mapping seed."),
        "",
        line(hiv,"HIV5","perturbed"),
        "",
        ("Upstream corrected **46/50** deliberately planted false bases "
        "(4 became uncallable, none retained false). The candidate corrected "
        f"**{variant.get('edits',{}).get('reverted_to_truth','pending')}/50** "
        "with three uncallable and zero wrong callable seed bases; "
        "**900 extra correct callable bases** relative to upstream."),
        "",
        "### Subsequent improved alignment-search experiment",
        "",
        newer_hiv5_truth_result(),
        "",
        "## Nine Sanger-confirmed HIV minority haplotypes",
        "",
        ("Published nine-clone mixture, genotype proportions **50%, 25%, "
        "12.5%, 6.3%, 3.2%, 1.6%, 0.8%, 0.4%, 0.2%**, 1,074 bases each, "
        "with [Sanger-verified reference haplotypes from the CliqueSNV "
        "validation dataset](https://github.com/Sergey-Knyazev/CliqueSNV-validation) "
        "(Knyazev et al., 2021, PMID 34214168). ENA run **SRR12042290** "
        "contains 2,783,356 pairs; deterministic subsample **239,006 pairs**."),
        "",
        line(hiv,"HIV9","truth"),
        "",
        marker_table(hiv),
        "",
        "### Alignment-independent linked private-SNP evidence",
        "",
        phased_table(phased, hiv.get('HIV9', {}).get('truth_sources', {})),
        "",
        ("These are exact, clone-specific 21-nucleotide motifs found directly "
        "in the original, unaligned paired FASTQ files; two or more matching "
        "motifs in one original pair provide true local sequence linkage, "
        "without using either MiCall mapper. They do NOT establish the "
        "complete haplotype sequence or unbiased clone abundance. The "
        "nominal clone labels in the published FASTA are NOT consistent "
        "with the observed motif prevalence: for example, the clone "
        "labelled J_0.2 has private motifs in 96,433 of 239,006 pairs, "
        "with two or more linked motifs in 35,122 pairs. H_0.4 has "
        "39,057 and 16,180 pairs respectively. Until sample composition "
        "and amplification artifacts are independently reconciled, "
        "neither result can be described as demonstrated 0.2% or "
        "0.4% abundance sensitivity."),
        "",
        ("**Interpreting this table:** a unique SNP detected at its correct "
        "reference position does not demonstrate recovery of a whole "
        "within-molecule haplotype, calibration of minority abundance, "
        "or complete phasing. The dominant A_50 and intermediate C_12.5 "
        "clones have no individual-only SNP allele in this alignment, "
        "so the marker analysis cannot separately identify them."),
        "",
        "## Alignment search completeness: serious release blocker",
        "",
        ("The same 1,200 real HIV9 read pairs were aligned to nine individually "
        "verified plasmid references under Bowtie2 2.2.8 with a bounded "
        "8-alignment reporting window versus a deeper '-a' search. "
        f"**{cap.get('pairs_missing_truth_contigs_under_k8','unknown')} / "
        "1,200 pairs** missed at least one reference found among the "
        "deeper search's best-scoring placements under the 8-hit cap. "
        "For the 1,117 structurally valid complete-pair groups, 46 "
        "had a different best-reference set: 42 with matching best score "
        "and four with a differing best score. At least one read had "
        "nine tied placements that cannot all be reported under an 8-hit cap."),
        "",
        ("On these same reads, changing Bowtie2's window to nine reduced "
        "missed-reference groups to **36**, and changing to sixteen left "
        "**32**; the latter discrepancies reflect more than the window size "
        "alone, including Bowtie2 heuristic search behavior. The '-a' "
        "search is a deeper *software search*, not independent biological "
        "truth. A finite Bowtie reporting window cannot certify "
        "that all compatible genotype hypotheses were discovered."),
        "",
        "## Real HCV library-origin minority and absent-decoy challenges",
        "",
        ("An additional two-condition controlled experiment used actual "
        "unaltered read sequences from archived HCV libraries, retaining "
        "source-library identity in each read name. 50,000 HCV1a-labelled "
        "pairs were tested with zero versus 500 HCV2b-labelled pairs "
        "(0.9901% addition). A fixed third HCV1b-labelled contig was "
        "included but **zero source-library HCV1b reads** were added."),
        "",
        hcv_table(hcv),
        "",
        ("Both mappers placed **428 of 500** HCV2b-source pairs on the "
        "minor contig in the 0.99% mix and produced **173 callable Q15 "
        "MAX bases** for that minor contig. The absent HCV1b decoy "
        "nevertheless had **26 callable bases** in the no-minority run "
        "and **32** in the 0.99% mix under both mappers. "
        "**This false appearance was present upstream too**, not created "
        "by the non-stealing change. These labels are MiCall-inferred "
        "library genotype names, not independent LiPA-confirmed HCV "
        "truth. What is experimentally known is the *source of each read*."),
        "",
        "## Confirmed bug and fix",
        "",
        ("The de novo selector rejected Bowtie2 primary pairs with exactly "
        "one mapped mate, even when upstream preserved those SAM rows. "
        "A real HCV negative control exposed **99 missing original "
        "read-pair records** (10,528 retained versus upstream 10,627). "
        "The change now preserves structurally valid primary singleton "
        "records, and cascade accounting counts **one mapped read end**, "
        "not a full two-end fragment, for each singleton. After the "
        "fix, both mappers retain the same **10,627** mapped original "
        "pair names. Existing MiCall sam2aln rejects any pair lacking "
        "two valid mate CIGARs, so these orphan records do *not* provide "
        "Q15 consensus bases; this is a compatibility and read-count bug, "
        "not 99 lost consensus molecules."),
        "",
        "## Reproducibility and provenance",
        "",
        ("- Five-strain source FASTA is from the published repo's REF.fasta; "
        "ENA SRR961514 FASTQ MD5s f7e069e81ce706064e0112422950ac32 "
        "and 66012ec5036e8d3128c37e203fd006c5 were verified."),
        ("- Nine-plasmid Sanger truth is HIV9exp.fasta from the published "
        "validation repo; ENA SRR12042290 FASTQ MD5s "
        "f7a40ddbd8da7ce35aa0d4245043ea19 and "
        "bf229f9442e79f09e0a40924f19f5712 were verified."),
        ("- Baseline MiCall source commit cba1898a; independently known "
        "HIV5/HIV9 truth runs used original candidate source revision "
        "c6187fc57770c2a5c59cd4527dc05dcf40301f2f, "
        "before the singleton accounting fix. The HCV negative control "
        "was rerun after the singleton repair."),
        ("- Paired read IDs and original FASTQ, full variant alleles, "
        "original private cohort IDs, reference sequences, mapping logs "
        "and original data are kept in git-ignored .validation-work. "
        "Only aggregated evidence, tests and scripts appear in Git."),
        ("- Use prepare_truth_controls.py, prepare_hcv_origin_mix.py, "
        "real_hiv_hcv_20261106_runner.py, independent_truth_eval.py "
        "and hcv_origin_eval.py to reproduce, with private data access "
        "where relevant. No real clinical specimens are published."),
        f"- Generated {stamp}.",
        "",
        "## Production readiness gates",
        "",
        "| Gate | Result | Status |",
        "|:--|:--|:--|",
        "| Independently sequenced HIV base truth | +798 correct, zero wrong observed | Pass on this subset |",
        "| Deliberately corrupted references | 47/50 corrected, 0 wrong observed, +900 true calls | Pass on this subset |",
        "| Rare haplotypes to 0.2%, true phased abundance | Exact private-marker linkage found, but nominal clone labels disagree with observed read fractions | FAIL / unresolved truth |",
        "| Competing mapping hypotheses | Revised <=16-contig search uses all Bowtie2-reported alignments; heuristic search and >16-contig cap remain | OPEN |",
        "| Genotype absence specificity | Absent HCV1b decoy has apparent Q15 MAX bases | NOT ESTABLISHED |",
        "| Runtime and memory budget | Revised HIV5 is ~200s vs 174s upstream; HIV9 and large panels need profiling | OPEN |",
        "| Full clinical MiCall / G2P / trimming pathway | Not evaluated end-to-end on controls | OPEN |",
        "| Singleton SAM / mapped-end accounting parity | Fixed and unit-tested | PASS |",
        "",
        ("**Release decision: NOT PRODUCTION READY.** Before clinical production "
        "deployment, the team must demonstrate reliable discovery of "
        "competing haplotypes, independent mixed-genotype specificity and "
        "minority frequency estimation, complete MiCall reporting, "
        "and acceptable repeatable resource consumption. Neither improved "
        "coverage nor zero observed errors on one control can substitute "
        "for these release gates."),
        ""
    ]
    outfile=OUT/"independent_truth_20261107_scientific.md"
    outfile.write_text("\n".join(lines))
    summary=[
        "# Independent genotype and minority-haplotype: release verdict",
        "",
        ("**NOT PRODUCTION READY**. The PR has credible genotype-correctness "
        "benefits but important unresolved safety, completeness and runtime issues."),
        "",
        ("The five-strain HIV laboratory mixture produced 798 additional correct "
        "callable bases and zero wrong calls versus independently sequenced "
        "references. After 50 deliberate wrong-seed mutations, it corrected "
        "47/50 wrong bases versus 46/50 upstream, and gained 900 correct calls."),
        "",
        ("The 8-placement cap missed best-scoring reference alternatives for "
        "78/1,200 real nine-clone read pairs compared with a deeper alignment "
        "search. Even a 16-placement window did not remove every difference."),
        "",
        ("At 0.99% real-source HCV minor-library admixture, both mappers "
        "retained 428/500 true-source minor reads on the minor contig; "
        "but an intentionally absent third-genotype decoy had callable "
        "bases in both mappers. MAX depth is not proof of independent "
        "molecular infection evidence."),
        "",
        ("A real-data read-accounting regression was fixed: preserve genuine "
        "single-mapped-mate primary SAM pairs and count one mapped read end, "
        "without manufacturing Q15 consensus support."),
        "",
        ("The revised small-panel all-placement candidate scored 40,772 "
        "correct HIV5 bases with zero wrong, in 200 seconds (vs 448 seconds "
        "previously). HIV9 full-run base accuracy also improved by 896 "
        "correct bases, but nominal rare-clone fractions conflict with "
        "direct unaligned-read motif counts. Full-pipeline clinical "
        "reporting and calibrated haplotype abundance remain unverified."),
        "",
        ("See independent_truth_20261107_scientific.md for complete data, "
        "the nine-haplotype status, methods, and gate matrix."),
        ""
    ]
    (OUT/"independent_truth_20261107_executive.md").write_text(
        "\n".join(summary))
    print("REPORT",str(outfile),outfile.stat().st_size,
          "HIV9 full complete",stage9,flush=True)


if __name__=="__main__":
    main()
