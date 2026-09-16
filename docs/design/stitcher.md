---
title: Contig Stitching in MiCall
---

MiCall performs post-assembly stitching after de novo assembly.
This document explains why there are two stitchers and where each
detailed design lives:

* [Referencefull stitcher](referencefull_stitcher.md)
* [Referenceless stitcher](referenceless_stitcher.md)

## The shared problem

De novo assemblers such as IVA or Haploflow do not always turn input
reads into a single contiguous sequence. They may return multiple
contigs: fragmented sequences that can represent adjacent or
overlapping portions of the same biological sequence, sometimes
encoding the same region more than once.

A post-assembly stitching/refinement step can therefore improve the
assembly by systematically arranging those contigs and resolving
discrepancies in overlapping regions.

MiCall has two different ways to perform that refinement.

## Referencefull stitcher

The referencefull stitcher
(`micall/utils/referencefull_contig_stitcher.py`,
`micall contig_stitcher with-references`)
is allowed to use a reference sequence and positions relative to that
reference.

That gives it valuable information unavailable from the contigs
alone. In particular, it can reason about ordering and adjacency
using the reference coordinate system and can often produce
substantially more complete assemblies.

This is useful and intentional.

However, its output is therefore reference-guided. If the biological
sequence genuinely differs structurally from the reference — for
example through a large deletion, inversion, rearrangement,
duplication, or other scrambled structure — reference-derived
ordering can potentially obscure or normalize that structure.

## Referenceless stitcher

The referenceless stitcher
(`micall/utils/referenceless_contig_stitcher.py`,
`micall contig_stitcher without-references`)
exists for a different purpose.

Its goal is:

> Improve the initial de novo assembly without using a reference
> sequence or reference-derived structural assumptions.

It relies only on evidence intrinsic to the sample, such as:

* the contig sequences themselves;
* sequence overlap between contigs;
* short-read evidence;
* other sample-derived linkage evidence actually available to the
  implementation.

Reference independence is a deliberate property of this output, not
an implementation deficiency.

The important use case is that a user may want an improved de novo
assembly while still preserving unusual biological structure exactly
as supported by the sample. Note that referenceless refinement is a
step applied after de novo assembly; "unstitched" output is not
automatically reference-free in the same sense.

## Why both exist

Neither stitcher is simply "better" than the other. They answer
different questions:

```text
referencefull:
    What assembly can we obtain when reference-derived structure
    is allowed as evidence?

referenceless:
    What improvement over the de novo assembly can be justified
    from sample-intrinsic evidence alone?
```

The referencefull result can be more complete. The referenceless
result has stronger reference-independence semantics. These are
complementary products.

Conceptually:

```text
de novo contigs (e.g. IVA, Haploflow)
        |
        +-- referencefull refinement -------> more complete,
        |    (reference + contigs + reads)     reference-guided assembly
        |
        +-- referenceless refinement --------> improved assembly with
             (contigs + reads,                 reference-independent
              no reference structure)         structure preserved
```

See the detailed designs for the algorithm behind each product:

* [Referencefull stitcher](referencefull_stitcher.md)
* [Referenceless stitcher](referenceless_stitcher.md)
