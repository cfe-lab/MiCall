---
title: Genome Coverage Maps
subtitle: Visualizing assembled sequences
---
When MiCall uses a denovo assembler, it can be confusing to figure out how the
assembled sequence relates to known references like HIV HXB2. The genome
coverage maps show which parts of each assembled contig map to which parts of
the known reference sequences.

Here's a basic example of a simple contig that was assembled from reads in the
NS5b region of Hepatitis C:

[![basic contig]][basic contig]

The coloured bars across the top represent the gene regions in a reference
sequence. The grey bars below represent the consensus sequence and the assembled
contig that the reads were mapped to and the consensus was built from. The
labels on the contig and consensus sequence show what reference was the best
BLAST match, as well as the *maximum* depth of reads that mapped to the contig.

The wavy green or blue region above the grey bars shows how the read depth
varies along a contig. The colour changes with depth, and darker blue means
deeper coverage.

Finally, the two arrows show where each part of the contig best matches the
known reference, according to the minimap2 search tool. MiCall will do its best to
line up the contig with its matching region in the reference, like this example,
but some of the more complicated sequences need these arrows to understand which
parts match which.

Often, a sample will assemble more than one contig, either because there are
multiple viruses in the same sample, or because the assembly is incomplete. This
example shows Hepatitis C and two separate contigs from HIV:

[![multiple contigs]][multiple contigs]

You can see the Hepatitis C reference across the top, followed by its contig.
Then the HIV reference starts the bottom section. The contig on the left is a
Nextera reaction in the RT gene that had to be assembled from fragments, and the
contig on the right is an amplicon reaction in V3LOOP that was detected because
there were many reads of the same length. (Notice the rectangular coverage.)
Below contig 3, there is a smaller consensus sequence without a contig. That
is a special case for MiCall's G2P analysis of V3LOOP sequences that do not get
assembled or mapped.

One other thing to note about this diagram is that the contigs vary in coverage
depth between 100 and 150 reads, but the coverage maps are all the same height.
The coverage maps show coverage changes within each contig, but can't be used to
compare coverage between contigs. The label on each contig shows the maximum
depth, and the colour changes with the depth.

Sometimes, a genome will include reversed sections. When that happens, the
direction of the arrows will show which way each section matched the reference.
Here's an example where minimap2 found three separate sections in a single contig:

[![inversion]][inversion]

Section 1.1 matched 5' LTR and gag in the forward direction, section 2.2 matched
RT in the reverse direction, and section 1.3 matched 3' LTR in the forward
direction. There are also some small, yellow sections on either side of section
1.2. Yellow means that section didn't match anything in the reference. You will
also sometimes see green sections that represent insertions and blank sections
that represent deletions.

Sometimes, a contig matches two sections of the reference, but minimap2 can't
tell exactly where the boundary is between the two sections. This happens when
the end of one section in the reference is an exact match for the start of the
other section. When this happens, the minimap2 matches overlap, and are
stacked on top of each other, like 1.1 and 1.2 in this diagram:

[![overlap]][overlap]

Because it's impossible to tell which match the overlapped bases belong to,
MiCall arbitrarily assigns the coverage to the one that comes first in the
contig. In this example, the two bases in gag would get the coverage, not the
two bases in nef.

Finally, there are often extra contigs that don't match any references, or only
match small sections within the contig. Those are displayed after the known
references as partial matches. If only a small part of a contig matched HIV, for
example, the contig would be labelled "HIV-partial". This example didn't match
any reference, so it's labelled "unknown-partial":

[![partial]][partial]

The grey bars at the top of the partial results section are each 500 bases long,
so you can estimate the length of each contig in that section.

[basic contig]: images/genome_coverage_basic.svg
[multiple contigs]: images/genome_coverage_multiple_contigs.svg
[inversion]: images/genome_coverage_inversion.svg
[overlap]: images/genome_coverage_hits_overlap.svg
[partial]: images/genome_coverage_partial.svg