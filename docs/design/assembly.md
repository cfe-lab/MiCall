---
title: De Novo Assembly
subtitle: Our understanding so far
---

So far, we've had the best assembly results with [IVA]. It's not perfect, though,
so this is a summary of what we've learned about how it works.

The best way to see the details is to change the [`iva` command] to write its
details to a log file like this:

```python
with open(os.path.join(tmp_dir, 'iva.log'), 'w') as log_file:
    run([IVA, '-vv', '--fr', joined_path, '-t', '2', iva_out_path],
        check=True,
        stdout=log_file,
        stderr=STDOUT)
```

That will show you the contig sequence as it gets assembled throughout the
process by writing its progress to the `iva.log` file.

Another option is to look in the scratch folder of a sample that was already
assembled, copy the `iva` command out of the assembly temp folder's `info.txt`
file, and add the `-vv` option when you run it again. You also have to change
the `iva_out` path to a folder that doesn't exist yet. For example, suppose this
is the command in `info.txt`:

    /git/MiCall/venv_micall/bin/iva --fr
    /data/run/micall-results/scratch/E1234_S1/assembly_ayvdjhcw/joined.fastq -t 2
    /data/run/micall-results/scratch/E1234_S1/assembly_ayvdjhcw/iva_out

You could rerun the assembly with this command to see all the assembly details:

    /git/MiCall/venv_micall/bin/iva -vv --fr
    /data/run/micall-results/scratch/E1234_S1/assembly_ayvdjhcw/joined.fastq -t 2
    /data/run/micall-results/scratch/E1234_S1/assembly_ayvdjhcw/iva_out2 > assembly.log 2>&1

Below, I'll walk through
a typical assembly, and explain what you see in the log file.

[IVA]: https://sanger-pathogens.github.io/iva
[`iva` command]: https://github.com/cfe-lab/MiCall/blob/283f9e50d6429f6a0d69714afa84b9b91862a6fa/micall/core/denovo.py#L213-L216


## Building the Seed
IVA runs [KMC], the k-mer counter to find the most common k-mer of the requested
length (usually 95). It logs the most common k-mer, and how many times it was
found.

    run kmc: ...
    Made new seed. kmer coverage 1998 and seed is ACGT...

Then it tries to extend the seed by going through each read to see if the last
95 bases match some part of the seed. If so, it stores the part of that read
that would extend the seed to the left. It does the same thing to the start of
all the reads and stores the part that would extend to the right.

Then it goes through all the possible extensions looking at the two most common
extensions to the left of each length from 50 down. If it finds an extension
that occurs at least 10 times, and at least 4 times more often than second place,
it uses it. Then it does the same thing to the right. Note that the "new length"
at the bottom of this section is actually the old length before the extensions
on left and right.

    Seed extension iteration 1
            k = 50 commonest two kmers: ('ACGTACGT...', 'CAGTACGT...') have frequency: (343, 567)
            k = 49 commonest two kmers: ('ACGTACG...', 'CAGTACG...') have frequency: (576, 609)
            k = 48 commonest two kmers: ('ACGTAC...', 'CAGTAC...') have frequency: (579, 620)
            k = 47 commonest two kmers: ('ACGTA...', 'CAGTA...') have frequency: (658, 696)
            k = 46 commonest two kmers: ('ACGT...', 'CAGT...') have frequency: (199, 1071)
            k = 50 commonest two kmers: ('ACGTACGT...', 'CAGTACGT...') have frequency: (132, 1111)
        Extend seed. new length=95. Bases added left:46. Bases added right:50
                     new seed: ACGTACGT...

It continues to extend the seed until it reaches the seed stop length (usually
720), or it stops making progress on both ends.

    Seed extension iteration 2
            ...
    Seed extension iteration 3
            ...
    Seed extension iteration 4
            ...
    Seed extension iteration 5
            ...
    Seed extension iteration 6
            ...
    Seed extension iteration 7
            k = 50 commonest two kmers: ('ACGTACGT...', 'CAGTACGT...') have frequency: (132, 1111)
            k = 50 commonest two kmers: ('ACGTACGT...', 'CAGTACGT...') have frequency: (132, 1111)
        Extend seed. new length=680. Bases added left:50. Bases added right:50
                     new seed: ACGTACGT...
        Extended seed OK.

If the seed extended to at least three quarters of the seed stop length, it's a
success and IVA will move on to the next step. If not, it records the failed
seed and starts again with the next most common k-mer as a new seed.

[KMC]: https://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about


## Extending with Read Pairs
Once the seed is big enough, it becomes a contig, and IVA switches to mapping
read pairs instead of looking for exact matches of the read ends.

    ______________________________ START ITERATION 1 ______________________________
    ---------------- iteration.1 start extension subiteration 0001 ----------------
        map reads	reads_1.fa	reads_2.fa
            map reads. Index:   smalt index -k 19 -s 11 iteration.1.1.1.map.map_index iteration.1.1.1.map.ref.fa
            map reads. Mapping: smalt map ...

The mapper is [smalt], and it runs in two steps: building an index of features
in the contig, and then mapping all the read pairs to find where they best align
with the contig. For this mapping step, the reads have to match at least half of
their bases to the contig. The most useful thing in this section is the file
name that shows which iteration, subiteration, and retry you're on.

After trying to map all the reads to the contig, IVA goes through all the
results looking for pairs that matched either end of the contig, and had a
soft-clipped section hanging off the end. It records the soft-clipped portions
in a collection at each end of the contig, and it writes all the useful read
pairs into temporary FASTA files. Useful read pairs fit one of the following:

* One read mapped, the other didn't. Hopefully, the other will map as the contig
    gets extended.
* Both reads mapped, and at least one of the pair has a soft-clipped portion off
    the end of the contig.

After going through all the mapped read pairs, IVA then tries to extend the left
end of the contig, similar to the way it extended the seed. Starting at length
100 or the longest soft-clipped region that was found, IVA counts up all the
possible extension sequences of that length. If the most common extension was
seen at least ten times, and at least four times more often than the second most
common extension, it gets selected. Otherwise, reduce the length by one, and try
again. Do the same thing for the right end.

        trying to extend left end ...
            k = 100 commonest two kmers: ('TCCTAGAG...', 'TCCTTGAG...') have frequency: (8, 33)
        trying to extend right end ...
            k = 100 commonest two kmers: ('ACGTACGT...', 'CAGTACGT...') have frequency: (40, 104)
            k = 99 commonest two kmers: ('ACGTACG...', 'CAGTACG...') have frequency: (41, 107)
            k = 98 commonest two kmers: ('ACGTAC...', 'CAGTAC...') have frequency: (41, 111)
            k = 97 commonest two kmers: ('ACGTA...', 'CAGTA...') have frequency: (31, 125)
        new left sequence: TCCTTGAG...
        new right sequence: CAGTA...
        extend contig seeded.00001	new_length:977	added_left:100	added_right:97

Now the new length really is the new length.

This step repeats up to four times, but it only tries mapping the useful read
pairs saved from the previous step.

        map reads	/.../tmp.extend.abc123/reads.0_1.fa	/.../tmp.extend.abc123/reads.0_2.fa
            map reads. Index:   smalt index -k 19 -s 11 iteration.1.1.2.map.map_index iteration.1.1.2.map.ref.fa
            map reads. Mapping: smalt map ...

You can see that the FASTA files are in a `tmp.extend` folder when IVA is
retrying the mapping with the useful read pairs from the previous step. It
retries up to five times with the same set of read pairs, or until it adds fewer
than 20 bases to the contig.

When it's finished its five retries or given up, IVA decides if the contig is
worth extending. A contig is worth extending if it made some progress in the
previous set of retries, or if it only retried once. If it's worth extending,
IVA starts a new subiteration.

    ---------------- iteration.1 start extension subiteration 0002 ----------------
        map reads	reads_1.fa	reads_2.fa
            map reads. Index:   smalt index -k 19 -s 11 iteration.1.2.1.map.map_index iteration.1.2.1.map.ref.fa
            map reads. Mapping: smalt map ...

You can see that it goes back to the full set of read pairs. It continues
extending the contig with subiterations and retries.

    ...
    ---------------- iteration.1 start extension subiteration 0003 ----------------
    ...
    ---------------- iteration.1 start extension subiteration 0004 ----------------
    ...

When it gets to subiteration 5, IVA does an extra mapping of all the read pairs
to the contig. This time, it requires the reads to match at least 90%, instead
of the 50% it uses when extending. Any reads that map as proper pairs are
discarded, because they match the contig that's already assembled. All the
remaining read pairs are used for the following subiterations. This filtering
gets repeated after every five subiterations.

    ---------------- iteration.1 start extension subiteration 0005 ----------------
        map reads	reads_1.fa	reads_2.fa
            map reads. Index:   smalt index -k 19 -s 11 /.../tmp.filter_reads.abc123/out.map_index /.../tmp.filter_reads.abc123/out.ref.fa
            map reads. Mapping: smalt map ...
        map reads	reads.subiter.5.reads_1.fa	reads.subiter.5.reads_2.fa
            map reads. Index:   smalt index -k 19 -s 11 iteration.1.5.1.map.map_index iteration.1.5.1.map.ref.fa
            map reads. Mapping: smalt map ...
    ...

This continues until there is no progress in a subiteration.

    ---------------- iteration.1 start extension subiteration 0014 ----------------
    ...
        new left sequence: 
        new right sequence: 
        extend contig seeded.00001	new_length:5839	added_left:0	added_right:0
        No bases added. Try trimming contigs

[smalt]: https://sourceforge.net/projects/smalt/

## Trimming
Now, IVA tries trimming the ends of the contig. It maps the reads and requires
90% again, then sorts, indexes, and builds a pileup with [samtools]. The pileup
reports how many reads are mapped to each position in the contig, or the
coverage. The coverage is tracked separately for forward and reverse reads. The
coverage is checked from the ends inward, until it meets the minimum coverage
(default 10) and the strand bias (minimum fraction from either direction,
default 0). The ends are trimmed and reported.

        map reads	reads_1.fa	reads_2.fa
            map reads. Index:   smalt index -k 19 -s 11 /.../tmp.trim_strand_biased_ends.abc123/out.map_index /.../tmp.trim_strand_biased_ends.abc123/out.ref.fa
            map reads. Mapping: smalt map ...
            map reads. sort:   samtools sort ...
            map reads. index:   samtools index /.../tmp.trim_strand_biased_ends.abc123/out.bam
    Trimming strand biased ends of contig seeded.00001 - good base range is 26 to 5839 from 5839 bases

Then another subiteration is tried.

    ---------------- iteration.1 start extension subiteration 0015 ----------------
    ...
        new left sequence: TCCTTGAG...
        new right sequence: 
        extend contig seeded.00001	new_length:5826	added_left:12	added_right:0
        map reads	reads_1.fa	reads_2.fa
            map reads. Index:   smalt index -k 19 -s 11 /.../tmp.trim_strand_biased_ends.abc123/out.map_index /.../tmp.trim_strand_biased_ends.abc123/out.ref.fa
            map reads. Mapping: smalt map ...
            map reads. sort:   samtools sort ...
            map reads. index:   samtools index /.../tmp.trim_strand_biased_ends.abc123/out.bam
    Trimming strand biased ends of contig seeded.00001 - good base range is 13 to 5826 from 5826 bases

If that doesn't extend at least as far as the section that was trimmed at the
end of the previous subiteration, it's not worth extending.

[samtools]: https://www.htslib.org/doc/samtools.html

## Making a New Seed
The reads are mapped against the contig with a 50% match required, and any
unmapped reads are used to make a new seed.

    _____________________________ Try making new seed _____________________________
        map reads	iteration.1.filtered_1.fa	iteration.1.filtered_2.fa
            map reads. Index:   smalt index -k 19 -s 11 /.../tmp.make_seed.abc123/out.map_index /.../tmp.make_seed.abc123/out.ref.fa
            map reads. Mapping: smalt map ...

IVA runs the k-mer counter to find the most common k-mer of the seed length, but
this time, it maps all the k-mers against the first seed and the contig. It then
chooses the most common k-mer that didn't map against the previous seed or
contig.

    run kmc: /.../tmp.run_kmc.abc123 kmc -fa -m32 -k95 ...
    First try of running kmc failed. Trying again with -m4 instead of -m32...
    run kmc: /.../tmp.run_kmc.abc123 kmc -fa -m4 -k95 ...
            map reads. Index:   smalt index -k 9 -s 1 /.../tmp.common_kmers.abc123/map.map_index /.../tmp.common_kmers.abc123/ref.fa
            map reads. Mapping: smalt map ...
    Made new seed. kmer coverage 1625 and seed is GACCTT...

The seed is extended in the same way as the first one. If the seed isn't at
least three quarters of the seed stop length, it is discarded, and a new seed
is built.

    Seed extension iteration 1
    ...
    Seed extension iteration 2
    ...
        Couldn't extend seed enough. That was attempt 1 of 10
    run kmc: ...
    First try of running kmc failed. Trying again with -m4 instead of -m32...
    run kmc: ...
            map reads. Index:   smalt index ...
            map reads. Mapping: smalt map ...
    Made new seed. kmer coverage 1554 and seed is CATG...
    Seed extension iteration 1
    ...
    Seed extension iteration 13
    ...
        Extend seed. new length=634. Bases added left:0. Bases added right:0
                     new seed: ACCTT...
        Extended seed OK.

## Extending Multiple Contigs
Now, IVA switches back to mapping reads, but onto both the contigs. It then uses
the soft clipped endings to try extending both contigs on the left and right
ends.

    ______________________________ START ITERATION 2 ______________________________
    ---------------- iteration.2 start extension subiteration 0001 ----------------
        map reads	iteration.1.filtered_1.fa	iteration.1.filtered_2.fa
            map reads. Index:   smalt index -k 19 -s 11 iteration.2.1.1.map.map_index iteration.2.1.1.map.ref.fa
            map reads. Mapping: smalt map ...
        trying to extend left end ...
            k = 100 commonest two kmers: ('TCCTAGAG...', 'TCCTTGAG...') have frequency: (8, 33)
            ...
            k = 13 commonest two kmers: ('TCCTA...', 'TCCTT...') have frequency: (2, 10)
        trying to extend right end ...
            k = 100 commonest two kmers: ('ACGTACGT...', 'CAGTACGT...') have frequency: (40, 104)
            ...
            k = 1 commonest two kmers: ('A', 'C') have frequency: (351, 551)
        new left sequence: TCCTT...
        new right sequence: 
        extend contig seeded.00001	new_length:5827	added_left:13	added_right:0
        trying to extend left end ...
            k = 100 commonest two kmers: ('TCCTAGAG...', 'TCCTTGAG...') have frequency: (8, 33)
            ...
            k = 47 commonest two kmers: ('TCCTA...', 'TCCTT...') have frequency: (3, 546)
        trying to extend right end ...
            k = 100 commonest two kmers: ('ACGTACGT...', 'CAGTACGT...') have frequency: (40, 104)
            ...
            k = 41 commonest two kmers: ('ACGTA...', 'CAGTA...') have frequency: (3, 223)
        new left sequence: TCCTT...
        new right sequence: CAGTA...
        extend contig seeded.00002	new_length:722	added_left:47	added_right:41

This continues with subiterations and sets of five retries, filtering out fully
mapped pairs after every five subiterations. When it stops making progress, IVA
trims the ends of the new contig for low coverage, and tries extending again.

If that makes no progress, it's time to see if the new contig is related to the
old contig.

## Combining Assembled Contigs
The contigs are compared to each other using [NUCmer] to find identical regions,
including reverse-complemented matches. If one contig is completely contained
inside the other, the smaller one is discarded. If two contigs have matching
ends, they are combined. If the matching ends are reverse complements of each
other, one of the contigs will be reverse complemented before combining the two.

[NUCmer]: http://mummer.sourceforge.net/manual/#nucmer

## More Seeds
The process continues with finding a seed, extending that seed, extending with
mapped reads, trimming the contigs, and combining the contigs, until it fails
to find or extend another seed.

## Reporting Contigs
IVA looks through all three reading frames in both directions for each contig,
and finds the longest open reading frame (without stop codons) for each contig.
It then names the contigs in order, with the longer open reading frames first.
It will also reverse complement contigs, if their longest open reading frame is
in the reverse direction.
