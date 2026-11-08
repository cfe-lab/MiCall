# Real HIV9 fixed-contig 17-reference remapping comparison

**This is a constructed, 17-authentic-reference competition panel, not a 17-contig clinical specimen or a full pipeline validation.**

The same 1,200 real paired HIV9 reads and 17 fixed authentic references (nine Sanger-verified HIV clones, five independently sequenced HIV strains and three archived HCV contigs) were mapped with two versions of the candidate. The former version capped large panels at -k 8; the new version permits all Bowtie2 reported hits for 17–32 references whose total length is at most 65,536 nt. FASTQ and contig SHA-256 checksums matched.

## Algorithm revision provenance

The previous branch was at ec50f74c and used the fixed eight-hit search for 17 contigs. The revised mapping source was subsequently committed unchanged as fad064abb503774103b9ac16f810c762285de97f. The private runner labels the two execution environments search-policy-before and search-policy-after. The exact original FASTQ and reference CSV SHA-256 values were compared locally but are not published to avoid linking private inputs. These runs were concurrent and not a controlled time trial.

## Independently sequenced HIV9 clone truth

| Clone | Before correct | After correct | Delta | Before wrong | After wrong |
|:--|--:|--:|--:|--:|--:|
| A_50 | 69 | 2 | -67 | 0 | 0 |
| B_25 | 0 | 0 | +0 | 0 | 0 |
| C_12.5 | 0 | 0 | +0 | 0 | 0 |
| D_6.3 | 0 | 0 | +0 | 0 | 0 |
| E_3.2 | 0 | 0 | +0 | 0 | 0 |
| F_1.6 | 0 | 0 | +0 | 0 | 0 |
| G_0.8 | 0 | 0 | +0 | 0 | 0 |
| H_0.4 | 558 | 633 | +75 | 0 | 0 |
| J_0.2 | 894 | 902 | +8 | 0 | 0 |
- **True HIV9 correct callable bases:** 1,521 before versus 1,537 after.
- **Incorrect callable bases:** 0 before versus 0 after.

## Alignment accounting and observed runtime

| Metric | Previous candidate | Resource-aware candidate |
|:--|--:|--:|
| Final mapped-mate placement rows | 3,828 | 4,210 |
| Saturated original read-pair groups | 1,167 | 0 |
| Total wall seconds | 10.54 | 16.19 |

**Nuance:** total correct bases improved but not monotonically on each clone. The best-sequence set changes when additional alignment alternatives are searched; the 1% and 0.2% sample frequencies in the published HIV9 benchmark are not independently established for this real-read subsample.

No variant correctness claim beyond the Sanger-verified contig-level truth is supported. Bowtie2 all-hit search is heuristic, and the larger-than-32 or >65,536 nt fallback retains a finite -k 8 limit. The 30-minute processing budget has not been established on broad 198-sample production panels.

**Clinical production decision: NO-GO.**
