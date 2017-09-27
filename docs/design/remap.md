---
title: Remap Design
description: How seeds are chosen
---

Choose which seeds to use as follows:
1. Count all the reads that the prelim_map step mapped to all the seeds,
  skipping any shorter than 50 bases, because they are often primers.
2. Look at how many reads mapped to each seed, and drop any that mapped
  fewer than ten.
3. If more than one seed from the same seed group (HCV-1a, 1b, 2, 3, 4, etc.)
  pass the threshold, select the one with the most reads. In case of a tie,
  choose the first seed, alphabetically.
4. Each iteration of remapping, look at the consensus sequence for each seed.
  If it has drifted closer to one of the other seeds than to its original seed,
  drop it. Drop one at a time, and recalculate distances, because two seeds can
  drift toward each other and both be candidates for dropping.
