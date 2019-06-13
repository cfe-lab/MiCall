---
title: Resistance Design
description: Rules for Quality Control
---
The main rules for the resistance calls are built into the `pyvdrm` library,
but this module has several rules for deciding whether the data is good enough
to make a resistance call.

The main flow is in the `report_resistance()` function. The basespace version
passes in a list of regions to report resistance on, the Kive version just
reports all regions.

The first quality check is for each region, in the `check_coverage()` function:

  * The average coverage across the region has to be at least 100.
  * The coverage on the start and end positions has to be at least 100. 

That rule has a special case for HCV-NS5b, where it is split in two. The main
sample gets checked for positions 1-228, and the MIDI sample gets checked for
positions 231-561.

The coverage at all resistance positions is checked in `read_aminos()`. All the
known resistance positions for the mapped genotype are are checked, and at
least one of a sample's gene regions must have coverage of at least 100 on all
resistance positions, or no report will be generated. There's a special case
for HCV-NS5b where coverage of at least 100 on all the resistance positions
below 231 will generate a report, and so will coverage of at least 100 on all
the resistance positions at 231 or above. This rule is applied separately for
each HCV genotype and for HIV.

Anything that doesn't meet these rules should write a row into the
`resistance_fail.csv` file.
