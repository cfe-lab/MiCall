---
title: MiCall
subtitle: Processing FASTQ data from an Illumina MiSeq
---

Maps all the reads from a sample against a set of reference sequences, then
stitches all the reads into consensus sequences and coverage maps.

A monitoring system regularly checks the file system for unprocessed runs,
transfers FASTQ.gz files to the cluster and executes the pipeline.

See the [list of steps and files][steps] for details of what the pipeline does.
The [admin] page describes how to look after the pipeline in Kive, and the
[local analysis] page describes how to get the docker version set up and run it
on your own data.

[steps]: /steps
[admin]: /admin
[local analysis]: /local_analysis.html

## Dual Licensing ##
Copyright (C) 2016, University of British Columbia

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, visit [gnu.org][gnu]. The source code for
this program is available from [github.com][github].

The program is also available for a fee under a more permissive license. For
example, if you want to run a changed version of the program on a network server
without publishing the changed source code, [contact us][contact] about
purchasing a license.

## Third Party Components ##
MiCall makes use of several open-source tools. Here is a list of tools with
their licenses.

Requests is distributed under the Apache 2.0 license.

Python 3 is distributed under the [Python 3 license][python].

Bowtie2 and Python-Levenshtein are distributed under the GNU General Public License (GPL).

Matplotlib is distributed under the [Matplotlib license][matplotlib].

Reportlab is distributed under the BSD license.

Pyyaml and Cutadapt are distributed under the MIT license.


[gnu]: https://www.gnu.org/licenses/
[github]: https://github.com/cfe-lab/MiCall
[contact]: mailto:micalldev@cfenet.ubc.ca
[python]: https://docs.python.org/3/license.html
[matplotlib]: https://matplotlib.org/users/license.html
