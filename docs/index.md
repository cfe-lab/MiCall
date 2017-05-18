# MiCall #
## Processing FASTQ data from an Illumina MiSeq ##
Maps all the reads from a sample against a set of reference sequences, then
stitches all the reads into consensus sequences and coverage maps.

A monitoring system regularly checks the file system for unprocessed runs,
transfers FASTQ.gz files to the cluster and executes the pipeline.

See the [list of steps and files][steps] for details of what the pipeline does.
The [admin] page describes how to look after the pipeline.

[steps]: http://cfe-lab.github.io/MiCall/steps
[admin]: http://cfe-lab.github.io/MiCall/admin

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

[gnu]: http://www.gnu.org/licenses/
[github]: https://github.com/cfe-lab/MiCall
[contact]: http://www.google.com/recaptcha/mailhide/d?k=01yIEKHPqcNIG1ecF-Khum4g==&c=SnES_swjOvb0d22WN4Q30vx9gKyzHNDkTxbY6Y1os4w=

