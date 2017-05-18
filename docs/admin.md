# Admin Tasks for the MiCall Pipeline #
## The MiCall Monitor ##

MiCall Monitor (or just Monitor) handles the automated processing of new MiSeq data
through the MiCall pipeline.  It periodically scans the `RAW_DATA` folder for data,
and when new data appears it interfaces with Kive to start the processing.
This folder is populated outside of MiCall:

* Runs get uploaded by the MiSeq to `RAW_DATA`.
* The `watch_runs.rb` script in that folder watches for the files to finish
    copying, and then creates a file named `needsprocessing` in the folder.
* The [MiseqQCReport][] scripts upload the QC data to QAI, and then create a 
`qc_uploaded` file.

The Monitor looks for folders with both of these flag files, and ignores ones
without.

[MiseqQCReport]: https://github.com/cfe-lab/MiSeqQCReport/tree/master/modules

### Hourly scan for new MiSeq runs ###

Every hour, Monitor looks for new data to process with the following procedure.

* Scan for folders that have a `needsprocessing` flag.  Any such folders get added 
to a list (in memory) of all run folders.  This distinguishes other random stuff or 
stuff that isn't ready to go from actual run folders.  The newest of these folders 
(going by the date in the folder name, *not* the filesystem creation date) is 
earmarked as The Newest Run, regardless of any processing that has already been 
done by older versions of MiCall.  Folders are visited in reverse chronological 
order based on the date in the folder name.

* Iterate through the run folder list in memory.  For each run folder, look for a 
results subfolder corresponding to the current version of MiCall (located in
`Results/version_X.Y`, where `X.Y` is the current version number) with a 
`doneprocessing` flag in it.  If this is identified, skip this folder (regardless 
of whether or not it's The Newest Run).

* When a run folder is found that does *not* have such a results folder, one of 
two things happens.  If this folder is The Newest Run, Monitor gets/creates a MiCall 
pipeline run for each sample, all at once (see "Get/Create Run" below).  All other 
folders have their samples added to an in-memory list of "Samples That Need Processing".
   
    * Get/Create Run: Monitor looks for the existence of the required datasets on Kive 
    (by both MD5 and filename) and creates them if they don't exist.  Then, Monitor 
    looks to see if this data is already being processed through the current version of 
    MiCall.  If not, it starts the processing.  Either way, the Kive processing task
    is added to an in-memory list of "Samples In Progress".
    
### Periodic load monitoring and adding new processing jobs ###

Every 30 seconds, Monitor checks on all Samples In Progress.

* If a MiCall sample is finished, and that is the last one from a given MiSeq run to 
finish, then all results for that run are downloaded into a subfolder specific to the 
current version of MiCall (located in `Results/version_X.Y` where `X.Y` is the current 
version number) in that run's corresponding folder in `RAW_DATA`, and a `doneprocessing` 
flag is added to that subfolder.  Post-processing (e.g. uploading stuff to QAI) occurs.

* Monitor keeps a specified lower limit of samples to keep active at any given time.  
If a MiCall sample finishes processing and the number of active samples dips below 
that limit, Monitor looks at its list of Samples That Need Reprocessing and starts 
the next one, moving it from that list to Samples In Progress.

### Ways to manipulate the order in which the Monitor processes data ###
Sometimes, you want to do unusual things. Here are a few scenarios we've run into.

#### You need to force Monitor to reprocess a given run ####
First, stop Monitor.  Remove the results subfolder for the current version of 
MiCall.  Restart Monitor.  On the next hourly scan, Monitor will handle this 
folder as if it had never been processed through the current version of MiCall.  That 
is, if it's The Newest Run, its samples will be immediately started and added to 
Samples In Progress, and if it's not The Newest Run, its samples will be added to 
Samples That Need Processing.

#### You need Monitor to skip a folder and handle an older one first ####
Stop Monitor.  Add an `errorprocessing` flag to the run's folder in `RAW_DATA`.  
This will make Monitor's next hourly scan believe that it's failed and should be 
skipped, and Monitor will move on to the next one.  Note though that this has no 
effect on which folder is The Newest Run: even if you're skipping The Newest Run, 
the next one down does not inherit the mantle.

* DO NOT delete the `needsprocessing` flag from a folder to try and keep Monitor from 
handling it.  This will cause Conan's scripts to reupload data to that folder.  

Restart Monitor.  After all samples from this run have finished processing, 
remove the fake `errorprocessing` flag you set (no need to stop and restart Monitor).

#### You need to stop what's currently running and handle an older run ####
Stop Monitor.  In Kive, stop the processing tasks that you need to clear out, 
but don't remove them; the progress you've already made can be reused later when 
revisiting these tasks later.  Now, do the manipulations you need to do as in the
above case to make Monitor deal with your desired run first.  Restart Monitor.

As in the above case, when you are ready to process the run you previously stopped,
you can remove the fake `errorprocessing` flag you created for that run, and Monitor
will then restart those processing tasks on its next hourly scan.  Kive will be able 
reuse the progress already made when you stopped them.
