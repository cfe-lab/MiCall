---
title: Admin Tasks for the MiCall Pipeline
description: Getting things done
---

## The MiCall Watcher ##

MiCall Watcher handles the automated processing of new MiSeq data
through the MiCall pipeline.  It periodically scans the `RAW_DATA` folder for data,
and when new data appears it interfaces with Kive to start the processing.
This folder is populated outside of MiCall:

* Runs get uploaded by the MiSeq to `RAW_DATA`.
* The `watch_runs.rb` script in that folder watches for the files to finish
    copying, and then creates a file named `needsprocessing` in the folder.

The Monitor looks for folders with this flag file, and ignores ones
without.

### Hourly scan for new MiSeq runs ###

Every hour, MiCall Watcher looks for new data to process with the following procedure.

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
two things happens.  If this folder is The Newest Run, MiCall Watcher gets/creates a MiCall 
pipeline run for each sample, all at once (see "Get/Create Run" below).  All other 
folders have their samples added to an in-memory list of "Samples That Need Processing".
   
    * Get/Create Run: MiCall Watcher looks for the existence of the required datasets on Kive 
    (by both MD5 and filename) and creates them if they don't exist.  Then, MiCall Watcher 
    looks to see if this data is already being processed through the current version of 
    MiCall.  If not, it starts the processing.  Either way, the Kive processing task
    is added to an in-memory list of "Samples In Progress".
    
### Periodic load monitoring and adding new processing jobs ###

Every 30 seconds, MiCall Watcher checks on all Samples In Progress.

* If a MiCall sample is finished, and that is the last one from a given MiSeq run to 
finish, then all results for that run are downloaded into a subfolder specific to the 
current version of MiCall (located in `Results/version_X.Y` where `X.Y` is the current 
version number) in that run's corresponding folder in `RAW_DATA`, and a `doneprocessing` 
flag is added to that subfolder.  Post-processing (e.g. uploading stuff to QAI) occurs.

* MiCall Watcher keeps a specified lower limit of samples to keep active at any given time.  
If a MiCall sample finishes processing and the number of active samples dips below 
that limit, MiCall Watcher looks at its list of Samples That Need Reprocessing and starts 
the next one, moving it from that list to Samples In Progress.

### Installing MiCall Watcher ###
Install the MiCall source code in a shared location:

    $ cd /usr/local/shared
    $ sudo git clone https://github.com/cfe-lab/MiCall.git

It should be run as a service, under its own user account, so first create the new user:

    # For CentOS:
    $ sudo useradd micall
    $ sudo passwd micall
    
    # For Ubuntu:
    $ sudo adduser micall

Log in as the micall user, and create a Python 3.6 virtual environment:

    $ cd ~
    $ python3.6 -m venv vmicall
    $ . vmicall/bin/activate
    (vmicall) $ cd /usr/local/share/MiCall
    (vmicall) $ pip install -r requirements-watcher.txt
    (vmicall) $ python micall_watcher.py --help

Look at the options you can give to the `micall_watcher.py` script when you
configure the service file in a later step.

Copy the logging configuration if you want to change any of the settings.

    $ cp micall_logging_config.py micall_logging_override.py

Read the instructions in the file, and edit the override copy. If the default
settings are fine, you don't need the override file.

Now configure the service using a systemd [service unit] configuration.
Here's an example configuration, in `/etc/systemd/system/micall_watcher.service`:

    [Unit]
    Description=micall_watcher
    
    [Service]
    ExecStart=/path/to/virtualenv/bin/python3.6 /path/to/MiCall/micall_watcher.py \
        --pipeline_version=8.0 --raw_data=/data/raw \
        --kive_server=bigbox --kive_user=micall_uploads \
        --micall_filter_quality_pipeline_id=100 --micall_main_pipeline_id=101 \
        --micall_resistance_pipeline_id=102 \
        --qai_server=smallbox --qai_user=micall_uploads
    Environment=MICALL_KIVE_PASSWORD=badexample MICALL_QAI_PASSWORD=worse
    User=micall
    
    # Allow the process to log its exit.
    KillSignal=SIGINT
    
    [Install]
    WantedBy=multi-user.target

The settings can either be given on the command line or set as
environment variables. Environment variables are a better option for
sensitive parameters like passwords, because the command line is visible to all
users. Make sure you reduce the read permissions on the `.service` file so
other users can't read it. The environment variable names are the same as the
command options, but they add a `MICALL_` prefix, if it's not already there.

Once you write the configuration file, you
have to enable and start the service. From then on, it will start automatically
when the server boots up.

    $ sudo systemctl daemon-reload
    $ sudo systemctl enable micall_watcher
    $ sudo systemctl start micall_watcher
    $ sudo systemctl status micall_watcher

[service unit]: https://www.freedesktop.org/software/systemd/man/systemd.service.html

### Ways to manipulate the order in which the MiCall Watcher processes data ###
Sometimes, you want to do unusual things. Here are a few scenarios we've run into.

#### You need to force MiCall Watcher to reprocess a given run ####
Remove the results subfolder for the current version of 
MiCall.  On the next hourly scan, MiCall Watcher will handle this 
folder as if it had never been processed through the current version of MiCall.  That 
is, if it's The Newest Run, its samples will be immediately started and added to 
Samples In Progress, and if it's not The Newest Run, its samples will be added to 
Samples That Need Processing.

#### You need MiCall Watcher to skip a folder and handle an older one first ####
Stop MiCall Watcher.  Add an `errorprocessing` flag to the run's folder in `RAW_DATA`.  
This will make MiCall Watcher's next hourly scan believe that it's failed and should be 
skipped, and MiCall Watcher will move on to the next one.  Note though that this has no 
effect on which folder is The Newest Run: even if you're skipping The Newest Run, 
the next one down does not inherit the mantle.

* DO NOT delete the `needsprocessing` flag from a folder to try and keep MiCall Watcher from 
handling it.  This will cause Conan's scripts to reupload data to that folder.  

Restart MiCall Watcher.  After all samples from this run have finished processing, 
remove the fake `errorprocessing` flag you set (no need to stop and restart MiCall Watcher).

#### You need to stop what's currently running and handle an older run ####
Stop MiCall Watcher.  In Kive, stop the processing tasks that you need to clear out, 
but don't remove them; the progress you've already made can be reused later when 
revisiting these tasks later.  Now, do the manipulations you need to do as in the
above case to make MiCall Watcher deal with your desired run first.  Restart MiCall Watcher.

As in the above case, when you are ready to process the run you previously stopped,
you can remove the fake `errorprocessing` flag you created for that run, and MiCall Watcher
will then restart those processing tasks on its next hourly scan.  Kive will be able 
reuse the progress already made when you stopped them.
