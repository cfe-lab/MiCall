---
title: Admin Tasks for the MiCall Pipeline
subtitle: Getting things done
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

### Installing base packages ###

MiCall is written in python, thus we need the following packages:

```shell
apt-get install -y python3 python3-venv git  # on Ubuntu & Debian.
```

### Installing MiCall Watcher ###

Install the MiCall source code in a shared location:

```shell
cd /usr/local/share
sudo git clone https://github.com/cfe-lab/MiCall.git
```

Create a Python virtual environment to run MiCall.

```shell
sudo python3 -m venv venv-micall
```

Configure micall logging, and then install micall package:

```shell

Copy the logging configuration if you want to change any of the settings.

```shell
sudo cp micall/utils/micall_logging_config.py micall/utils/micall_logging_override.py
sudo emacs micall/utils/micall_logging_override.py
sudo venv-micall/bin/pip install ./MiCall[watcher]
```

Read the instructions in the file, and edit the override copy. If the default
settings are fine, you don't need the override file.

Depending on how you configured the logging, you'll probably need to create a
log folder and grant access to the micall user.

```shell
sudo mkdir /var/log/micall
sudo chown micall:micall /var/log/micall
```

MiCall watcher should be run as a service, under its own user account, 
so first create the new user:

```shell
sudo useradd --system micall
sudo su micall  # switch to micall account.
. venv-micall/bin/activate # activate the virtual environment.
```

Test that everything is installed with the right permissions:

```shell
micall watcher --help
```

Look at the options you can give to the `watcher` script when you
configure the service file in the next step.

Now configure the service using a systemd [service unit] configuration.
Here's an example configuration, in `/etc/systemd/system/micall_watcher.service`:

```toml
[Unit]
Description=micall_watcher
    
[Service]
ExecStart=/usr/local/share/venv-micall/bin/micall watcher
EnvironmentFile=/etc/micall/micall.conf
User=micall
    
# Allow the process to log its exit.
KillSignal=SIGINT
    
[Install]
WantedBy=multi-user.target
```

Micall watcher accepts multiple settings which can be passed
directly as command line arguments, or as environment variables.
Environment variables are a better option for sensitive parameters like passwords,
because the command line is visible to all users.
Environment variables go in the configuration file listed in the
`EnvironmentFile=` setting. In this example, it's `/etc/micall/micall.conf`

```shell
exit # logout from "micall" account.
sudo mkdir /etc/micall
sudo emacs /etc/micall/micall.conf
sudo chmod 600 /etc/micall/micall.conf
```

Make sure you reduce the read permissions on the `.conf` file so
other users can't read it. The environment variable names are the same as the
command options, but they add a `MICALL_` prefix, if it's not already there.
To list all the available options, run `micall watcher --help`.
Below is the example config:

```shell
# This is an example of /etc/micall/micall.conf
# You can add comment lines that start with #
MICALL_KIVE_SERVER=https://kive.example.com
MICALL_KIVE_USER=kiveusername
MICALL_KIVE_PASSWORD=kivepassword

MICALL_QAI_SERVER=https://qai.example.com
MICALL_QAI_USER=qaiuser
MICALL_QAI_PASSWORD=qaipassword

MICALL_RAW_DATA=/data/raw

MICALL_MAIN_PIPELINE_ID=100
MICALL_FILTER_QUALITY_PIPELINE_ID=101
MICALL_RESISTANCE_PIPELINE_ID=102
```

Don't put the environment variables directly in the `.service` file, because
its contents are visible to all users with `systemctl show micall_watcher`.

Once you write the configuration file, you
have to enable and start the service. From then on, it will start automatically
when the server boots up.

```shell
sudo systemctl daemon-reload
sudo systemctl enable micall_watcher
sudo systemctl start micall_watcher
sudo systemctl status micall_watcher
```

If the service fails to start, look for detailed messages in the log file, in
`/var/log/syslog`, or in `/var/log/messages`.

[service unit]: https://www.freedesktop.org/software/systemd/man/systemd.service.html

### Restarting the MiCall Watcher ###
If you installed it as a service as described above, then it's easy:

```shell
sudo systemctl restart micall_watcher
```

Don't launch the `micall watcher` command on its own, or the service will run
won't know that it's running. That can end up running two copies of the watcher
process, and it gets confused.

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

### Purging Old Files
The RAW_DATA drive occasionally gets full, and we go through purging extra files
from old runs. One of the biggest sets of files is BCL files that you can find
in a run under `Data/Intensities/BaeCalls/L001/*/*.bcl`.

You can see how much space they take within a run folder:

```shell
find -name "*.bcl" -print0 | du -ch --files0-from -
```

We usually keep the last year's worth of BCL files around, so to delete all the
BCL files from before May 2022, we ran this command in the runs folder:

```shell
find */Data/Intensities/BaseCalls/L001 -name "*.bcl" -not -newer 220527_M04401_0226_000000000-K5YRD/SampleSheet.csv -print -delete
```
