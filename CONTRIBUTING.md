# Contributing to the MiCall Project

If you like this project and want to make it better, help out. You could report
a bug, or pitch in with some development work.

## Bug Reports and Enhancement Requests

Please create issue descriptions [on GitHub][issues]. Be as specific as possible.
Which version are you using? What did you do? What did you expect to happen? Are
you planning to submit your own fix in a pull request?

[issues]: https://github.com/cfe-lab/MiCall/issues

## Development
To see how all the tools should be installed, follow the steps in `Dockerfile`
and `dev.dockerfile`. If you prefer, you can run your development environment
under docker, as described in `dev.dockerfile`. The same installation steps are
also listed in the `Singularity` file. The docker file runs with debian, and the
singularity file runs with CentOS.

If you want to see what's currently being worked on, check out the active tasks
in our [milestones].

[milestones]: https://github.com/cfe-lab/MiCall/milestones

### Python
Check that Python is already installed.

    python --version

We have tested with Python 3.8.

### BaseSpace
Set up the [native apps virtual machine][bsvm], and configure a shared folder
called MiCall that points to the source code. To get the shared folder working,
you'll probably need to update the VBox guest additions and add the basespace
user to the `vboxsf` group. Make sure you have a developer account on
illumina.com.

Use the `docker_build.py` script to build a Docker image and push it to
BaseSpace. If you add `-t vX.Y`, it will add a tag to the Docker image. If you
add `-a <agent id>`, it will launch the spacedock tool to process samples as a
local agent. You can also set the `BASESPACE_AGENT_ID` environment variable so
you don't have to supply it every time. You can get the agent id from the Form
Builder page on BaseSpace.

    sudo /media/sf_MiCall/docker_build.py -a abcde12345

[bsvm]: https://developer.basespace.illumina.com/docs/content/documentation/native-apps/setup-dev-environment

### Test data
If you want to run `micall_watcher.py`, you have to set up data folders for raw
data and for the working folders. You'll also need to set up the QAI project
and the MiseqQCReport so you can download QC data and upload results.

1. Create a data folder somewhere on your workstation, like ~/data. Create
   subdirectories called miseq and RAW_DATA. Add folders RAW_DATA/MiSeq/runs.
2. Connect to the shared drive [using CIFS][cifs] and mount
   smb://192.168.68.144/RAW_DATA as /media/RAW_DATA.
3. Navigate down to /media/RAW_DATA/MiSeq/runs, pick a recent folder, and make
   sure it has a file named needsprocessing.
4. Copy SampleSheet.csv to a sample run folder under your local
   RAW_DATA/MiSeq/runs folder.
5. Navigate down to Data\Intensities\BaseCalls, and copy a few of the .fastq.gz
   files to your sample run folder under Data/Intensities/BaseCalls.
5. Copy the Interop folder and the files RunInfo.xml and runParameters.xml.
6. Open settings.py for editing.
7. Point `home` at your local data/miseq folder.
8. Point `rawdata_mount` at your local RAW_DATA folder.
9. Set the Oracle connection information to a test database where you can upload
   sequence data. (https://git-int.cfenet.ubc.ca/wscott/oracleserver)
10. Run the Ruby console for QAI and `LabMiseqRun.import('01-Jan-2000')` for the
    date of your sample run.
11. Upload the projects to a micall pipelines in QAI, use `micall.utils.projects_upload` to create a new pipeline in QAI
11. Run micall_watcher.py, it does need arguments. Look up the container app ids from Kive, check the Kive server URL and ports as well as QAI server and port

[cifs]: https://wiki.ubuntu.com/MountWindowsSharesPermanently

### Looking at SAM files
When you don't understand the pipeline's output, it can be helpful to look at
the raw reads in a sequence viewer like [Tablet][tablet]. Run the `micall_docker`
script on a run folder or a single sample, like this:

    python micall_docker.py folder --debug_remap --all_projects --keep_scratch /path/to/run

The options tell it to write the debug files, use all projects, and save the
scratch folder that holds all the debug files. Look through the scratch folders
under the run folder to find
the one for the sample you're interested in. The remap step writes the mapping
results as `debug_remapX_debug.sam` and `debug_remapX_debug_ref.fasta`, where
`X` is the remapping iteration number. You should be able to open an assembly
in Tablet using those two files. If the SAM file contains multiple regions,
you'll probably have to sort it with the `micall/utils/sort_sam.py` script. That
same script can convert `prelim.csv` into a SAM file.

If you want to understand the de novo assembly process, read through the
[assembly page].

[tablet]: https://ics.hutton.ac.uk/tablet/
[assembly page]: https://cfe-lab.github.io/MiCall/design/assembly.html

### GitHub Web Site
Most of the time, you can change the web site content just by editing the
markdown files in the `docs` folder. However, you may occasionally need to dig
into the page templates or do more serious work. If that happens, you can test
out the web site locally before publishing it.

1. Install Ruby 2.6, preferably with [Ruby Version Manager].

       rvm install 2.6
       rvm use 2.6

2. Install the gems for the web site.

       cd MiCall/docs
       gem install bundler
       bundle install

3. Serve the web site.

       bundle exec jekyll serve

What changes might you want to make? The web site is based on the
[Bulma Clean Theme], so read through the documentation there to see if it
already has the feature you want. Usually, the advanced features require you
to write files in the `docs/_data` folder or add settings to the front matter
at the top of a markdown file.

If you have to add a new feature to the web site, you can override one of the
files in the theme by copying it into the `docs/_includes` folder, and making
changes there. Consider offering it back to the theme project as a pull request,
because any files you override won't get automatic improvements from the
original theme project.

[Ruby Version Manager]: https://rvm.io/rvm/install
[Bulma Clean Theme]: https://github.com/chrisrhymes/bulma-clean-theme

### Releases
This section assumes you already have a working server up and running, and you
just want to publish a new release. If you're setting up a new server, follow
similar steps to setting up a development workstation. Follow these steps:

1. Check that all the issues in the current milestone are closed, and make sure
    the code works in your development environment. Run all the unit
    tests as described above, process the microtest data set with
    `release_test_microtest.py`.
2. Check if the kiveapi package needs a new release by looking for new commits.
    Make sure you tested with the latest version.
3. Determine what version number should be used next.
4. Use the `projects_dump.py` script for the previous version and compare
    `projects.json` to check that the projects match, or that the differences
    were intended. Test the `projects_upload.py` script with your updated project
    files in your local test QAI.
5. Check the history of the HIV and HCV rules files in the `micall/resistance`
    folder. If they have changed, create a new display file in the `docs` folder
    and upgrade the version numbers in the `genreport.yaml` file and
    `asi_algorithm.py`.
5. Check the history of the `micall.alignment` folder. If it has changed since
    the last release, then update the version number in `setup.py`.
5. Update the change notes in the Singularity file, and commit those changes.
6. [Create a release][release] on Github. Use "vX.Y" as the tag, where X.Y
    matches the version you used in QAI. If you have to redo
    a release, you can create additional releases with tags vX.Y.1, vX.Y.2, and
    so on. Mark the release as pre-release until you finish deploying it.
7. Rebuild the Singularity image, and upload it to your local Kive server.
    Process the microtest data.
7. Upload the Singularity image to the Kive test server, and record the
    ids of the new apps.
8. Process all the samples from test_samples.csv on the Kive test server, and
    run the micall_watcher service on a VirtualBox. Use the
    `release_test_*.py` scripts to compare the results of the new release with
    the previous version. Also run the internal scripts `miseq_gen_results.rb`
    and `miseq_compare_results.rb` to look for differences. Get the comparison
    signed off to begin the release process.
8. Upload the Singularity image to the main Kive server, and
    record the id of the new apps.
8. Upload the pipeline definitions to QAI, using the `projects_upload.py`
    script. There is no need to create the new pipeline version in QAI beforehand,
    the script will do this for you - just remember to update the `Order by` field
    afterwards.
8. Stop the micall_watcher service on the main Kive server after you check that
    it's not processing any important runs.

        ssh user@server
        tail /var/log/micall/micall.log
        sudo systemctl stop micall_watcher

9. Get the code from Github into the server's environment.

        ssh user@server
        cd /usr/local/share/MiCall
        git fetch
        git checkout tags/vX.Y

10. Look for changes in [`micall_watcher.py`'s `parse_args()` function][parse_args].
    Either look at the blame annotations at the link above, or review the
    changes in the new release. If there are new or changed settings, adjust
    the configuration in `/etc/systemd/system/micall_watcher.service` or
    `/etc/micall/micall.conf`.
11. Update the container app ids and pipeline version number in
    `/etc/systemd/system/micall_watcher.service`. If you change the configuration, reload it:

        sudo systemctl daemon-reload

12. Check that the kiveapi package is the same version you tested with. If not,
    do a Kive release first.

        cd /usr/local/share/Kive
        /usr/local/share/venv-micall/bin/pip show kiveapi
        cat api/setup.py

13. Start the micall_watcher service, and tail the log to see that it begins
    processing all the runs with the new version of the pipeline.

        sudo systemctl start micall_watcher
        sudo systemctl status micall_watcher
        tail -f /var/log/micall/micall.log

    If the log doesn't help, look in `/var/log/messages` on CentOS or
    `/var/log/syslog` on Ubuntu.

14. Launch the basespace virtual machine, and build a new Docker image
    from GitHub. Tag it with the release number, and push it to the Illumina
    repository. You might have to log in to docker before running this.

        cd /media/sf_micall
        sudo python3 docker_build.py -t vX.Y

16. Edit the `callbacks.js` in the form builder, and add the `:vX.Y` tag to the
    `containerImageId` field.
17. Activate the new revisions in the form builder and the report builder.
17. Submit the new revision of the MiCall app for review, and ask our contact
    at Illumina to set the price to zero for the list of test users.
17. Tag the same docker image and push it to docker hub. Unfortunately, the old
    version of docker that comes with the basespace virtual machine
    [can't log in] to docker hub, so you'll have to save it to a tar file and
    load that into your host system's version of docker.

        ssh basespace@localhost -p2222
        cd /media/sf_micall
        sudo docker tag docker.illumina.com/cfe_lab/micall:vX.Y cfelab/micall:vX.Y
        sudo docker save cfelab/micall:vX.Y >micall-vX.Y.tar
        exit
        sudo docker load <micall-vX.Y.tar
        sudo docker push cfelab/micall:vX.Y
        rm micall-vX.Y.tar

18. Send an e-mail to users describing the major changes in the release.
19. Close the milestone for this release, create one for the next release, and
    decide which issues you will include in that milestone.
20. When the release is stable, check that it's included on the [Zenodo] page.
    If you included more than one tag in the same release, the new tags have
    not triggered Zenodo versions. Edit the release on GitHub, copy the
    description text, update the release, then click the Delete button. Then
    create a new release with the same description, and that will trigger a
    Zenodo version.

[release]: https://help.github.com/categories/85/articles
[parse_args]: https://github.com/cfe-lab/MiCall/blame/master/micall_watcher.py
[Zenodo]: https://doi.org/10.5281/zenodo.2644171
[can't log in]: https://www.docker.com/blog/registry-v1-api-deprecation/
