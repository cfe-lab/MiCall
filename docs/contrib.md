---
title: Contributing to the MiCall Project
---

If you like this project and want to make it better, help out. You could report
a bug, or pitch in with some development work.

## Bug Reports and Enhancement Requests

Please create issue descriptions [on GitHub][issues]. Be as specific as possible.
Which version are you using? What did you do? What did you expect to happen? Are
you planning to submit your own fix in a pull request?

[issues]: https://github.com/cfe-lab/MiCall/issues

## Development

The easiest way to start developing MiCall is by using DevContainers.

1. **Open Project**:
    - If you're using Visual Studio Code on your local machine, open the MiCall project folder and select the "Reopen in Container" prompt to initialize the DevContainer environment. Make sure you have the necessary DevContainer extension installed beforehand, available [here](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers).
    - For a web-based development environment, you can develop directly on GitHub using GitHub Codespaces by navigating to the MiCall repository on GitHub and selecting "Code" > "Open with Codespaces" to launch a pre-configured environment.

2. **Dependency Installation**: All required dependencies will be automatically installed whether you are using a local DevContainer or GitHub Codespace.
3. **Verification**: To ensure that the environment is correctly configured, execute `pytest` within the DevContainer or Codespace. All tests should pass, indicating that the setup is successful.

### Local install

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

```shell
python --version
```

We have tested with Python `3.12`.

### Test data
If you want to run `micall watcher`, you have to set up data folders for raw
data and for the working folders. You'll also need to set up the QAI project
and the MiseqQCReport so you can download QC data and upload results.

1. Create a data folder somewhere on your workstation, like ~/data. Create
   subdirectories called miseq and RAW_DATA. Add folders RAW_DATA/MiSeq/runs.
2. Connect to the shared drive [using CIFS][cifs] and mount
   `smb://192.168.68.144/RAW_DATA` as `/media/RAW_DATA`.
3. Navigate down to /media/RAW_DATA/MiSeq/runs, pick a recent folder, and make
   sure it has a file named needsprocessing.
4. Copy SampleSheet.csv to a sample run folder under your local
   RAW_DATA/MiSeq/runs folder.
5. Navigate down to Data\Intensities\BaseCalls, and copy a few of the .fastq.gz
   files to your sample run folder under Data/Intensities/BaseCalls.
6. Copy the Interop folder and the files RunInfo.xml and runParameters.xml.
7. Open settings.py for editing.
8. Point `home` at your local data/miseq folder.
9. Point `rawdata_mount` at your local RAW_DATA folder.
10. Set the Oracle connection information to a test database where you can upload
   sequence data. (https://git-int.cfenet.ubc.ca/wscott/oracleserver)
11. Run the Ruby console for QAI and `LabMiseqRun.import('01-Jan-2000')` for the
    date of your sample run.
12. Upload the projects to a micall pipelines in QAI, use `micall.utils.projects_upload` to create a new pipeline in QAI
13. Run `micall watcher`, it does need arguments. Look up the container app ids from Kive, check the Kive server URL and ports as well as QAI server and port

[cifs]: https://wiki.ubuntu.com/MountWindowsSharesPermanently

### Looking at SAM files
When you don't understand the pipeline's output, it can be helpful to look at
the raw reads in a sequence viewer like [Tablet][tablet]. Run the `micall analyze`
command on a run folder or a single sample, like this:

```shell
micall analyze folder --debug_remap --all_projects --keep_scratch /path/to/run
```

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

```shell
# Enter MiCall's docs directory.
cd docs/
# Build and run the documentation website.
docker build --tag micall-documentation .
docker run \
       --volume .:/docs \
       --publish 4000:80 \
       --rm \
       --name micall-documentation \
       --interactive --tty \
       micall-documentation
# Goto http://localhost:4000
```

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

### Updating HIVdb

When a new version of the Stanford HIV Drug Resistance Database is released,
follow these steps to update MiCall:

1. **Download the new rules file**: Visit [https://hivdb.stanford.edu/](https://hivdb.stanford.edu/)
   and download the latest XML algorithm file. The file can typically be found
   in the algorithm/downloads section or by searching for "ASI XML" or "algorithm XML".

2. **Note the version details**: Record the version number (e.g., 9.8) and the
   modification date from the XML file or the HIVdb website.

3. **Add the new XML file**: Place the downloaded XML file in `micall/resistance/`
   with the naming convention `HIVDB_X.X.xml` (e.g., `HIVDB_9.8.xml`).

4. **Update `HIV_RULES_PATH`**: In `micall/resistance/resistance.py`, update
   the path to point to the new XML file:
   ```python
   HIVDB_VERSION = 'X.X'
   HIV_RULES_PATH = os.path.join(os.path.dirname(__file__), f'HIVDB_{HIVDB_VERSION}.xml')
   ```

5. **Update version in genreport**: In `micall/resistance/genreport.yaml`,
   update the version number and modification date in the `generated_by_text`
   field (look for "Stanford HIVdb"):
   ```yaml
   generated_by_text: >
     Generated by MiCall {} on Illumina BaseSpace using Stanford HIVdb
     X.X, modified on YYYY-MM-DD.
   ```

6. **Handle new drugs** (if applicable): If the new HIVdb version introduces
   new drugs, add them to the appropriate drug class in `genreport.yaml`.
   For example, to add a new INSTI drug:
   ```yaml
   INSTI:
     - [BIC, Bictegravir]
     - [NEW, New Drug Name]  # Add new drug here
     - [DTG, Dolutegravir]
   ```
   Use the format: `[ABBREVIATION, Full Drug Name]`. The order determines how
   drugs appear in the report.

7. **Check pyvdrm compatibility**: Verify the current pyvdrm version in
   `pyproject.toml` supports the new HIVdb version. The pyvdrm library parses
   the HIVdb XML format. If compatibility issues arise:
   - Check the [pyvdrm repository](https://github.com/cfe-lab/pyvdrm) for updates
   - Update the version in `pyproject.toml` if needed:
     ```toml
     "pyvdrm @ git+https://github.com/cfe-lab/pyvdrm.git@vX.X.X"
     ```

8. **Update tests**: Modify test files to reflect the new version and any
   algorithm behavior changes:
   - `micall/tests/test_resistance.py`: Update version number assertions
   - `micall/tests/test_asi_algorithm.py`: Update expected resistance scores
     if the algorithm produces different results

9. **Run tests**: Execute the test suite to ensure everything works:
    ```shell
    pytest micall/tests/test_resistance.py
    pytest micall/tests/test_asi_algorithm.py
    pytest  # Run all tests
    ```

10.  **Manual verification**: Process a sample dataset through the pipeline and
    review the generated resistance reports to ensure output is correct and
    formatting is preserved.

**Note on CA (Capsid) region**: Some HIVdb versions include the CA region for
Capsid Inhibitors. If present and not yet supported in MiCall, the region is
automatically skipped by the `get_algorithm_regions()` function in
`resistance.py`. Decide whether to add CA support based on clinical needs.

### Releases
This section assumes you already have a working server up and running, and you
just want to publish a new release. If you're setting up a new server, follow
similar steps to setting up a development workstation. Follow these steps:

1. Check that all the issues in the current milestone are closed, and make sure
    the code works in your development environment. Run all the unit
    tests as described above, process the microtest data set with
    `micall release_test_microtest`.
2. Check if the kiveapi package needs a new release by looking for new commits.
    Make sure you tested with the latest version.
3. Determine what version number should be used next.
4. Use the `micall projects_dump` script for the previous version and compare
    `projects.json` to check that the projects match, or that the differences
    were intended. Test the `micall projects_upload` script with your updated project
    files in your local test QAI.
5. Check the history of the HIV and HCV rules files in the `micall/resistance`
    folder. If they have changed, create a new display file in the `docs` folder
    and upgrade the version numbers in the `genreport.yaml` file and
    `asi_algorithm.py`.
5. Check the history of the `micall.alignment` folder. If it has changed since
    the last release, then update the version number in `setup.py`.
5. Update the change notes in the Singularity file, and commit those changes.
6. [Create a release][release] on Github. Use "vX.Y.Z" as the tag, where X.Y
    matches the version you used in QAI. If you have to redo
    a release, you can create additional releases with tags vX.Y.Z.1, vX.Y.Z.2, and
    so on. Mark the release as pre-release until you finish deploying it.
7. Rebuild the Singularity image, and upload it to your local Kive server.
    Process the microtest data.
7. Upload the Singularity image to the Kive test server, and record the
    ids of the new apps.
8. Process all the samples from test_samples.csv on the Kive test server, and
    run the `micall_watcher` service on a VirtualBox. Use the
    `micall release_test_*` scripts to compare the results of the new release with
    the previous version. Also run the internal scripts `miseq_gen_results.rb`
    and `miseq_compare_results.rb` to look for differences. Get the comparison
    signed off to begin the release process.
8. Upload the Singularity image to the main Kive server, and
    record the id of the new apps.
9. Upload the pipeline definitions to QAI, using the `micall projects_upload`
    command. There is no need to create the new pipeline version in QAI beforehand,
    the command will do this for you - just remember to update the `Order by` field
    afterwards.
10. Stop the watcher service on the main Kive server after you check that
    it's not processing any important runs.

    ```shell
    ssh user@server
    tail /var/log/micall/micall.log
    sudo systemctl stop micall_watcher
    ```

11. Get the code from Github into the server's environment.

    ```shell
    ssh user@server
    cd /usr/local/share/MiCall
    git fetch
    git checkout tags/vX.Y.Z
    ```

12. Look for changes in [`micall/monitor/watcher.py`'s `parse_args()` function][parse_args].
    Either look at the blame annotations at the link above, or review the
    changes in the new release. If there are new or changed settings, adjust
    the configuration in `/etc/systemd/system/micall_watcher.service` or
    `/etc/micall/micall.conf`.
13. Update the container app ids and pipeline version number in
    `/etc/systemd/system/micall_watcher.service`. If you change the configuration, reload it:

    ```shell
    sudo systemctl daemon-reload
    ```

14. Check that the kiveapi package is the same version you tested with. If not,
    do a Kive release first.

    ```shell
    cd /usr/local/share/Kive
    /usr/local/share/venv-micall/bin/pip show kiveapi
    cat api/setup.py
    ```

15. Start the watcher service, and tail the log to see that it begins
    processing all the runs with the new version of the pipeline.

    ```shell
    sudo systemctl start micall_watcher
    sudo systemctl status micall_watcher
    tail -f /var/log/micall/micall.log
    ```

    If the log doesn't help, look in `/var/log/messages` on CentOS or
    `/var/log/syslog` on Ubuntu.

16. Upload a BaseSpace app. See [basespace_upload.md](./basespace_upload.md) for instructions on how to do this.

17. Send an e-mail to users describing the major changes in the release.
18. Close the milestone for this release, create one for the next release, and
    decide which issues you will include in that milestone.
19. When the release is stable, check that it's included on the [Zenodo] page.
    If you included more than one tag in the same release, the new tags have
    not triggered Zenodo versions. Edit the release on GitHub, copy the
    description text, update the release, then click the Delete button. Then
    create a new release with the same description, and that will trigger a
    Zenodo version.

[release]: https://help.github.com/categories/85/articles
[parse_args]: https://github.com/cfe-lab/MiCall/blame/master/micall/monitor/watcher.py
[Zenodo]: https://doi.org/10.5281/zenodo.2644171
[can't log in]: https://www.docker.com/blog/registry-v1-api-deprecation/
[docker hub]: https://hub.docker.com/orgs/cfelab/members
