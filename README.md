
This project contains example data for testing MiCall.

To use it, run

```shell
cd data
docker run --rm -v .:/data cfelab/micall folder --project HIVB --skip trim.censor --denovo --keep_scratch . results
```

# Download

Open [this link](https://codeload.github.com/cfe-lab/MiCall/legacy.zip/refs/heads/example-inputs) to download this repository as a ZIP archive.

If you are having issues when downloading the sample inputs, follow [these instructions](https://www.wikihow.com/Download-a-GitHub-Folder).

# Fast Data

Additionally, this repository contains inputs that are much smaller. Use these inputs if you don't want to wait about an hour for the `data/` inputs to be processed.

To use the smaller inputs, run:

```shell
cd fastdata
docker run --rm -v .:/data cfelab/micall folder --project HIVB --skip trim.censor --denovo --keep_scratch . results
```

These commands should take about 1 minute to complete.
