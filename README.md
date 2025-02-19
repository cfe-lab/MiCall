
This project contains example data for testing MiCall.

To use it, run

```shell
cd data
micall micall_docker folder --project HIVB --skip trim.censor --denovo --keep_scratch . results
```

Or, with docker version of MiCall:

```shell
cd data
docker run --rm -v .:/data cfelab/micall folder --project HIVB --skip trim.censor --denovo --keep_scratch . results
```

# Download

If you are having issues when downloading the sample inputs, follow [these instructions](https://www.wikihow.com/Download-a-GitHub-Folder).
