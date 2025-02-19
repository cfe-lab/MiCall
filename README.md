
This project contains example data for testing MiCall.

To use it, run

```shell
cd data
micall micall_docker folder --project HIVB --skip trim.censor --denovo --keep_scratch . results
```

Or, with docker version of MiCall:

```shell
cd data
docker run --rm -v .:/data micall folder --project HIVB --skip trim.censor --denovo --keep_scratch . results
```
