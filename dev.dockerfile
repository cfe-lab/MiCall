# Dockerized version of MiCall development tools.
# To build this image, first build the production image in Dockerfile, probably
# using the docker_build.py script. Use the docker images command to see that
# docker.illumina.com/cfe_lab/micall:latest is the image you want to use, then
# build this image, using a command like this:
#
#     docker build -t micall:dev --file dev.dockerfile .
#
# To test out the image, run the test suite, with a command like this:
#
#     docker run --rm -it --entrypoint pytest -w /opt/micall \
#     --volume ~/git/micall:/opt/micall micall:dev
#
# That lets you edit the source code on your host system, but run it under
# docker with all the tools installed for you.

FROM docker.illumina.com/cfe_lab/micall:latest

## Add the dev packages.
COPY pyproject.toml /opt/micall/
RUN pip install -r .[dev]
