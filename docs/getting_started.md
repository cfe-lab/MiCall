---
title: Getting Started
subtitle: Running MiCall on its own
---

If you want to process a few FASTQ files or MiSeq run folders through MiCall,
this page will walk you through installing Docker and using it to run MiCall.
If you want to set up a more automated process, then it might be worth taking
the extra time to set up the tools described on the [admin] page.

## Install Docker
Docker is a tool that makes it easier for us to share MiCall with you. Instead
of listing the many tools that MiCall is built on and telling you how to install
them, you only have to install one: Docker. You can think of Docker as running
a smaller computer inside your computer, and MiCall comes with a script that
tells Docker exactly how to set up that smaller computer with all the tools that
MiCall needs.

The first step is to follow the instructions for installing [Docker]. You should
be able to install it on just about any kind of computer.

## Build the Docker Image
Once you've installed Docker and successfully run the `hello-world` container,
it's time to set up MiCall. First, you have to get the source code. You can
either clone the GitHub repository, or download the source archive from the
latest [release].

Here's what the commands might look like for cloning the GitHub repository and
checking out version `vX.Y`.

    git clone https://github.com/cfe-lab/MiCall.git
    cd MiCall
    git checkout tags/vX.Y

Here's what the commands might look like for downloading the source archive and
unpacking it.

    wget https://github.com/cfe-lab/MiCall/archive/vX.Y.zip
    unzip vX.Y.zip
    rm vX.Y.zip
    cd MiCall-X.Y

Once you have the source code, you can build the docker image with a command
like this:

    sudo docker build --tag micall:X.Y .

If you're not running Docker on Linux, the command may look slightly different.
See the [Docker build] tutorial for more details.

## Run Docker
To check that the container built successfully, try this command:

    sudo docker run micall:X.Y --help

That should display a help message that lists all the commands MiCall
understands. You probably want to start with either the `folder` or the `sample`
command, and read the help messages for them. These two commands will display
those help messages:

    sudo docker run micall:X.Y folder --help
    sudo docker run micall:X.Y sample --help

One slight annoyance with running MiCall under Docker is that all the result
files are owned by the root user. If you don't like that, use the `chown`
command to change their owner.

[admin]: admin.md
[Docker]: https://docs.docker.com/get-started/
[release]: https://github.com/cfe-lab/MiCall/releases
[Docker build]: https://docs.docker.com/get-started/part2/
