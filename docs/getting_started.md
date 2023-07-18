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

## Pull the Docker Image
Once you've installed Docker and successfully run the `hello-world` container,
it's time to install MiCall. Choose a version to install by looking at the
[Docker Hub tags]. The tags that include "dev" have not been validated yet, so
you probably want the latest tag without "dev". Once you've chosen a version,
you can pull and run it with a command like this, where `vX.Y` is the tag you've
chosen:

    sudo docker run cfelab/micall:vX.Y --help

That will take a while to download, and then it will display a help message that
lists all the commands MiCall
understands. You probably want to start with either the `folder` or the `sample`
command, and read the help messages for them. These two commands will display
those help messages:

    sudo docker run cfelab/micall:vX.Y folder --help
    sudo docker run cfelab/micall:vX.Y sample --help

One slight annoyance with running MiCall under Docker is that all the result
files are owned by the root user. If you don't like that, use the `chown`
command to change their owner.

## Docker Tips
The help messages will be easier to read if you pass through the `$COLUMNS`
environment variable.

    sudo docker run -e COLUMNS=$COLUMNS cfelab/micall:vX.Y folder --help

Docker will clean up after itself if you use the `--rm` option.

    sudo docker run --rm -e COLUMNS=$COLUMNS cfelab/micall:vX.Y folder --help

If you forget to clean up, you can see a list of all the left-over containers,
then either remove a container by name, or use the prune command to remove all
the exited containers.

    $ sudo docker ps -a
    CONTAINER ID        IMAGE                      COMMAND                  CREATED             STATUS                      PORTS               NAMES
    135421f7c094        cfelab/micall:v7.13.dev0   "python /opt/micall/…"   35 seconds ago      Exited (0) 31 seconds ago                       gallant_gauss
    7567aae23a09        cfelab/micall:v7.13.dev0   "python /opt/micall/…"   45 seconds ago      Exited (0) 40 seconds ago                       strange_greider
    $ sudo docker rm gallant_gauss
    gallant_gauss
    $ sudo docker container prune -f
    Deleted Containers:
    7567aae23a094615f251720e8161a7d0958844f89b357030e05f15ed3cc572eb
    
    Total reclaimed space: 984.4kB
    $

[admin]: admin.md
[Docker]: https://docs.docker.com/get-started/
[Docker Hub tags]: https://hub.docker.com/r/cfelab/micall/tags
