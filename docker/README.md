# Running neonMicrobe R package in Docker Containers

**Why run Docker?**

No need to install anything other than the [Docker Run Engine](https://docs.docker.com/get-docker/), the package and its dependencies will be installed inside the Docker container image. 

**Rocker Project:**

This docker container's base image is the [Rocker Project's](https://www.rocker-project.org/) R base image tag `rocker/r-base:4.0.4`. The `neonMicrobe` Dockerfile in this directory details the dependencies to install `neonMicrobe` on a [debian](https://www.debian.org/) system.

## Command Line Interface

**Running neonMicrobe with Docker Locally**

*Download Docker Container image to local system* 

`docker pull rbartelme/neonmicrobe:v1`

*Example: Docker Run Statement*

`docker run --rm -it -v $(pwd):/work rbartelme/neonmicrobe:v1 bash`

*`docker run` argument meaning in above statement:* 

* `--rm` remove container from local machine after scripts finish. 
    * [docker container cleanup best practices](https://dzone.com/articles/docker-clean-after-yourself)
* `-it` run docker interactively.
* `-v $(pwd):/work` copy contents of current working directory into the container image in the directory `/work`

### Running on HPC with Singularity

---

## Running on CyVerse Discovery Environment 2.0



### RStudio Application

### CyVerse Data Store

If you would like to store your data on CyVerse for running analyses using the CyVerse Data Store, please look at the guide [here.](https://learning.cyverse.org/projects/data_store_guide/en/latest/)