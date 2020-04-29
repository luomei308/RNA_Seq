# Run RNA-seq data

A Nextflow implementation of the Tuxedo Suite of Tools Workflow is based on the 2016 Nature Protocols publication: ["Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown"](http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html).

> nextflow>=0.22.0
test: nextflow: version 18.10.1 build 5003
          java: java version "1.8.0_181"

Please learn about the basic knowledge for docker (https://www.runoob.com/docker/docker-tutorial.html)

## More details about parameters
Check the [Tuxedo-NF README.md](https://github.com/evanfloden/tuxedo-nf/blob/master/README.md)

## Quick start
Make sure you have all the required dependencies listed in the last section.

Install the Nextflow runtime by running the following command:

``` shell
$ curl -fsSL get.nextflow.io | bash
$ docker pull dockerluom/tuxedo-nf
```

you should do in the following: 
  1. cope all files into your dir
  2. configure the parameters in "nextflow.config"
  3. include your fq files, genome gtf, index files.
  4. change the max tasks by change the maxForks in xxx.nf file

When done, you can launch the pipeline execution by entering the command shown below:

``` shell
$ nextflow tuxedo.nf -with-docker dockerluom/tuxedo-nf -c nextflow.config (non-trim)
$ nextflow tuxedo_trim.nf -with-docker dockerluom/tuxedo-nf -c nextflow.config (including trim steps)
```


