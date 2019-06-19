# Genomics Workshop - NextFlow, Docker and Bioconda on Azure


## Pre-requisites
1. Set up Avere vFXT
    - Launch the Avere vFXT using ARM or in the Azure portal
    - Download https://aka.ms/genomicsworkshop/sample_data as genomics_workshop_data.tar
    - The default export mount on Avere is /msazure
    - Make a subdir in the export point: /msazure/data
    - Expand the genomics_workshop_data.tar into /msazure/data so that the data is in it shows as /msazure/data/genomics_workshop_data

    _As an alternate to Avere, use the NFS server `/shared` in the cluster below. But the paths in the nextflow scripts and the instructions have to be modified._

2. A CycleCloud installation

    - Setup CycleCloud in the Avere vFXT subnet
    - Clone/Download the HLS-workshop CycleCloud project (https://github.com/jermth/hls-workshop)
    - Upload the hls-workshop project in CycleCloud
    - Import the cluster template (cyclecloud import_template Genomics-Workshop -f Genomics-Workshop.cluster.txt)

3. Launch the Genomics-Workshop cluster

## Preamble

 1. The cluster template sets up an SGE cluster
 2. You'll need to change the Avere IPAddress in the cluster setup, but the other settings don't have to be set
    - Leave the Avere mount point as /avere -- the example scripts below work on that assumption
 3. This cluster has a couple of cluster-init projects defined that does the following:
    - For installing docker runtime on each node
    - Pulling a couple of biocontainer images for the examples
    - Creating an /etc/profile.d script that stages the example netflow files and symlinking
 4. The exercises illustrate the following:
    - How to install packages using Conda. This is all in the user's home directories and available throughout the cluster.
    - How to run docker images.
    - How to do the the above through nextflow
    - How to use nextflow on the SGE

	
## Exercise 1: Datasets 

 1. Log into the Cluster headnode
 2. On first login, the `~/genomics-workshop/` directory is staged for the user
 3. That directory has a symlink to the human genome index and sample Fastqs
    - bowtie2_index/grch38/
    - NA12787/
 4. Verify the Avere mount there: 

        ls /avere/data/genomics_workshop_data
 6. The `~/genomics-workshop/` directory also carries two example nextflow files.

## Exercise 2: Application install using Conda

1. Install conda. This is done for the user account, so run the following as the user:

    ```
    # Get the installer
    wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh

    # Run the installer
    bash Miniconda2-latest-Linux-x86_64.sh  -b

    # Set up the path and environment
    echo "source ~/miniconda2/etc/profile.d/conda.sh" >> ~/.bashrc
    echo "export PATH=~/miniconda2/bin:$PATH" >> ~/.bashrc
    source ~/miniconda2/etc/profile.d/conda.sh
    export PATH=~/miniconda2/bin:$PATH
    ```

2. Install tools:

    ```
    # Set up the conda "channels"
    conda config --add channels bioconda
    conda config --add channels conda-forge

    # install bioinformatics tools
    conda install -y bowtie2 samtools bcftools htop glances nextflow
    ```

3. Run a local bowtie job with the conda tools

    ```
    cd ~/genomics-workshop
    mkdir exercise2

    bowtie2 -x bowtie2_index/grch38/grch38_1kgmaj -1 NA12787/SRR622461_1.fastq.gz -2 NA12787/SRR622461_2.fastq.gz  -u 100000 -S exercise2/exercise2.sam

    # view the file
    head -n 50 exercise2/exercise2.sam

    Convert the sam file into bam and sort it
    samtools view -bS exercise2/exercise2.sam | samtools sort > exercise2/exercise2.bam

    Generate a VCF from the bam
    bcftools mpileup --no-reference exercise2/exercise2.bam > exercise2/exercise2.vcf

    # view the VCF
    head -n 50 exercise2/exercise2.vcf
    ```

## Exercise 3: Docker
1. Run the same set of bowtie and samtools commands, but using Docker containers

    ```
    # view the docker images on this VM
    docker image list

    # Run bowtie2 again, but now using the container instead of the conda install binary

    cd ~/genomics-workshop
    mkdir docker_results
    docker run -u $(id -u ${USER}):$(id -g ${USER})  -v ~/genomics-workshop/docker_results:/results -v /avere/data/genomics_workshop_data:/data biocontainers/bowtie2:v2.3.1_cv1  bowtie2 -x bowtie2_index/grch38/grch38_1kgmaj -1 NA12787/SRR622461_1.fastq.gz -2 NA12787/SRR622461_2.fastq.gz  -u 100000 -S /results/exercise3.sam

    # if you have another window, run docker ps to see the container running

    docker run -u $(id -u ${USER}):$(id -g ${USER})  -v ~/genomics-workshop/docker_results:/results -v /avere/data/genomics_workshop_data:/data biocontainers/samtools:v1.7.0_cv4 /bin/bash -c "samtools view -bS /results/exercise3.sam | samtools sort -o /results/exercise3.bam"
    ```

## Exercise 4: Nextflow with conda tools
1. Open and look at the bowtie.local.nf file. This nf file describes a workflow. The processes defined run the same set of commands in exercise2. Now however, the Nextflow workflow manager chains together the outputs.

2. Run the same set of applications locally using nextflow:
	
        nextflow run bowtie2.local.nf

## Exercise 5: Nextflow with docker containers
- Nextflow works with docker conainers as well. Look at the differences between bowtie2.container.nf and bowtie2.local.nf

        nextflow run bowtie2.container.nf

## Exercise 6: Nextflow with docker, on SGE
- Finally, use Nextflow with the SGE scheduler by specifying the execution engine in the config file. 
- Using the same NF flow file, the docker tasks are now submitted to the SGE compute nodes. 

        nextflow -c nextflow.sge.config run bowtie2.container.nf
- Verify that the job is in queue by using the `qstat -f` command
- Notice also that the autoscaler kicks in and scales up a compute node.
