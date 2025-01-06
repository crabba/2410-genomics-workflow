# 2410-genomics-workflow
Nextflow workflow for fastqc to mutect2

## Introduction
This Nextflow workflow is a step by step migration of the FASTQ to Mutect2 workflow initially developed in Galaxy. Each step is implemented as a Docker container and runs as an independent process with its own local file system.  The Nextflow engine coordinates the running and resubmission of tasks, task order dependency, and the piping of files between tasks.

Pre-existing public Docker containers have been used where they are available, and when not available, Dockerfiles are provided. Not all input parameters have been ported from the Galaxy workflow into Nextflow; some modification and validation of each task and the overall workflow will be required to generate results identical to those created by Galaxy.  

The workflow may be tested on a local instance (virtual machine), with the Nextflow engine, Docker, and each Docker container all running on the same machine.  Small modifications have been added to allow the same workflow to run on AWS, with input and output data stored on Amazon S3, and all compute infrastructure automatically created upon demand for the duration of the run. The AWS-specific modifications in no way alter the function of the workflow and should result in identical outputs.

Usage
```
nextflow run main.nf --fastq 'H06HDADXX130110.1.*_{1,2}.fastq.gz' -params-file params.json
```

## Parameters file

Nextflow accepts a JSON parameters file at the command line with the `-params-file` flag. While Nextflow can read a nested JSON file, AWS HealthOmics requires a flat (single level) JSON parameters file.  This repository contains a flat JSON file where each key has a unique name, including the step number.

In the sample `params.json` file, the following fields differ when running on HealthOmics rather than locally using Nextflow

**Input files**
When running locally using Nextflow and Docker, input files may come from either the local file system or S3.  Note that when sourced from S3, Nextflow will download the files from S3 and re-run all steps, even when the `-resume` flag is used.  This is because the file is downloaded to ensure it is identical with the previous run, and in doing so will gain a new modification date, which is one of the conditions for Nextflow caching.

```
# Input files come from S3 (optional for local execution, required for HealthOmics)
    "fastq_1": "s3://MY-BUCKET/MY-PATH/R1.fq.gz",
    "ref_prefix": "s3://MY-BUCKET/ref/hg38/",    
# Input files come from local file system (optional for local execution)
    "fastq_1": "R1.fq.gz",
    "ref_prefix": "ref/hg38",
```

**Container images**
When running locally using Nextflow and Docker, input files may come from either the local Docker deployment, or from ECR.
```
# Container images come from AWS ECR (optional for local execution, required for HealthOmics)
    "omics_ecr_prefix": "426997564500.dkr.ecr.us-east-1.amazonaws.com/",
# Input files come from local Docker (optional for local execution)
    "omics_ecr_prefix": ""
```

**Omics switch**
This switch controls whether the Nextflow engine writes results to the special directory name `/mnt/workflow/pubdir`
```
    "onOmics": true,
```

## Test Environment

This workflow may be tested on any Linux machine with the following software installed:
* Docker
* Nextflow
* Containers as listed in the parameters file

## Migration from Galaxy to Nextflow

This section describes the creation of a Nextflow workflow that replicates the original Galaxy workflow.  This workflow can be run and debugged locally and does not include HealthOmics-specific attribues.

### Identify Step Tools

Examine the Galaxy workflow diagram and its accompanying `.ga` file. Tasks have been assigned numbers in the `.ga` file and it is convenient to retain these numbers, and the application they run, throughout the migration to Nextflow. Example: As Fastq To Bam is step 4 in Galaxy, this step in Nextflow is called `fastq_to_bam_04` since process and variable names may not start with a digit.

Within the step, the relevant tool is found in the `tool_id` and `tool_version` lines:

```
"tool_id": "kcctools.ohsu.edu/repos/onwuzu/fg_fastq_to_bam/fg_FastqToBam/2.2.1+galaxy1",
"tool_version": "2.2.1+galaxy1",
```

### Locate Matching Tool Images

Docker images matching the tool and version number may be found in public repositories such as Docker Hub and BioContainers.  Alternatively a Docker image can be built from a Dockerfile (see `fgbio` example in repository). 

### Extract Tool Parameters

In each step, the line `tool_state` contains the inputs used in this step.  In the `.ga` file, this file is non-standard JSON:

```
"tool_state": "{\"__job_resource\": {\"__job_resource__select\": \"no\", 
...
\"umi_tag\": \"\", \"__page__\": null, \"__rerun_remap_job_id__\": null}",
```

To edit this into a standard JSON document:
* Delete the quote and comma marks surrounding the top level brackets at the beginning and end of the line:
    * `"tool_state": "{` -> `"tool_state": {`
    * `null}",` -> `null}`
* Search and replace any occurrences of `\"\"` to `""` (otherwise Visual Studio will incorrectly delete extra characters in the following replacement step)
* Select and replace all occurrences of `\"` to `"`

The document should now look similar to the following:

```
"tool_state": {"__job_resource": {"__job_resource__select": "no", 
...
"umi_tag": "", "__page__": null, "__rerun_remap_job_id__": null}

```

And can be reformatted:

```
"tool_state": {
    "__job_resource": {
        "__job_resource__select": "no",
...
    "umi_tag": "",
    "__page__": null,
    "__rerun_remap_job_id__": null
}
```

Applicable parameters may be found at the top level or in fields such as `optional` and `inputs`.

### Identify Matching Tool Parameters

Generate a list of parameters from the replacement tool by running the tool container.

```
$ docker run -it broadinstitute/gatk:latest /bin/bash
# java -jar /gatk/gatk.jar
USAGE:  <program name> [-h]

Available Programs:
-------------------------------------------------
Base Calling:                                    
    CheckIlluminaDirectory (Picard)              
    CollectIlluminaBasecallingMetrics (Picard)   
...    
```

Build a command string for each process to replicate the parameters used by Galaxy.  Parameter values can be hard coded or interpolated from definitions in the parameter file. Note that the syntax for Nextflow variables vs shell variables depends upon whether a `shell` vs `script` block is used, and single vs double quotes.

```
process bwa_mem2_08 {
    input:
    path reference_sequence
    
    shell:
    '''
    bwa-mem2 \
    mem \
    -t !{params.steps_08_bwa_mem2_threads} \
    !{reference_sequence} \
    ...
    '''
}
```

### Running Nextflow Locally

Install Nextflow and run from the command line.  Nextflow will launch each process as a Docker container and will run the provided script within the container environment. File system mounting is handled automatically, and each process runs within its own directory, so there is no requirement to use distinct file names between processes.  Result files are written to the working directory (default name `work`). 

Symbolic links are used to avoid copying files between directories.  Most applications will treat a symbolic link the same as a file, but some applications will follow a symbolic link back to the location of the original file.  This can be an issue when files are implicitly assumed to exist in the same directory, such as looking for a `.fasta.fai` file in the same directory as a `.fasta` file.  For this reason, all files implicitly assumed to co-exist must be defined as Nextflow channels and passed as a parameter to a process, even if that file is not explicitly called by the shell script.  This particularly applies to files associated with indexes (`.fai`, `.ann`, `.amb` etc).

**Running Workflows**

Run Nextflow workflows from the command line, passing the parameter file:

```
nextflow run main.nf -resume -with-report -with-trace -params-file params-local-flat.json
```

The report file (HTML) and trace file (text) contain details on each step of the run including per-process CPU and memory usage, and execution time.  This information can be used to tune the requested CPU and memory for each step.

## Running Nextflow on AWS HealthOmics

**Storage**
Nextflow can natively read from and write to S3. A file may be defined locally as `foo_r1.fastq.gz` or on S3 as `s3://mybucket/mypath/foo_r1.fastq.gz` with no change in behaviour. 

Intermediate files are written on a shared high performance file system in a directory specific to the process, and will not be persisted when the run is completed. To save required files to the output S3 location, define the `publishDir` directive with the value `/mnt/workflow/pubdir`.  See [Writing workflows in Nextflow](https://docs.aws.amazon.com/omics/latest/dev/workflow-languages.html#workflow-languages-nextflow) in the HealthOmics documentation. 

**Compute**

All compute nodes are created upon demand when a job is run. An instance (virtual machine) to host each Docker container is created for the duration of a process.  The instance is sized to accommodate the requested `cpus` and `memory` requested for each process.  The instance vCPU count will be a power of 2 (2, 4, 8...) and the memory will be 2, 4, or 8 GB per vCPU.  Both `cpus` and `memory` requested will be satisfied by the instance.

The same container image may be used for local testing, and for production runs in HealthOmics. HealthOmics requires that the container be stored in AWAS ECR (Elastic Container Registry), and the container ARN be defined in the `container` directive of each process.

See: [Compute and memory requirements for HealthOmics tasks](https://docs.aws.amazon.com/omics/latest/dev/memory-and-compute-tasks.html) and [HealthOmics workflow definition requirements](https://docs.aws.amazon.com/omics/latest/dev/workflow-defn-requirements.html)

**Parameters**
Nextflow accepts a single JSON parameter file of arbitrary depth, allowing nested statements to be referenced in the Nextflow script.  The following is a valid parameter file when not running in HealthOmics:

`config.json`
```
{
    steps: {
        foo: {
            "mykey": "myvalue"
        }
    }
}
```

This parameter value is read in the Nextflow source file:

`main.nf`
```
myvar = params.steps.foo.mykey
```

AWS HealthOmics requires that the Nextflow parameter file be of depth 1, as each key in the file results in a HealthOmics parameter which may also be supplied as a console input when running from the AWS console.

`config.json`
```
{
    "steps_foo_mykey": "myvalue"
}
```

`main.nf`
```
myvar = params.steps_foo_mykey
```

**Creating a Workflow**

Follow the process described in [Creating private workflows in HealthOmics](https://docs.aws.amazon.com/omics/latest/dev/workflows-setup.html).  AWS CLI commands to create a HealthOmics private workflow are shown below.

```
# Define private variables
export NEXTFLOW_WORKFLOW_URI=s3://MY-BUCKET/MY-PATH
export PROJECT_NAME=MY-PROJECT-NAME

# Define variables for use in the commands
export TIME_NOW=$(date +%y%m%d-%H%M%S)
export CREATE_WORKFLOW_INPUT_JSON=aws-omics-create-workflow-input-${TIME_NOW}.json
export CREATE_WORKFLOW_OUTPUT_JSON=aws-omics-create-workflow-output-${TIME_NOW}.json

# Zip the Nextflow file and upload to S3
rm main.zip
zip main.zip main.nf
aws s3 cp main.zip ${NEXTFLOW_WORKFLOW_URI}/

# Create the HealthOmics workflow from a populated workflow template file
export NEXTFLOW_PARAM_DESCRIPTION=$(cat params-healthomics-description.json)
# Populates: NEXTFLOW_PARAM_DESCRIPTION, NEXTFLOW_WORKFLOW_URI
envsubst < aws-omics-create-workflow-template.json  > $CREATE_WORKFLOW_INPUT_JSON
aws omics create-workflow --cli-input-json file://${CREATE_WORKFLOW_INPUT_JSON} > $CREATE_WORKFLOW_OUTPUT_JSON
export OMICS_WORKFLOW_ID=$(cat $CREATE_WORKFLOW_OUTPUT_JSON | jq -r .id)
```

Start the HealthOmics run from a populated template file, which includes the ID of the workflow created above, and the contents of the parameters definition file.

```
# Define private variables
export OMICS_CACHE_URI=s3://MY-BUCKET/MY-PATH
export OMICS_OUTPUT_URI=s3://MY-BUCKET/MY-PATH
export OMICS_ROLE_ARN=arn:aws:iam::000000000000:role/MY-ROLE-NAME

# Create a cache (one-time operation)
# Populates: OMICS_CACHE_URI
envsubst < aws-omics-create-run-cache-template.json > aws-omics-create-run-cache-input.json
aws omics create-run-cache --cli-input-json file://aws-omics-create-run-cache-input.json > aws-omics-create-run-cache-output.json
export OMICS_CACHE_ID=$(cat aws-omics-create-run-cache-output.json | jq -r .id)

# Submit run
export TIME_NOW=$(date +%y%m%d-%H%M%S)
export START_RUN_INPUT_JSON=aws-omics-start-run-input-${TIME_NOW}.json
export START_RUN_OUTPUT_JSON=aws-omics-start-run-output-${TIME_NOW}.json
export OMICS_PARAMS_INPUTS=$(cat params-healthomics.json)
# Populates: OMICS_CACHE_ID, TIME_NOW, OMICS_OUTPUT_URI, OMICS_PARAMS_INPUTS, OMICS_ROLE_ARN, OMICS_WORKFLOW_ID
envsubst < aws-omics-start-run-template.json >  $START_RUN_INPUT_JSON
aws omics start-run --cli-input-json file://${START_RUN_INPUT_JSON} > $START_RUN_OUTPUT_JSON
```

**Workflow analytics**

Detailed metrics are recorded from each run, including CPU and memory utilization.  After the run is complete, the Omics Run Analyzer from [AWS HealthOmics Tools](https://github.com/awslabs/amazon-omics-tools) will generate a CSV report showing the allocated resources, percent utilization of each, and estimted cost of each task.

**Security**
HealthOmics, as with all AWS managed services, runs in an AWS-owned account and requires permissions to access resources in any other AWS account. In this example, HealthOmics requires permissions to:

* Read from the input S3 bucket
* Write to the output S3 bucket
* Write to CloudWatch logs
* Retrieve container images from Elastic Container Registry

These permissions are granted through a set of policies attached to a role, which is supplied as an input when a job is submitted to HealthOmics. The role explicitly states trust in the HealthOmics service, which can then assume the role and exercise the privileges granted in the attached policies.

A suitable role is described in the [AWS HealthOmics End to End Workshop](https://catalog.workshops.aws/amazon-omics-end-to-end/en-US), this role may be adapted and modified to suit.
