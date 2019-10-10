# gatk3-4-rnaseq-germline-snps-indels
**GATK3 and this workflow is now longer supported, this repo is intended for legacy purposes. Pleas visit [gatk4-rnaseq-germline-snps-indels](https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels) for the latest version of the workflow**  

Workflows for processing RNA data for germline short variant discovery with GATK (v3+v4) and related tools 

### Requirements/expectations :
 - BAM 

### Output :
 - A BAM file and its index.
 - A VCF file and its index. 
 - A Filtered VCF file and its index. 

 Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
 For program versions, see docker containers.

### Important Note :
- The workflow allows users to opt to use GATK4 for all tools instead of the default combination between 
GATK3 and GATK4, but note the workflow has not been validated to run using all GATK4 tools. 
- Runtime parameters are optimized for Broad's Google Cloud Platform implementation. 
- For help running workflows on the Google Cloud Platform or locally please
view the following tutorial [(How to) Execute Workflows from the gatk-workflows Git Organization](https://software.broadinstitute.org/gatk/documentation/article?id=12521)

### LICENSING :
 This script is released under the WDL source code license (BSD-3) (see LICENSE in
 https://github.com/broadinstitute/wdl). Note however that the programs it calls may
 be subject to different licenses. Users are responsible for checking that they are
 authorized to run all programs before running this script. Please see the docker
 page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
 licensing information pertaining to the included programs. 
