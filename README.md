# gatk3-4-rnaseq-germline-snps-indels
 Workflows for processing RNA data for germline short variant discovery with GATK (v3+v4) and related tools 

### Requirements/expectations :
 - BAM 

### Output :
 - A BAM file and its index.
 - A VCF file and its index. 
 - A Filtered VCF file and its index. 

 Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
 For program versions, see docker containers.

 LICENSING :
 This script is released under the WDL source code license (BSD-3) (see LICENSE in
 https://github.com/broadinstitute/wdl). Note however that the programs it calls may
 be subject to different licenses. Users are responsible for checking that they are
 authorized to run all programs before running this script. Please see the docker
 page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
 licensing information pertaining to the included programs. 
