workflow RNAseq {

	File inputBam
	String sampleName = basename(inputBam,".bam")

	File refFasta
	File refFastaIndex
	File refDict
	
	File wgsCallingIntervalList

	Array[File] knownVcfs
	Array[File] knownVcfsIndices

    File dbSnpVcf
    File dbSnpVcfIndex

	## Inputs for STAR
	Int? readLength
	File? zippedStarReferences
	File annotationsGTF
  
  	## Optional user optimizations
  	Int? haplotypeScatterCount
    Int scatterCount = select_first([haplotypeScatterCount, 6])

    Int? preemptible_tries
    Int preemptible_count = select_first([preemptible_tries, 3])

	call RevertSam {
		input:
			input_bam = inputBam,
			base_name = sampleName + ".reverted",
			sort_order = "queryname",
			preemptible_count = preemptible_count
	}

	call SamToFastq {
		input:
		unmapped_bam = RevertSam.output_bam,
		base_name = sampleName,
		preemptible_count = preemptible_count
	}

	if (!defined(zippedStarReferences)) {

		call StarGenerateReferences { 
			input:
			ref_fasta = refFasta,
			ref_fasta_index = refFastaIndex,
			annotations_gtf = annotationsGTF,
			read_length = readLength,
			preemptible_count = preemptible_count
		}
	}

	File starReferences = select_first([zippedStarReferences,StarGenerateReferences.star_genome_refs_zipped,""])

	call StarAlign { 
		input: 
			star_genome_refs_zipped = starReferences,
			fastq1 = SamToFastq.fastq1,
			fastq2 = SamToFastq.fastq2,
			base_name = sampleName + ".star",
			read_length = readLength,
			preemptible_count = preemptible_count
	}

	call MergeBamAlignment {
		input: 
			unaligned_bam = RevertSam.output_bam,
			star_bam = StarAlign.output_bam,
			base_name = ".merged",
			ref_fasta = refFasta,
			ref_dict = refDict,
			preemptible_count = preemptible_count
	}

	call MarkDuplicates {
		input:
			input_bam = MergeBamAlignment.output_bam,
			base_name = sampleName + ".dedupped",
			preemptible_count = preemptible_count
	}

	call SplitNCigarReads {
		input:
			input_bam = MarkDuplicates.output_bam,
			input_bam_index = MarkDuplicates.output_bam_index,
			base_name = sampleName + ".split",
			ref_fasta = refFasta,
			ref_fasta_index = refFastaIndex,
			ref_dict = refDict,
			interval_list = wgsCallingIntervalList,
			preemptible_count = preemptible_count
	}

	call BaseRecalibrator {
		input:
			input_bam = SplitNCigarReads.output_bam,
			input_bam_index = SplitNCigarReads.output_bam_index,
			recal_output_file = sampleName + ".recal_data.csv",
  			dbSNP_vcf = dbSnpVcf,
  			dbSNP_vcf_index = dbSnpVcfIndex,
  			known_indels_sites_VCFs = knownVcfs,
  			known_indels_sites_indices = knownVcfsIndices,
  			ref_dict = refDict,
  			ref_fasta = refFasta,
  			ref_fasta_index = refFastaIndex,
  			preemptible_count = preemptible_count
	}

	call ApplyBQSR {
		input:
			input_bam = SplitNCigarReads.output_bam,
			input_bam_index = SplitNCigarReads.output_bam_index,
			base_name = sampleName + ".aligned.duplicates_marked.recalibrated",
			ref_fasta = refFasta,
			ref_fasta_index = refFastaIndex,
			ref_dict = refDict,
			recalibration_report = BaseRecalibrator.recalibration_report,
			preemptible_count = preemptible_count
	}

	call ScatterIntervalList {
		input:
			interval_list = wgsCallingIntervalList,
      		scatter_count = scatterCount,
      		preemptible_count = preemptible_count
	}

    scatter (i in range(ScatterIntervalList.interval_count)) {
		call HaplotypeCaller {
			input:
				input_bam = ApplyBQSR.output_bam,
				input_bam_index = ApplyBQSR.output_bam_index,
				base_name = sampleName + ".hc",
				interval_list = ScatterIntervalList.out[i],
				ref_fasta = refFasta,
				ref_fasta_index = refFastaIndex,
				ref_dict = refDict,
				preemptible_count = preemptible_count
		}
	}

	call MergeVCFs {
    	input:
      		input_vcfs = HaplotypeCaller.output_gvcf,
      		input_vcfs_indexes = HaplotypeCaller.output_gvcf_index,
      		output_vcf_name = sampleName + ".g.vcf.gz",
      		preemptible_count = preemptible_count
  	}

	call VariantFiltration {
		input:
			input_vcf = MergeVCFs.output_vcf,
			input_vcf_index = MergeVCFs.output_vcf_index,
			base_name = sampleName + ".variant_filtered.vcf.gz",
			ref_fasta = refFasta,
			ref_fasta_index = refFastaIndex,
			ref_dict = refDict,
			preemptible_count = preemptible_count
	}

	output {
		File recalibrated_bam = ApplyBQSR.output_bam
		File recalibrated_bam_index = ApplyBQSR.output_bam_index
		File merged_gvcf = MergeVCFs.output_vcf
		File merged_gvcf_index = MergeVCFs.output_vcf_index
		File variant_filtered_gvcf = VariantFiltration.output_vcf
		File variant_filtered_gvcf_index = VariantFiltration.output_vcf_index
	}
}


#NOTE: assuming aggregated bams & paired end fastqs
task SamToFastq {
	File unmapped_bam
	String base_name

	Int preemptible_count

	command <<<
	 	java -jar /usr/gitc/picard.jar \
	 	    SamToFastq \
	 	    INPUT=${unmapped_bam} \
	 	    VALIDATION_STRINGENCY=SILENT \
	 	    FASTQ=${base_name}.1.fastq.gz \
	 	    SECOND_END_FASTQ=${base_name}.2.fastq.gz
	>>>

	output {
		File fastq1 = "${base_name}.1.fastq.gz"
    	File fastq2 = "${base_name}.2.fastq.gz"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
		memory: "4 GB"
		disks: "local-disk " + sub(((size(unmapped_bam,"GB")+1)*5),"\\..*","") + " HDD"
		preemptible: preemptible_count
	}
}

task StarGenerateReferences {
    File ref_fasta
    File ref_fasta_index
    File annotations_gtf
    Int? read_length  ## Should this be an input, or should this always be determined by reading the first line of a fastq input

    Int? num_threads
    Int threads = select_first([num_threads, 8])
    Int preemptible_count

    command <<<
        mkdir STAR2_5

        STAR \
            --runMode genomeGenerate \
            --genomeDir STAR2_5 \
            --genomeFastaFiles ${ref_fasta} \
            --sjdbGTFfile ${annotations_gtf} \
            ${"--sjdbOverhang "+(read_length-1)} \
            --runThreadN ${threads}

        ls STAR2_5

        tar -zcvf star-HUMAN-refs.tar.gz STAR2_5
    >>>

    output {
        Array[File] star_logs = glob("*.out")
        File star_genome_refs_zipped = "star-HUMAN-refs.tar.gz"
    }

    runtime {
        docker: "ruchim/star:v1"
        disks: "local-disk 80 HDD"
        cpu: threads
        memory: "25 GB"
        preemptible: preemptible_count
    }
}


task StarAlign {
    File star_genome_refs_zipped
    File fastq1
    File fastq2
    String base_name
    Int? read_length

    Int? num_threads
    Int threads = select_first([num_threads, 8])
    Int? star_mem_max_gb
    Int star_mem = select_first([star_mem_max_gb, 45])
    #Is there an appropriate default for this?
    Int? star_limitOutSJcollapsed
    Int preemptible_count

    command <<<
        set -e

        tar xvzf ${star_genome_refs_zipped}

        STAR \
            --genomeDir STAR2_5 \
            --runThreadN ${threads} \
            --readFilesIn ${fastq1} ${fastq2} \
            --readFilesCommand "gunzip -c" \
            ${"--sjdbOverhang "+(read_length-1)} \
            --outSAMtype BAM SortedByCoordinate \
            --twopassMode Basic \
            --limitBAMsortRAM ${star_mem+"000000000"} \
            --limitOutSJcollapsed ${default=1000000 star_limitOutSJcollapsed} \
            --outFileNamePrefix ${base_name}.
    >>>

    output {
        File output_bam = "${base_name}.Aligned.sortedByCoord.out.bam"
        File output_log_final = "${base_name}.Log.final.out"
        File output_log = "${base_name}.Log.out"
        File output_log_progress = "${base_name}.Log.progress.out"
        File output_SJ = "${base_name}.SJ.out.tab"
    }

    runtime {
        docker: "ruchim/star:v1"
        disks: "local-disk " + sub(((size(fastq1,"GB")+size(fastq2,"GB")*10)+30),"\\..*","") + " HDD"
        memory: (star_mem+1) + " GB"
        cpu: threads
        preemptible: preemptible_count
    }
}

task MergeBamAlignment {

    File ref_fasta
    File ref_dict

    File unaligned_bam
    File star_bam
    String base_name

    Int preemptible_count
    #Using default for max_records_in_ram
 
    command <<<
        java -jar /usr/gitc/picard.jar \
            MergeBamAlignment \
            REFERENCE_SEQUENCE=${ref_fasta} \
            UNMAPPED_BAM=${unaligned_bam} \
            ALIGNED_BAM=${star_bam} \
            OUTPUT=${base_name}.bam \
            INCLUDE_SECONDARY_ALIGNMENTS=false \
            PAIRED_RUN=False \
            VALIDATION_STRINGENCY=SILENT
    >>>
 
    output {
        File output_bam="${base_name}.bam"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
        disks: "local-disk " + sub(((size(unaligned_bam,"GB")+size(star_bam,"GB")+1)*5),"\\..*","") + " HDD"
        memory: "4 GB"
        preemptible: preemptible_count
    }
}

task MarkDuplicates {

 	File input_bam
 	String base_name

 	Int preemptible_count

 	command <<<
 	    java -jar /usr/gitc/picard.jar \
 	        MarkDuplicates \
 	        INPUT=${input_bam} \
 	        OUTPUT=${base_name}.bam  \
 	        CREATE_INDEX=true \
 	        VALIDATION_STRINGENCY=SILENT \
 	        METRICS_FILE=${base_name}.metrics
 	>>>

 	output {
 		File output_bam = "${base_name}.bam"
 		File output_bam_index = "${base_name}.bai"
 		File metrics_file = "${base_name}.metrics"
 	}

	runtime {
		disks: "local-disk " + sub(((size(input_bam,"GB")+1)*3),"\\..*","") + " HDD"
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
		memory: "4 GB"
		preemptible: preemptible_count
	}
}

## Not validated in GATK4 
task SplitNCigarReads {

	File input_bam
	File input_bam_index
	String base_name
	File interval_list

	File ref_fasta
	File ref_fasta_index
	File ref_dict

	Int preemptible_count

    command <<<
    	java -jar /usr/gitc/GATK35.jar \
    		-T SplitNCigarReads \
    		-R ${ref_fasta} \
    		-I ${input_bam} \
    		-o ${base_name}.bam \
    		-rf ReassignOneMappingQuality \
    		-RMQF 255 \
    		-RMQT 60 \
    		-U ALLOW_N_CIGAR_READS
    >>>

 	output {
 		File output_bam = "${base_name}.bam"
 		File output_bam_index = "${base_name}.bai"
 	}

    runtime {
    	disks: "local-disk " + sub(((size(input_bam,"GB")+1)*5 + size(ref_fasta,"GB")),"\\..*","") + " HDD"
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
		memory: "4 GB"
    	preemptible: preemptible_count
    }
}


task BaseRecalibrator {

    File input_bam
    File input_bam_index
    String recal_output_file

    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    Int preemptible_count

    command <<<
        /gatk/gatk-launch --javaOptions "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R ${ref_fasta} \
            -I ${input_bam} \
            --useOriginalQualities \
            -O ${recal_output_file} \
            -knownSites ${dbSNP_vcf} \
            -knownSites ${sep=" -knownSites " known_indels_sites_VCFs}
    >>>

    output {
        File recalibration_report = recal_output_file
    }

    runtime {
        memory: "6 GB"
        disks: "local-disk " + sub((size(input_bam,"GB")*3)+30, "\\..*", "") + " HDD"
        docker: "broadinstitute/gatk:4.beta.6"
        preemptible: preemptible_count
    }
}


task ApplyBQSR {

    File input_bam
    File input_bam_index
    String base_name
    File recalibration_report

    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int preemptible_count

    command <<<
        /gatk/gatk-launch \
            --javaOptions "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
            -XX:+PrintGCDetails -Xloggc:gc_log.log \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
            ApplyBQSR \
            --addOutputSAMProgramRecord \
            -R ${ref_fasta} \
            -I ${input_bam} \
            --useOriginalQualities \
            -O ${base_name}.bam \
            -bqsr ${recalibration_report}
    >>>

    output {
        File output_bam = "${base_name}.bam"
        File output_bam_index = "${base_name}.bai"
    }

    runtime {
        memory: "3500 MB"
        disks: "local-disk " + sub((size(input_bam,"GB")*4)+30, "\\..*", "") + " HDD"
        preemptible: preemptible_count
        docker: "broadinstitute/gatk:4.beta.6"
    }
}

task HaplotypeCaller {

	File input_bam
	File input_bam_index
	String base_name

	File interval_list

  	File ref_dict
  	File ref_fasta
  	File ref_fasta_index

  	Int preemptible_count

	command <<<
		java -jar /usr/gitc/GATK35.jar \
		    -T HaplotypeCaller \
		    -R ${ref_fasta} \
		    -I ${input_bam} \
		    -L ${interval_list} \
		    -dontUseSoftClippedBases \
		    -stand_call_conf 20.0 \
		    -o ${base_name}.vcf.gz
	>>>

    output {
        File output_gvcf = "${base_name}.vcf.gz"
        File output_gvcf_index = "${base_name}.vcf.gz.tbi"
    }

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
		memory: "3 GB"
		disks: "local-disk " + sub((size(input_bam,"GB")*2)+30, "\\..*", "") + " HDD"
		preemptible: preemptible_count
	}
}

task VariantFiltration {

	File input_vcf
	File input_vcf_index
	String base_name

  	File ref_dict
  	File ref_fasta
  	File ref_fasta_index

  	Int preemptible_count

	command <<<
		/gatk/gatk-launch \
		    VariantFiltration \
			-R ${ref_fasta} \
			-V ${input_vcf} \
			-window 35 \
			-cluster 3 \
			-filterName "FS" \
			-filter "FS > 30.0" \
			-filterName "QD" \
			-filter "QD < 2.0" \
			-O ${base_name}
	>>>

	output {
    	File output_vcf = "${base_name}"
    	File output_vcf_index = "${base_name}.tbi"
	}

	runtime {
		docker: "broadinstitute/gatk:4.beta.6"
		memory: "3 GB"
		disks: "local-disk " + sub((size(input_vcf,"GB")*2)+30, "\\..*", "") + " HDD"
		preemptible: preemptible_count
	}
}

task MergeVCFs {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name

    Int? disk_size = 5
    Int preemptible_count

    # Using MergeVcfs instead of GatherVcfs so we can create indices
    # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
    command <<<
        java -Xms2000m -jar /usr/gitc/picard.jar \
            MergeVcfs \
            INPUT=${sep=' INPUT=' input_vcfs} \
            OUTPUT=${output_vcf_name}
    >>>

    output {
        File output_vcf = output_vcf_name
        File output_vcf_index = "${output_vcf_name}.tbi"
    }

    runtime {
        memory: "3 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
        preemptible: preemptible_count
    }
}

task ScatterIntervalList {

    File interval_list
    Int scatter_count

    Int preemptible_count

    command <<<
        set -e
        mkdir out
        java -Xms1g -jar /usr/gitc/picard.jar \
            IntervalListTools \
            SCATTER_COUNT=${scatter_count} \
            SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
            UNIQUE=true \
            SORT=true \
            INPUT=${interval_list} \
            OUTPUT=out

        python3 <<CODE
        import glob, os
        # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
        intervals = sorted(glob.glob("out/*/*.interval_list"))
        for i, interval in enumerate(intervals):
          (directory, filename) = os.path.split(interval)
          newName = os.path.join(directory, str(i + 1) + filename)
          os.rename(interval, newName)
        print(len(intervals))
        CODE
    >>>

    output {
        Array[File] out = glob("out/*/*.interval_list")
        Int interval_count = read_int(stdout())
    }

    runtime {
        disks: "local-disk 1 HDD"
        memory: "2 GB"
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
        preemptible: preemptible_count
    }
}

task RevertSam {
    File input_bam
    String base_name
    String sort_order

    Int preemptible_count

    command <<<
        java -jar /usr/gitc/picard.jar \
        	RevertSam \
        	INPUT=${input_bam} \
        	OUTPUT=${base_name}.bam \
        	VALIDATION_STRINGENCY=SILENT \
        	ATTRIBUTE_TO_CLEAR=FT \
        	ATTRIBUTE_TO_CLEAR=CO \
        	SORT_ORDER=${sort_order}
    >>>

    output {
        File output_bam = "${base_name}.bam"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
        disks: "local-disk " + sub(((size(input_bam,"GB")+1)*5),"\\..*","") + " HDD"
        memory: "4 GB"
        preemptible: preemptible_count
    }
}
