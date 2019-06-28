GENOME= 'hg38'

rule fastqc

rule multi-fastqc

rule trimmomatic

rule trim-muti-fastqc

rule bwa:
	input:
		"R1={sample}.R1.fastq.gz"
		"R2={sample}.R2.fastq.gz"
	output:
		"bwa_out/{sample}.bam"
	shell:
		"bwa mem -M -t 32 GENOME input.R1 input.R2 | samtools view -Sb - > {output}"
#seems name sorted is better than coord sorted. the piard tool can handle both name or coord sorted bam, but with different behavior
rule samtool_sort:
	input:
		"bwa_out/{sample}.bam"
	output:
		"sorted/{sample}.bam"
	shell:
		"samtools sort -T sorted/{wildcards.sample} "
		"-O bam {input} > {output}"

rule remove_dup
	intput:
		"sorted/{sample}.bam"
	output:
		bam="dedup/{sample}.dedup.bam"
		matrix="dedup/{sample}_rmdup_metrics"
	shell:
		"${gatk4} MarkDuplicates "
		"-I {input} -ASO coordinate --CREATE_INDEX true "
		"-O {output} --REMOVE_DUPLICATES true "
		"-M output.matrix --TMP_DIR /scratch "

rule recalibrate
	input:
		"dedup/{sample}.dedup.bam"
	output:
		"{sample}_bsqr.table"
	shell:
		"${gatk4} --java-options "-Xmx70G" BaseRecalibrator "
		"-I {input} -R ${genome_fasta} --known-sites ${dbsnp} --known-sites ${MILLS_INDELS} --known-sites ${kG_snps} "
		"-O {output} --tmp-dir $tempdir
rule first
	input:
		"{sample}_bsqr.table"
	output:
		"{sample}_recall.bam"
	shell:
		"${gatk4} --java-options "-Xmx70G" ApplyBQSR "
		"-R ${genome_fasta} -I dedup/{sample}.dedup.bam --bqsr-recal-file {input} "
		"-O Sample_${f}_recal.bam --tmp-dir ${tmpDir}"
rule hyplotype_variant
	input:
		"{sample}_recall.bam"
	output:
		"{sample}_HC.g.vcf.gz"
	shell:
		"${gatk4} --java-options "-Xmx70g" HaplotypeCaller "
		"-R ${genome_fasta} -I {input} -O {output} --emit-ref-confidence GVCF "
		"--dbsnp ${dbsnp} --tmp-dir ${tmpDir}"

rule consolidate
	intput:
		"{sample}_HC.g.vcf.gz"
	output:

	shell:
		"${gatk4} --java-options "-Xmx100g" GenomicsDBImport "
		"--genomicsdb-workspace-path Strassmann_V4_HC_calls_database_Apri2019_${k} --intervals ${k} --tmp-dir ${tmpDir} ${fns} 2>&1|tee >./genomicsDBimport_${k}.log &

rule joint_call
	input:


	shell:
		${gatk4} --java-options "-Xmx100g" GenotypeGVCFs -R ${genome_fasta} -V gendb://Strassmann_V4_HC_calls_database_Apri2019_${k} -G StandardAnnotation --use-new-qual-calculator -O Strassmann_V4_HC_rawcalls_Apri2019_${k}.vcf --tmp-dir ${tmpDir}








