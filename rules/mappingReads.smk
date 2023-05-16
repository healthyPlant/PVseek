"""
    mapping reads to plant virus database, then identify virus stains according to the aligned reads
    use the sequence with maximum reads in a taxon as a reference for the virus
"""
__author__ = 'Xiaojun Hu <xiaojun.hu@usda.gov>'

workDir = config["workDir"]
logDir = workDir + "/" + config["run_info"]["log"]
mapDir = workDir + "/" + config["run_info"]["map"]
reportDir = workDir + "/" + config["run_info"]["report"]
trimDir = workDir + "/" + config["run_info"]["trim"]
platform = config["platform"] # Illumina_sRNA
mappingTool = config["mappingTool"]
bowtie_param = config["bowtie_param"]
bowtie2_param = config["bowtie2_param"]
minimap2_param = config["minimap2_param"]
bwa_param = config["bwa_param"]
viralRefInfo = config["viralRefInfo"] #/db/viralRef.info.txt
virusTaxon = config["virusTaxon"]
virusDb = config["virusDb"]
readNumThd = config["readNumThd"] #minimu read number for a virus to be selected
genomeCovThd = config["genomeCovThd"] #minimum genome coverage for a virus to be selected


if platform == "Illumina_pe":
  rule pairRead_map:
      """
      mapping paired-end reads to viral references
      """
      input:
        reads1 = trimDir + "/{sample}_R1.trimmed.fastq.gz",
        reads2 = trimDir + "/{sample}_R2.trimmed.fastq.gz",
      output:
        bam = mapDir + "/{sample}.bam",
      message:
        '''--- mapping {wildcards.sample} reads to viral reference database.'''
      threads: config["number_of_threads"]
      log:
        logDir + "/mapping/{sample}.map2Ref.log"
      shell:
        """
        touch {output.bam}
        if [ "{mappingTool}" == "bowtie2" ]; then
          bowtie2_index={virusDb}.1.bt2
          #if the index file does not exist, make one
          if [ ! -s $bowtie2_index ]; then
            bowtie2-build --threads {threads} {virusDb} {virusDb}
          fi
          bowtie2 -p {threads} -q {bowtie2_param} -x {virusDb} -1 {input.reads1} -2 {input.reads2} | samtools view -Sb -F 4 | samtools sort - -o {output.bam} 2>{log}
          samtools index {output.bam}
        else
          bwa_index={virusDb}.bwt
          if [ ! -s $bwa_index ]; then
            bwa index {virusDb}
          fi
          bwa mem {bwa_param} -t {threads} {virusDb} {input.reads1} {input.reads2} | samtools view -Sb -F 4 - | samtools sort -o {output.bam} 2> {log}
          samtools index {output.bam}
        fi
        """
elif platform == "HiPlex" or platform == "Illumina_se":
  rule singleRead_map:
      """
      mapping HiPlex or RNA-seq single reads to viral references
      """
      input:
        reads = trimDir + "/{sample}.trimmed.fastq.gz",
      output:
        bam = mapDir + "/{sample}.bam",
      message:
        '''--- mapping {wildcards.sample} reads to viral reference database.'''
      threads: config["number_of_threads"]
      log:
        logDir + "/mapping/{sample}.map2Ref.log"
      shell:
        """
        touch {output.bam}
        if [ "{mappingTool}" == "bowtie2" ]; then
          bowtie2_index={virusDb}.1.bt2
          #if the index file does not exist, make one
          if [ ! -s $bowtie2_index ]; then
            bowtie2-build --threads {threads} {virusDb} {virusDb}
          fi
          bowtie2 -p {threads} -q {bowtie2_param} -x {virusDb} -U {input.reads} | samtools view -Sb -F 4 | samtools sort - -o {output.bam} 2>{log}
          samtools index {output.bam}
        else
          bwa_index={virusDb}.bwt
          if [ ! -s $bwa_index ]; then
            bwa index {virusDb}
          fi
          bwa mem {bwa_param} -t {threads} {virusDb} {input.reads} | samtools view -Sb -F 4 - | samtools sort -o {output.bam} 2> {log}
          samtools index {output.bam}
        fi
        """
elif platform == "Illumina_sRNA":
  rule sRNA_map:
      """
      mapping small RNA reads to viral references
      """
      input:
        reads = trimDir + "/{sample}.trimmed.fastq.gz",
      output:
        bam = mapDir + "/{sample}.bam",
      message:
        '''--- mapping {wildcards.sample} reads to viral reference database.'''
      threads: config["number_of_threads"]
      log:
        logDir + "/mapping/{sample}.map2Ref.log"
      shell:
        """
        bowtie_index={virusDb}.1.ebwt
        #if the index file does not exist, make one
        if [ ! -s $bowtie_index ]; then
          bowtie-build --threads {threads} {virusDb} {virusDb}
        fi
        bowtie -p {threads} -q {bowtie_param} -x {virusDb} -S {input.reads} | samtools view -Sb -F 4 | samtools sort - -o {output.bam} 2>{log}
        samtools index {output.bam}
        """

elif platform == "Nanopore":
  rule nanopore_map:
      """
      mapping Nanopore reads to viral references
      """
      input:
        reads = trimDir + "/{sample}.trimmed.fastq.gz",
      output:
        bam = mapDir + "/{sample}.bam",
      message:
        '''--- mapping {wildcards.sample} reads to viral reference database.'''
      threads: config["number_of_threads"]
      log:
        logDir + "/mapping/{sample}.map2Ref.log"
      shell:
        """
        if [[ ! -s ${virusDb}.mmi ]]; then
          minimap2 -d {virusDb}.mmi {virusDb}  # indexing the reference database
        fi
        minimap2 {minimap2_param} {virusDb}.mmi {input.reads} | samtools view -Sb -F 4 | samtools sort - -o {output.bam} 2>{log}
        samtools index {output.bam}
        """
else:
  sys.exit("Error: invalid sequencing platform. Must be 'Illumina_pe', 'Illumina_se', 'Illumina_sRNA', 'HiPlex', or 'Nanopore'")


rule identify_virus:
    """
    Identify virus strain based on reads aligned to the virus. The reference is the sequence with maximum reads in the taxon
    """
    input:
      bam = mapDir + "/{sample}.bam",
    output:
      mappedRead = mapDir + "/{sample}.mappedRead.txt",
      coverage = mapDir + "/{sample}.coverage.txt",
      done = touch(logDir + "/checkPoint/{sample}.mapping.done")
    message:
      '''--- Identify viruses in {wildcards.sample}.'''
    params:
      outPrefix = mapDir + "/{sample}.mappedRef"
    shell:
      """
      touch {output.mappedRead}
      touch {output.coverage}
      if [ -s {input.bam} ]; then
        samtools view {input.bam} | cut -f3 | sort | uniq -c | sort -nr > {output.mappedRead}
        #identify the depth at each locus from a bam file.
        samtools depth {input.bam} > {output.coverage}
      fi
      #get reference with highest mapped reads in a taxon
      python {scripts_dir}/getMappedRef.py -m {output.mappedRead} -d {output.coverage} -t {virusTaxon} -r {viralRefInfo} -n {readNumThd} -g {genomeCovThd} -o {params.outPrefix}
      """
