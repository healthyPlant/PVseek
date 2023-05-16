"""
    Raw reads are cleaned by fastp
"""
__author__ = 'Xiaojun Hu <xiaojun.hu@usda.gov>'

workDir = config["workDir"]
rawReadDir = workDir + "/" + config["run_info"]["raw"]
logDir = workDir + "/" + config["run_info"]["log"]
trimDir = workDir + "/" + config["run_info"]["trim"]
reportDir = workDir + "/" + config["run_info"]["report"]
input_format = config["input_format"]
strand1 = config["strand1"]
strand2 = config["strand2"]
seq_type = config["seq_type"]

samples = config["samples"]

fastp_params = ""
if platform == "Illumina_pe":
  fastp_params = config["fastp_pe_param"]
elif platform == "Illumina_se":
  fastp_params = config["fastp_se_param"]
elif platform == "Illumina_sRNA":
  fastp_params = config["fastp_srna_param"]
elif platform == "HiPlex":
  fastp_params = config["fastp_hiplex_param"]
elif platform == "Nanopore":
  fastp_params = config["fastp_nanopore_param"]
else:
  sys.exit("Error: invalid sequencing platform. Must be 'Illumina_pe', 'Illumina_se', 'Illumina_sRNA', 'HiPlex', or 'Nanopore'")
  
rawFiles = glob.glob(rawReadDir + "/*." + input_format)
true001 = 0
for rawFile in rawFiles:
  if "_001." in rawFile:
    true001 = 1
    break

# Define input files
def raw_inputs(wildcards):
    if true001 and strand1:  
      return rawReadDir + "/{sample}_" + strand1 + "_001." + input_format
    elif strand1:
      return rawReadDir + "/{sample}_" + strand1 + "." + input_format
    else:
      return rawReadDir + "/{sample}." + input_format

def raw_pe_inputs(wildcards):
    if true001 and strand2:  
      return expand(rawReadDir + "/{sample}_{strand}_001." + input_format, strand=[strand1,strand2], sample=wildcards.sample)
    elif strand2:
      return expand(rawReadDir + "/{sample}_{strand}." + input_format,  strand=[strand1,strand2], sample=wildcards.sample)
    else:
      sys.exit("Error: invalid sequencing type parameter. ")


if seq_type == 'se':
  rule run_fastp_se:
      """
      Trim adapters and low quality reads by fastp
      """
      input:
        reads = raw_inputs #rawReadDir + "/{sample}_" + strand + "_001." + input_format, 
      output:
        reads = trimDir + "/{sample}.trimmed.fastq.gz",
        json = trimDir + "/{sample}.fastp.json",
        html = trimDir + "/{sample}.fastp.html",
      message:
        '''--- Trim {wildcards.sample} reads using fastp.'''
      log:
        logDir + "/fastp/{sample}.log"
      params: 
        fastp_params,
      threads: 8 #config["number_of_threads"] #config.get("number_of_threads", 1)
      shell:
        """     
        fastp --thread {threads} {fastp_params} -i {input.reads} -o {output.reads} --json {output.json} --html {output.html} > {log} 
        """
elif seq_type == 'pe':
  rule run_fastp_pe:
      """
      Trim adapters and low quality reads by fastp
      """
      input:
        reads = raw_pe_inputs #rawReadDir + "/{sample}_" + strand + "_001." + input_format,
      output:
        reads1 = trimDir + "/{sample}_R1.trimmed.fastq.gz",
        reads2 = trimDir + "/{sample}_R2.trimmed.fastq.gz",
        json = trimDir + "/{sample}.fastp.json",
        html = trimDir + "/{sample}.fastp.html",
      message:
        '''--- Trim {wildcards.sample} reads using fastp.'''
      log:
        logDir + "/fastp/{sample}.log"
      params: 
        fastp_params,
      threads: 8 #config["number_of_threads"] #config.get("number_of_threads", 1)
      shell:
        """     
        fastp --thread {threads} {fastp_params} -i {input.reads[0]} -I {input.reads[1]} -o {output.reads1} -O {output.reads2} --json {output.json} --html {output.html} > {log} 
        """

rule sum_readStats:
  """
  Summary raw read stats 
  """
  input:
    expand(trimDir + "/{sample}.fastp.json", sample=samples),
  output:
    reportDir + "/qcReadNumber.txt"
  message:
    '''--- Summary raw and trimmed read stats.'''
  shell:
    """
    python {scripts_dir}/sumFastp.py -w {trimDir} -o {output}
    """
