"""
    Snakemake main file 
"""
__author__ = 'Xiaojun Hu <xiaojun.hu@usda.gov>'

import os,subprocess,itertools,glob,re,sys
import datetime
#import snakemake.workflow as workflow
from snakemake.utils import report
import time

start = time.time()

#SETUP PROCESS-RELATED CONFIGURATION 
try:
	CONFIGFILE = str(workflow.overwrite_configfiles[0])
except:
	CONFIGFILE = str(workflow.overwrite_configfile[0])

#config["workDir"] = os.getcwd()
#setup workdir dynamically
#rules will then be expressed relative to the workdir:
workdir: config["workDir"]  #can be changed at the snakemake command line
fastqDir = config["fastqDir"]  #can be changed at the snakemake command line

#Get relative position of Snakefile from wd
SNAKEFILE = workflow.snakefile
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)
CORES = workflow.cores #note workflow.cores work only on version 5.11.2 above
#print(CORES)
config["number_of_threads"] = CORES  #change threads number by command paramter --cores/-j 16
config["snakefile"] = SNAKEFILE
config["snakefile_dir"] = SNAKEFILE_DIR
#print(SNAKEFILE)

#print(workflow.__dict__) #print all attributes
#Establish snakefile and environment dictionaries
snakefiles_dir = os.path.abspath(os.path.join(SNAKEFILE_DIR, "rules"))		#directory for additional snakefiles
environments_dir = os.path.abspath(os.path.join(SNAKEFILE_DIR, "envs"))	#directory for conda environment yaml files
scripts_dir = os.path.abspath(os.path.join(SNAKEFILE_DIR, "scripts")) 			#directory for extra scripts used in workflow

#running folders
rawReadDir = config["run_info"]["raw"]
logDir = config["run_info"]["log"]
mapDir = config["run_info"]["map"]
reportDir = config["run_info"]["report"]
reportFile = reportDir + "/report.txt" 
seq_type = config["seq_type"]

#if os.path.isfile(reportFile):
#    os.remove(reportFile)

## Check whether input fastq files are compressed or not
if  "input_format" in config.keys():
	input_format = config["input_format"]
else:
	input_format = "fastq.gz"

#setup soft link or copy fastq files to rawReadDir
if fastqDir:
	#set up soft link for fastq files if fastq folder is surpplied
	if not os.path.exists(rawReadDir):
		os.system("ln -sf " + fastqDir + " " + rawReadDir)
	# Check if given path is link
	if os.path.exists(rawReadDir) and os.path.islink(rawReadDir):
		print(fastqDir, " soft link is built." )
	else:
		# soft link broken, copy files
		if not os.path.isdir(rawReadDir): #check a folder exists
			os.mkdir(rawReadDir) #make a folder
		os.system("cp " + fastqDir + "/*." + input_format + " " + rawReadDir + "/")  
			
else:
	if not os.listdir(rawReadDir): 
		print(rawReadDir + " is empty. Exit!") 
		sys.exit()

#get sample names from rawReadDir
SAMPLES = [os.path.basename(x) for x in glob.glob(rawReadDir + '/*.' +  input_format)]
#remove 'Undetermined_*'
SAMPLES = [x for x in SAMPLES if not x.startswith('Undetermined_')]
#print(SAMPLES)
#remove '_R1_001.fastq.gz'
SAMPLES = [x.replace("." + input_format,"") for x in SAMPLES]
if config["strand1"] != '':
	SAMPLES = [re.sub(r"(.*)_%s.*" % config["strand1"], "\\1", x) for x in SAMPLES]

if seq_type == 'pe':
	SAMPLES = [re.sub(r"(.*)_%s.*" % config["strand2"], "\\1", x) for x in SAMPLES]

SAMPLES = list(set(SAMPLES))  #remove duplicated from a list using set

print("Samples: %s"  % ", ".join(SAMPLES))
config["samples"] = SAMPLES  #save sample names in a global viriable
#print("config samples: ", config['samples'])

if len(SAMPLES) == 0:
	sys.exit("No sample files found. Exit!")

#Run rules
include: os.path.join(snakefiles_dir, "mappingReads.smk")
include: os.path.join(snakefiles_dir, "cleanReads.smk")
include: os.path.join(snakefiles_dir, "makeReport.smk")

rule all:
	input: #targets
		reportFile,
		reportDir + "/qcReadNumber.txt",
	message: "Rule all"
	shell: "echo Job done    `date '+%Y-%m-%d %H:%M'`"	


# The onsuccess handler is executed if the workflow finished without error.
onsuccess:
	#check filtered results, if no viruses, report them
	noViralSamples = []
	for sample in SAMPLES:
		filename = mapDir + "/" + sample + ".mappedRef.filtered.txt"
		if os.stat(filename).st_size == 0:
			noViralSamples.append(sample)
	failNum = len(noViralSamples)

	print("\n############################################################################")
	print("# PV-seek finished without errors ")
	if failNum > 0:
		print("# ", failNum, " samples without viruses found: ", ", ".join(noViralSamples))
	print("# Plase check the result in the folder: " + reportDir)
	print("############################################################################")
	print("\nRunning time in minutes: %s\n" % round((time.time() - start)/60,1))

# Else, the onerror handler is executed.
onerror:
	print("\n\n#####################\n# An error occurred #\n#####################\n\n")
	print("\nRunning time in minutes: %s\n" % round((time.time() - start)/60,1))

# onstart handler will be executed before the workflow starts. Note that dry-runs do not trigger any of the handlers
onstart:
	print("Running PV-seek for samples : %s"  % ", ".join(SAMPLES))