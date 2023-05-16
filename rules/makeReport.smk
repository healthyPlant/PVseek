"""
    Make a final report 
"""
__author__ = 'Xiaojun Hu <xiaojun.hu@usda.gov>'

workDir = config["workDir"]
mapDir = workDir + "/" + config["run_info"]["map"]
reportDir = workDir + "/" + config["run_info"]["report"]
logDir = workDir + "/" + config["run_info"]["log"]
badSeqIds = config["badSeqIds"]
reportFile = reportDir + "/report.txt"
samples = config["samples"]
coverage = config["coverage"]
mappingTool = config["mappingTool"]
platform = config["platform"]
virusDb = config["virusDb"]
virusTaxon = config["virusTaxon"]
viralRefInfo = config["viralRefInfo"]

rule generate_report:
    """
    Generate the final report
    """
    input:
	    mappingDone = expand(logDir + "/checkPoint/{sample}.mapping.done", sample=samples)
    output:
	    reportFile, #reportDir + "/report.txt"
    message:
	    '''--- Generate the final report.'''
    priority: -100
    shell:
        """
        python {scripts_dir}/summarizeMapping.py -w {workDir} -f {badSeqIds} -o {output}
        #draw coverage graph
        if [[ {coverage} == "yes" ]]; then
            if [[ {platform} != "HiPlex" ]]; then
                python {scripts_dir}/batchPlotViralCoverage.py -w {workDir} -d {virusDb} -s {scripts_dir} -m {mappingTool} -o {workDir}/mapping/coverage
            fi
        fi
        """
