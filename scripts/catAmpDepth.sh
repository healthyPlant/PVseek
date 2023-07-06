#!/usr/bin/env bash
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 05/08/2023
#This program concatenates all virus amplicon depths from different samples 
#Usage: bash path/to/PVseek/scripts/catAmpDepth.sh [primer name file] [amplicon graph folder] [project name] [PVseek work dir] 
######################################################################

primerFile="$1" #db/Pomes_HiPlex_virus_primerName.txt
ampFolder="$2" #~/test_PVseek/MonsterplexRun99/amplicon #which has *.ampDepth.txt
project="$3" #MonsterplexRun99
workdir="$4" #~/test_PVseek/MonsterplexRun99

echo $# arguments 
if [ "$#" -le 2 ]; then
    echo "Usage: bash /path/to/PVseek/scripts/catAmpDepth.sh [primer name file] [amplicon graph folder] [project name] [PVseek work dir] "
    echo "Example: bash /path/to/PVseek/scripts/catAmpDepth.sh /path/to/PVseek/db/Pomes_HiPlex_virus_primerName.txt /path/to/test_PVseek/MonsterplexRun99/amplicon MonsterplexRun99 /path/to/test_PVseek/MonsterplexRun99"
    exit
fi

while read -r line; do
    name=$(echo "$line" | tr -d '\n') 
    echo "Virus: $name"
    #concatenate all files
    cat $ampFolder/*.$name.*ampDepth.txt > $workdir/$project.$name.ampDepth.txt
    #Check whether the command is valid or invalid
    if [ $? -ne 0 ]; then
        echo "Virus: $name is not in $project."
    fi

    #only edit not empty files
    if [ -s $workdir/$project.$name.ampDepth.txt ]; then
        #just keep the first line header and delete others
        sed -i '1!{/^Sample/d;}'  $workdir/$project.$name.ampDepth.txt
    else
        rm -rf $workdir/$project.$name.ampDepth.txt
    fi
done < "$primerFile"