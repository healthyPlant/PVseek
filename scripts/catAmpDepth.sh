#!/usr/bin/env bash
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 05/08/2023
#This program concatenates all virus amplicon depths from different samples 
#Usage: bash catAmpDepth.sh [path/to/script name] [primer name file] [project name] [PVseek work dir] 
######################################################################

primerFile="$1" #db/Pomes_HiPlex_virus_primerName.txt
project="$2" #MonsterplexRun99
workdir="$3" #/ppq/data0/test_PVseek/MonsterplexRun99

echo $# arguments 
if [ "$#" -le 2 ]; then
    echo "Usage: bash catAmpDepth.sh [path/to/script name] [primer name file] [project name] [PVseek work dir] "
    echo "Example: bash catAmpDepth.sh /my/PVseek/scripts/catAmpDepth.sh /my/PVseek/db/Pomes_HiPlex_virus_primerName.txt MonsterplexRun99 /my/test_PVseek/MonsterplexRun99"
    exit
fi

while read -r line; do
    name=$(echo "$line" | tr -d '\n') 
    echo "Virus: $name"
    #concatenate all files
    cat $workdir/amplicon/*.$name.*ampDepth.txt > $workdir/$project.$name.ampDepth.txt
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