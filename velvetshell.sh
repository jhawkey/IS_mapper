#!/bin/bash
#running VelvetOptimiser through the IS mapping pipeline, assembly done in temp folder and final contigs
#file renamed and copied to output folder
#velvet folder deleted at end of pipeline run

#if there are six arguments given
if [ $# -eq 6 ]; then
	#run VelvetOptimiser
	cd "$1" && VelvetOptimiser.pl -s "$2" -e "$3" -f " -short -bam $4 " -d "$5"

	#if this runs correctly then copy the contigs file to its new location and name
	if [ $? -eq 0 ]; then
		cp "$5"/contigs.fa "$6"

	#otherwise report exit status
	else
		echo "VelvetOptimiser failed with exit status $?"
		exit $?
	fi

#if six arguments aren't given, show the usage and report back no zero exit status
else
	echo "usage: velvetshell.sh velvetOutputDir sValue eValue bamFile finalVelvetDir newContigFile"
	exit -1
fi

#all ran correctly, report 0 exit status
exit 0
