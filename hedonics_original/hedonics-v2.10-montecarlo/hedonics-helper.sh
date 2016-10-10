#!/bin/bash

#
# This script will provide functions to help with the hedonics program
#


#
# Make a bunch of directories to contain each run. Takes one argument the number of directories to create
#
make-seed-dirs () {
	$num_of_dirs=$1
	for dir_number in `seq -w 1 $num_of_dirs`;
	do
		DIR_NAME="randSeed$dir_number"
		mkdir -p $DIR_NAME
		cd $DIR_NAME
		cp ~/bin/hedonics .
		cp ~/etc/input.txt .
		sed -i 's|randSeedStart 10|randSeedStart $dir_number|g' input.txt
		sed -i 's|randSeedEnd 10|randSeedStart $dir_number|g' input.txt
	done
}


#
# execute the code giving one argument, the path of the top level directory of the simulation that contains the randSeed directories.
# 
execute-seed-dirs () {
	cd $1
	for seed_dir in `ls -1 randSeed*`;
	do
		cd $1/$seed_dir
		./hedonics input.txt
	done
}


#
# collect the results and put them in a single file to analyze. 
# takes one argument the name of the file to put the results in.
#
collect-results () {

	OUTPUT=$1
	for seed_dir in `ls -1 randSeed*`;
	do
		cat $seed_dir/data1.txt >> $OUPUT
	done
}


#
# Change the value for the firmNum in the input files within 
# the randSeed directories
# takes one argument the new value to change firmNum to. 
#

change-firm-num () {
	newValue=$1
	for seed_dir in `ls -1 randSeed*`;
	do
		firmNumber=`grep firmNum input.txt`
		sed -i "s/$firmNumber/firmNum $newValue/g" input.txt
	done
}
