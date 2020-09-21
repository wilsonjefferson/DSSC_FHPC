#		F.H.P.C. Assingment 2

#		@file prefix_sum_timer.cc
#		@author Pietro Morichetti
#		@date 17/12/2019
#		@version 1.1 

# COLLECT TIMES AND ERRORS
# this script takes in input a integer value (called the "discriminator") to distinguish 
# among different execution of itself, an execution example is: ./script_elapsed 1. The output
# is a ./file_<discriminator>.txt that store the elapsed time and the percentage error.

typeset -i tmp=${1}
typeset -i procs=1

if [ -e ./file_${tmp}.txt ] ; then # WARNING: all the .txt with this name will be overwritten!
 rm ./file_${tmp}.txt
fi

while (( procs<=20 )) # all the threads subset
do
	perf stat -r 100 -e cycles ./a.out ${procs} 2> ./a.txt    # various perf statistic
	
	elapsedtime=$(cat ./a.txt | grep seconds | tr --squeeze-repeats ' ' | cut -d" " -f2)
	percent_error=$(cat ./a.txt | grep seconds | tr --squeeze-repeats ' ' |  cut -d" " -f8 | cut -d"%" -f1)
	
	row="${procs} ${elapsedtime} ${percent_error}" # build the row
	echo ${row} >> file_${tmp}.txt # put the row into the file
	
	rm ./a.txt # delete support file
 (( procs++ ))
done
