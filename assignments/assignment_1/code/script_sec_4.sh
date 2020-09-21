#		F.H.P.C Assignment 1
#		@file script_sec_4.sh		
#		@author Pietro Morichetti
#		@date 08/11/2019
#		@version 1.1 


# COLLECT ALL THE DATA FROM MPI_SUM PROGRAM
# If compare some "invisible" character, it means that you are running 
# the script with the wrong <<characters format>>; so UTF8 - linux end line

typeset -i N=${1} # dimension of the problem
typeset -i N_var=${1}
typeset -i tmp=${2} # integer value to discern several files
typeset -i procs=1 # start from one process
typeset -i check=0
typeset -i i=1

if [ -e ./file_strong_${tmp}.txt ] ; then # delete the old file
 rm ./file_strong_${tmp}.txt
fi

if [ -e ./file_weak_${tmp}.txt ] ; then # delete the old file
 rm ./file_weak_${tmp}.txt
fi

while (( procs<=16 )) # run from 1 to 16 processes
do
	while (( i<=2 )) # make two times to collect strong and weak data with the same conditions
	do
		if (( check==0 )) # first run in the strong way
		then
				echo ${N} > ./number_operation.txt # overwrite the size of the problem
				#echo "N = ${N}"	
				/usr/bin/time mpirun -np ${procs} ./mpi_sum.x 2> ./a.txt
		else # second run in the weak way
				echo ${N_var} > ./number_operation.txt 
				/usr/bin/time mpirun -np ${procs} ./mpi_sum.x 2> ./a.txt
				N_var=$((${N_var}*(${procs}+1)/${procs})) # increases the size in proportion to the increase in processes
				#echo "N_var = ${N_var}"	
		fi
 	
 	 #cat ./a.txt
   echo "-------------------"	 # divider the strong to the weak data
   
   # collect all the data from the output in the matrix type
	 elapsedtime=$(cat ./a.txt | grep elapsed | cut -d" " -f3 | cut -d"e" -f1 | cut -d":" -f2)
	 usertime=$(cat ./a.txt | grep user | cut -d"u" -f1)
	 systemtime=$(cat ./a.txt | grep system | cut -d" " -f2 | cut -d"s" -f1)
	 cpuwork=$(cat ./a.txt | grep CPU | cut -d" " -f4 | cut -d"%" -f1)
	 
	 # general row of the matrix
	 row="${procs} ${usertime} ${systemtime} ${elapsedtime} ${cpuwork}"
	 
	 if ((check==0)) # write the row into the correct file and change the check
	 then
				echo ${row} >> file_strong_${tmp}.txt
				check=1
		else
				echo ${row} >> file_weak_${tmp}.txt
				check=0
	 fi
	 (( i++ )) 
 	 #echo "check = ${check}"
 done
 rm ./a.txt # remove the temporal file
 (( procs++ ))
 i=1
done
