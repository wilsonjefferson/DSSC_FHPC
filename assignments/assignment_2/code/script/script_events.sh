#		F.H.P.C. Assingment 2

#		@file prefix_sum_timer.cc
#		@author Pietro Morichetti
#		@date 17/12/2019
#		@version 1.1 

# COLLECT A SPECIFIC SUBSET OF PERF EVENTS
# This script takes in input the number of repetitions needed to restrict
# the measurement errors, which are done during the execution of the program.
# Moreover, it takes the integer value as the discriminator of the file .txt
# that collects all the data. The output is a .txt file and contains several
# hardware events.

typeset -i repetitions=${1}
typeset -i tmp=${2}
typeset -i check_1=0 # first events encountered parameter
typeset -i check_2=0 # type of analysis parameter
typeset -i procs=1

if (( repetitions==0 )) ; then # default number of repetitions
	repetitions=20
fi

while (( procs<=20 )) # all the subset of threads
do
	case ${check_2} in
	0)
			perf stat -r ${repetitions} -e cpu-cycles:u,instructions:u,branch-instructions:u,branch-misses:u,cache-misses:u,cache-references:u,bus-cycles:u ./a.out ${procs} 2> ./a.txt ;;
	1)
			perf stat -r ${repetitions} -e L1-dcache-loads:u,L1-dcache-load-misses:u,L1-dcache-stores:u,L1-dcache-store-misses:u ./a.out ${procs} 2> ./a.txt ;;
	2)
			perf stat -r ${repetitions} -e LLC-loads:u,LLC-load-misses:u,LLC-stores:u,LLC-store-misses:u ./a.out ${procs} 2> ./a.txt ;;
  3)
			perf stat -r ${repetitions} -e dTLB-loads:u,dTLB-load-misses:u,dTLB-stores:u,dTLB-store-misses:u ./a.out ${procs} 2> ./a.txt ;;
	esac
	
	if (( check_1==0 )) ; then # new events encountered, catch the names
		cat a.txt | grep - | sed '$d' | tr --squeeze-repeats ' ' > c.txt
		sed -i 's/<not supported>/NA/g' c.txt # if not exist, then replace with 'NA'
		cat c.txt | cut -d" " -f3 > ./b.txt
		fields=$(paste -s -d" " ./b.txt)
	fi
  
  # catch the events value
  cat a.txt | grep - | sed '$d' | tr --squeeze-repeats ' ' > c.txt
	sed -i 's/<not supported>/NA/g' c.txt # if not exist then replace with 'NA'
	cat c.txt | cut -d" " -f2 > ./b.txt
  values=$(paste -s -d" " ./b.txt)
  
  # paste different field, from different analys, in different files
  if (( check_2==0 )) && (( check_1==0 )) ; then
		echo ${fields} > ./e_${tmp}_1.txt
		check_1=1
	elif (( check_2==1 )) && (( check_1==0 )) ; then
		echo ${fields} > ./e_${tmp}_2.txt
		check_1=1
	elif (( check_2==2 )) && (( check_1==0 )) ; then
		echo ${fields} > ./e_${tmp}_3.txt
		check_1=1
	elif (( check_2==3 )) && (( check_1==0 )) ; then
		echo ${fields} > ./e_${tmp}_4.txt
		check_1=1
  fi
  
  # paste different events value, from different analysis, in specific files
  if (( check_2==0 )) && (( check_1==1 )) ; then
		echo ${values} >> ./e_${tmp}_1.txt
	elif (( check_2==1 )) && (( check_1==1 )) ; then
		echo ${values} >> ./e_${tmp}_2.txt
	elif (( check_2==2 )) && (( check_1==1 )) ; then
		echo ${values} >> ./e_${tmp}_3.txt
	elif (( check_2==3 )) && (( check_1==1 )) ; then
		echo ${values} >> ./e_${tmp}_4.txt
  fi
	
	# if the analysis of a subset of events is finished then change subset
	if (( procs==20 )) && (( check_2==0 )) ; then
		check_1=0
		procs=1
		check_2=1
	elif (( procs==20 )) && (( check_2==1 )) ; then
		check_1=0
		procs=1
		check_2=2
	elif (( procs==20 )) && (( check_2==2 )) ; then
		check_1=0
		procs=1
		check_2=3
	elif (( procs==20 )) && (( check_2==3 )) ; then
		check_1=0
		procs=1
		check_2=4
	else
		(( procs++ ))
	fi
done

# merge the different files in a columns view
paste -d" " ./e_${tmp}_1.txt ./e_${tmp}_2.txt ./e_${tmp}_3.txt ./e_${tmp}_4.txt > ./a.txt
sed -i 's/,/./g' ./a.txt
cat a.txt | column -t > ./file_events_tot_${tmp}.txt
rm ./a.txt ./b.txt ./c.txt ./e_${tmp}_1.txt ./e_${tmp}_2.txt ./e_${tmp}_3.txt ./e_${tmp}_4.txt
