#		F.H.P.C. Assingment 2

#		@file prefix_sum_timer.cc
#		@author Pietro Morichetti
#		@date 17/12/2019
#		@version 1.1 

# COLLECT TRANSFER RATE AND TWO PERF EVENTS
#	This script is used to collect information about the transfer rate
# of specific code's portion (MUST put some timer into the program
# and the program have to figure out only the GB/s in putput!),
# and collect data from cpu-cycle and number of instruction through
# perf calls. The Output of script is three .txt that contain the 
# mean value on to all repetition for each subset of threads.

typeset -i repetitions=${1}
typeset -i tmp=${2}
typeset -i repet=1
typeset -i count=0 # avarage parameter
typeset -i total=0 # avarage parameter
typeset -i procs=1

if [ -e ./timer_frequency_${tmp}.txt ] ; then # WARNING: delete the file if it already exists!
 rm ./timer_frequency_${tmp}.txt
fi

if [ -e ./timer_cycle_${tmp}.txt ] ; then # WARNING: delete the file if it already exists!
 rm ./timer_cycle_${tmp}.txt
fi

if [ -e ./timer_ins_${tmp}.txt ] ; then # WARNING: delete the file if it already exists!
 rm ./timer_ins_${tmp}.txt
fi

if (( repetitions==0 )) ; then # default number of repetitions
        repetitions=10
fi

while (( procs<=20 )) # for each subset of threads
do
        while (( repet<=repetitions )) # repet the analysis several times
        do
          perf stat -o ./a.txt -e cpu-cycles:u,instructions:u ./a.out ${procs} > ./l.txt # catch two events and the output that is the transfer rate

          frequence=$(cat ./l.txt)
          if [ "$frequence" = "inf" ] ; then # catch failing test -> repet again
                echo "failed."
          else
                # store the frequency
                echo $frequence >> ./e.txt # store the frequncy of the actual analysis
                ((repet++))  # go on with repetitions                           
				          echo "success."

                # catch the events values
                cat a.txt | grep :u | tr --squeeze-repeats ' ' | cut -d" " -f2 > ./d.txt
                values=$(paste -s -d" " ./d.txt)
                echo ${values} >> ./f.txt
          fi
        done
					
					# average the frequencies of all the repetition of this subset of threads
        awk '{ total += $1; count++ } END { printf("%.3f\n", total/count) }' e.txt >> ./timer_frequency_${tmp}.txt
        count=0
        total=0
					
					 # average the cycles of all the repetition of this subset of threads
        sed -i 's/,//g' ./f.txt
        awk '{ total += $1; count++ } END { printf("%.3f\n", total/count) }' f.txt >> ./timer_cycle_${tmp}.txt
        count=0
        total=0
					 # average the instructions of all the repetition of this subset of threads
        awk '{ total += $2; count++ } END { printf("%.3f\n", total/count) }' f.txt >> ./timer_ins_${tmp}.txt

        rm f.txt e.txt
        count=0
        total=0
        repet=1
        (( procs++ )) # consider new subset of threads
done
rm ./a.txt ./d.txt ./l.txt     