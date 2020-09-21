# MACHINE SPECIFICATIONS

## test machine (local environment)
Acer Aspire E 15, E5-571G-597D, Intel Core i5-4210U, 1.70 GHz
https://www.acer.com/ac/en/LA/content/model/NX.MLZST.001

## effective machine
Ulisse

# COMMENTS ON THE REPORT

- all analysis and graphs have been done with all cores
  avaiable, therefore the information on the subset with 
  20 threads is misrepresented because at least one core 
  is used by the S.O.
- In the exercise_0*_datasets/events there are some events
  without integer value, instead they have a 'NA' parameter
  i.e. that Ulisse's perf cannot perfom these type of analysis
  (at the THIS moment).

# EXERCISE 0

## compile
 gcc cleaned_t_b_one.c -std=c99 -fopenmp
 gcc cleaned_t_b_all.c -std=c99 -fopenmp

 - optional question
   c++ cleaned_t_b_all_option_1.cc -fopenmp
   c++ cleaned_t_b_all_option_2.cc -fopenmp

 - IPC vs Transfer rate
   c++ cleaned_t_b_one_timer.cc -std=c++0x // de-comment some timers to count another region of code
   c++ cleaned_t_b_all_timer.cc -std=c++0x // de-comment some timers to count another region of code

## analysis
 ./script_events.sh <number of repetition> <discriminator>
 ./script_timer.sh <number of repetition> <discriminator>

# EXERCISE 3

## compile
 serial = c++ prefix_sum.cc -std=c++0x
 parallel = c++ prefix_sum.cc -std=c++0x -fopenmp

 - IPC vs Transfer rate
   serial = c++ prefix_sum_timer.cc -std=c++0x // de-comment some timers to count another region of code
   parallel = c++ prefix_sum_timer.cc -std=c++0x -fopenmp // de-comment some timers to count another region of code

## analysis
./script_events.sh <number of repetition> <discriminator>
./script_timer.sh <number of repetition> <discriminator>
