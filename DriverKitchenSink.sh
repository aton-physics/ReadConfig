#!/bin/bash
#$ -S /bin/bash
#$ -o /stratt/andrew/ReadConfig/output
#$ -wd /stratt/andrew/ReadConfig
#$ -b y
#$ -j y
#$ -N KitchenSink
#$ -m ae
#$ -M atatat123@gmail.com
#$ -cwd
#$ -v LD_LIBRARY_PATH="/stratt/vale/toor/lib:/stratt/vale/toor/gcc62/lib64"
#


task_num=$SGE_TASK_ID

num_degrees=$(awk -v var=$task_num 'NR==var{print $2"degrees"}' inputfiles/input.data)

/stratt/andrew/ReadConfig/Trajectory.exe -c
/stratt/andrew/ReadConfig/Trajectory.exe -s | xz > \
"/stratt/andrew/ReadConfig/Trajectory/$num_degrees/Trajectory$SGE_TASK_ID.xz"
unxz -c "/stratt/andrew/ReadConfig/Trajectory/$num_degrees/Trajectory$SGE_TASK_ID.xz" | ./Statics.exe
unxz -c "/stratt/andrew/ReadConfig/Trajectory/$num_degrees/Trajectory$SGE_TASK_ID.xz" | ./Corr.exe

