#!/bin/bash
#$ -S /bin/bash
#$ -o /stratt/andrew/ReadConfig/output
#$ -wd /stratt/andrew/ReadConfig
#$ -b y
#$ -j y
#$ -N ZipGenCfg
#$ -m ae
#$ -M atatat123@gmail.com
#$ -cwd
#$ -v LD_LIBRARY_PATH="/stratt/vale/toor/lib:/stratt/vale/toor/gcc62/lib64"
#

task_num=$SGE_TASK_ID

num_degrees=$(awk -v var=$task_num 'NR==var{print $2"degrees"}' inputfiles/input.data)

/stratt/andrew/ReadConfig/Trajectory.exe -c #| xz >  \
#"/stratt/andrew/ReadConfig/Trajectory/$num_degrees/Trajectory$SGE_TASK_ID.xz"
