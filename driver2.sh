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

unxz -c "/stratt/andrew/ReadConfig/Trajectory/75degrees/Trajectory$SGE_TASK_ID.xz" | ./ReadConfig.exe 
