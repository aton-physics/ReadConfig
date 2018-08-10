#!/bin/bash
#$ -S /bin/bash
#$ -o /home/aton/ReadConfig/output
#$ -wd /home/aton/ReadConfig
#$ -b y
#$ -j y
#$ -N ZipGenCfg
#$ -m ae
#$ -M atatat123@gmail.com
#$ -cwd
#$ -v LD_LIBRARY_PATH="/stratt/vale/toor/lib:/stratt/vale/toor/gcc62/lib64"
#

/home/aton/ReadConfig/Trajectory.exe -s | xz > "/home/aton/ReadConfig/Trajectory/75degrees/Trajectory$SGE_TASK_ID.xz"
