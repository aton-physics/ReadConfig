#!/bin/bash
#$ -S /bin/bash
#$ -o /stratt/andrew/ReadConfig/output
#$ -wd /stratt/andrew/ReadConfig
#$ -b y
#$ -j y
#$ -N NoSave
#$ -m ae
#$ -M atatat123@gmail.com
#$ -cwd
#$ -v LD_LIBRARY_PATH="/stratt/vale/toor/lib:/stratt/vale/toor/gcc62/lib64"
#

./Trajectory.exe -s | ./ReadConfig.exe
