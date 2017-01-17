#!/usr/bin/python
import os,inspect
homedir=os.getcwd()
callfile=os.path.abspath(os.path.dirname(inspect.stack()[0][1]))
print callfile
filelist=os.listdir(callfile)
sourcename=[elem for elem in filelist if elem[-3:]=="f90"]
sourcename.sort()
sourcename=sourcename[-1]
sourcefile=callfile+"/"+sourcename
print sourcefile
execute="jq_exe"
os.system("/share/apps/openmpi-intel/default/bin/mpif90 -o "+homedir+"/"+execute+" -O3 "+sourcefile)
jobfile=homedir+"/job.sh"
f=open(jobfile,"w")
f.write("#!/bin/sh\n"+"#PBS -N JQ_"+homedir.split("/")[-1]+"\n")
f.write("#PBS -q production\n")
f.write("#PBS -l select=16:ncpus=1:mpiprocs=1\n")
f.write("#PBS -l place=free\n")
f.write("#PBS -V\n")
f.write("#PBS -r n\n")
f.write("#PBS -o "+homedir+"/Output\n")
f.write("#PBS -e "+homedir+"/Error\n")
f.write("cd "+homedir+"\n")
f.write("cat $PBS_NODEFILE > "+homedir+"/nodes\n")
f.write("/share/apps/openmpi-intel/default/bin/mpirun  ./"+execute) 
f.close()
os.system("qsub "+homedir+"/job.sh")




