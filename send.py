#!/usr/bin/python
import random
import os, sys
# IsCluster=True
IsCluster=False
sourcedir="./program/"
execute="XXZ.exe"
Dim = 3
Latticename = "Pyrochlore"
#Latticename = "Cubic"
Nk = 6
homedir=os.getcwd()
filelist=os.listdir(sourcedir)
sourcename=[elem for elem in filelist if elem[0:3]=="XXZ" and elem[-3:]=="f90"]
sourcename.sort()
sourcename=sourcename[-1]
os.system("gfortran "+sourcedir+"/"+sourcename+" -O3 -o "+homedir+"/"+execute)
infilepath=homedir+"/infile"
outfilepath=homedir+"/outfile"
jobfilepath=homedir+"/jobfile"
inlist=open(homedir+"/inlist","r")
i=0
nblck=int(inlist.readline().split(":")[1])
nsample=int(inlist.readline().split(":")[1])
nsweep=int(inlist.readline().split(":")[1])
nsave=int(inlist.readline().split(":")[1])
ntoss=int(inlist.readline().split(":")[1])
isload=int(inlist.readline().split(":")[1])
seed=-int(random.random()*1000)
if(os.path.exists(infilepath)!=True):
    os.system("mkdir "+infilepath)
if(os.path.exists(outfilepath)!=True):
    os.system("mkdir "+outfilepath)
if(os.path.exists(jobfilepath)!=True):
    os.system("mkdir "+jobfilepath)
for eachline in inlist:
        i+=1
        para=eachline.split()
        for j in range(int(para[-1])):
            seed-=1
            infile="_in"+str(i)+"_"+str(j)
            outfile="_out"+str(i)+"_"+str(j)
            jobfile="_job"+str(i)+"_"+str(j)+".sh"
            f=open(infilepath+"/"+infile,"w")
            item=para[0:-1]
            item.append(str(ntoss))
            item.append(str(nsample))
            item.append(str(nsweep))
            item.append(str(nsave))
            item.append(str(seed))
            item.append(str(nblck))
            item.append(str(isload))
            item.append("_".join(para[2:-1])+"_"+str(j)+".cnf")
            item.append("static_corr_"+str(i)+"_"+str(j)+".txt")
            item.append("mid_hs_sqa0_"+str(i)+"_"+str(j)+".txt")
            for k in range(Nk):
                item.append("\n")
                item.append("corr_k"+str(k+1)+"_"+str(i)+"_"+str(j)+".txt")
            for k in range(Nk):
                item.append("\n")
                item.append("corr_k"+str(k+1)+"_tau_"+str(i)+"_"+str(j)+".txt")
            stri=" ".join(item)
            f.write(str(Dim)+"\n")
            f.write(Latticename+" "+stri)
            f.close()

            if IsCluster==False:
                os.system("./"+execute+" < "+infilepath+"/"+infile+" > "+outfilepath+"/"+outfile+" &")
                #os.system("./"+execute+" < "+infilepath+"/"+infile)
            else:
                with open(jobfilepath+"/"+jobfile, "w") as fjob:
                    fjob.write("#!/bin/sh\n"+"#PBS -N "+jobfile+"\n")
                    fjob.write("#PBS -o "+homedir+"/Output\n")
                    fjob.write("#PBS -e "+homedir+"/Error\n")
                    fjob.write("#PBS -l walltime=200:00:00\n")
                    fjob.write("echo $PBS_JOBID >>"+homedir+"/id_job.log\n")
                    fjob.write("cd "+homedir+"\n")
                    fjob.write("./"+execute+" < "+infilepath+"/"+infile+" > "+outfilepath+"/"+outfile)

                os.system("qsub "+jobfilepath+ "/"+jobfile)
                os.system("rm "+jobfilepath+ "/"+jobfile)
print("Jobs manage daemon is ended")
sys.exit(0)
