#!/usr/bin/python
import os

FileName="./pyrochlo.dat"
Head="pyrochlo"
Parameters=[1,4,5,6,7]
DataRank={15:"M2", 25:"C"}
#DataRank={22:"Spatial_Sus",23:"Total_Sus",24:"Corr"}
#DataRank={10:"VBS",11:"Neel",27:"Q_VBS",28:"Q_Neel"}

############   Read data from file and cut them into blocks   ####################

f=open(FileName,"r")
DataRankList=DataRank.keys()
print DataRankList
print DataRank.values()
#Divide the data file into several data blocks
Blocks=[elem.strip() for elem in f.read().split(Head)]
Blocks=filter(lambda x:x!="",Blocks)
#print Blocks[0]
f.close()

IndexMap=[]
for line in Blocks[0].split("\n"):
    try:
        if(int(line.split()[0]) in DataRankList):
            IndexMap.append(int(line.split()[0]))
    except:
        pass
#print IndexMap

#############  Extract parameters and variables from data blocks    ###############

ParaList=[]
WeightList=[]
DataList=[]
for block in Blocks:
    #extract the parameters list from all the data blocks
    blocklist=block.split("\n")
    #print blocklist[0].split()[-1]
    ParaList.append(tuple([float(blocklist[0].split()[i-1]) for i in Parameters]))
    #WeightList.append(float(blocklist[0].split()[-1]))
    WeightList.append(1.0)
    DataListTemp=[]
    #extract the datas asked by DataRank from all the data blocks
    
    for line in blocklist[1:]:
        try:
            if(int(line.split()[0]) in DataRankList):
                DataListTemp.append([float(line.split()[1]),float(line.split()[2])])
        except:
            pass
    DataList.append(DataListTemp)
#print ParaList
ParaSet=sorted(list(set(ParaList)))
#print ParaSet

#############   Merge data blocks  and create DataMap     ##################

#DataMap={Parameters:Data}
DataMap={}

for paras in ParaSet:
    Index=[]
    DataValue=[]
    datas=[]
    for i in range(0,len(ParaList)):
        if(ParaList[i]==paras):
            Index.append(i)
    WeightTemp=[WeightList[j] for j in Index]
    DataTemp=[DataList[j] for j in Index]
    for i in range(0,len(DataTemp[0])):
        mean=0.0
        std=0.0
        for j in range(0,len(DataTemp)):
            mean+=DataTemp[j][i][0]
            std+=DataTemp[j][i][1]**2*WeightTemp[j]
        mean=mean/len(DataTemp)
        std=(std/len(DataTemp)/sum(WeightTemp))**0.5
        datas.append([mean,std])
    #print datas,DataTemp
    DataMap[paras]=datas
    #print DataMap


###########    Output DataMap as the format             ###################
#os.remove("./*.ext")

Switch=True   # Classify data by variables
#Switch=False  #Classify data by parameters and variables


SubParameter=[]
for elem in ParaSet:
    SubParameter.append(elem[0])
SubParameter=list(set(SubParameter))

os.system("rm *.ext")
for para in SubParameter:
    SubDataMap={}
    ParaInFileList=[]
    for elem in ParaSet:
        if(Switch or elem[0]==para):
            #if(elem[0]==elem[1]):     #add some filter here
                paras=(elem[0],elem[1],elem[2],elem[3])
                SubDataMap[paras]=DataMap[elem]

##    for i in range(0,len(DataRank.values())):
##        if(not Switch):
##            OutputName=DataRank.values()[i]+"_"+str(para)+".ext.dat"
##        else:
##            OutputName=DataRank.values()[i]+".ext"
##        f=open(OutputName,"w")
##        for elem in sorted(SubDataMap.keys()):
##            OutStr="   ".join([str(e) for e in elem])+"   "+str(SubDataMap[elem][i][0])+"   "+str(SubDataMap[elem][i][1])+"\n"
##            f.write(OutStr)
##        f.close()
    #print IndexMap
    for i in range(0,len(IndexMap)):
        if(not Switch):
            OutputName=DataRank[IndexMap[i]]+"_"+str(para)+".ext.dat"
        else:
            OutputName=DataRank[IndexMap[i]]+".extcom.dat"
        f=open(OutputName,"w")
        for elem in sorted(SubDataMap.keys()):
            OutStr="   ".join([str(e) for e in elem])+"   "+str(SubDataMap[elem][i][0])+"   "+str(SubDataMap[elem][i][1])+"\n"
            f.write(OutStr)
        f.close()

