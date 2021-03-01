import sys
import os
import subprocess
import argparse
from collections import OrderedDict
from multiprocessing import Pool

class RunParallelBcftools:

    def __init__(self,bam,reference,idx,cores,bins):
        self.bam=bam
        self.outDir=bam.replace(".txt","")
        self.reference=reference
        self.idx=idx
        self.cores=cores
        self.bins=bins
        self.inputList=[]
        self.outputDict=OrderedDict()
        self.withinChrmDict=OrderedDict()
        self.globalLog=open(self.bam+"GlobalLog.txt","w")


    def ReadIdx(self):
        newDirCommand=("mkdir ./bcftools_"+self.outDir)
        subprocess.call([newDirCommand],shell=True)
        with open(self.idx) as source:
            for line in source:
                a=line.strip().split()
                self.chrm=a[0]
                newSubDirCommand=("mkdir ./bcftools_"+self.outDir+"/"+a[0])
                subprocess.call([newSubDirCommand],shell=True)
                chrmBins=int(a[1])//int(self.bins)
                self.start=1-int(self.bins)
                self.end=0
                self.count=0
                if int(a[1])%int(self.bins)!=0:
                    chrmBins+=1
                self.outputDict[self.chrm]=int(chrmBins)
                for i in range(chrmBins):
                    self.count+=1
                    self.start+=int(self.bins)
                    self.end+=int(self.bins)
                    if self.end>int(a[1]):
                        self.end=int(a[1])
                    tmpList=[self.reference,self.bam,self.chrm,self.start,self.end,self.count]
                    self.inputList.append(tmpList[:])
                    del tmpList[:]
        return self.inputList

    def MergeWithinChrm(self):
        for chrom in self.outputDict:
            self.withinChrmDict[chrom]=[]
            self.directory="./bcftools_"+self.outDir+"/"+chrom
            for i in range(self.outputDict[chrom]):
                fileName=chrom+"_"+str(i+1)+".bcf"
                filePath=os.path.join(self.directory,fileName)
                if os.path.isfile(filePath) and os.access(filePath,os.R_OK):
                    self.globalLog.write(fileName+" exists"+"\n")
                    indexCommand=("bcftools index "+filePath)
                    subprocess.call([indexCommand],stderr=self.globalLog,shell=True)
                    self.withinChrmDict[chrom].append(filePath)
                else:self.globalLog.write(fileName+" does not exist"+"\n")
        return self.withinChrmDict

    def RunParallel(self):
        self.inputList=self.ReadIdx()
        self.cores=round(len(self.inputList)/int(self.cores))
        acrossMerge=[]
        with Pool(processes=len(self.inputList)) as pool:
            pool.map(RunBcftools,self.inputList,int(self.cores))
        self.withinChrmDict=self.MergeWithinChrm()
        for chrm in self.withinChrmDict:
            mergeCommand=("bcftools concat "+" ".join(self.withinChrmDict[chrm])+" -Ob -o "+"./bcftools_"+self.outDir+"/"+chrm"_merged.bcf")
            subprocess.call([mergeCommand],stderr=self.globalLog,shell=True)
            indexCommand=("bcftools index "+"./bcftools_"+self.outDir+"/"+chrm+"_merged.bcf")
            acrossMerge.append("./bcftools_"+self.outDir+"/"+chrm+"_merged.bcf")
            subprocess.call([indexCommand],stderr=self.globalLog,shell=True)
        genomeMergeCommand=("bcftools concat "+" ".join(acrossMerge)+" -Ob -o "+"./bcftools_"+self.outDir+"/"+self.outDir+"_meged.bcf")
        subprocess.call([genomeMergeCommand],stderr=self.globalLog,shell=True)

def RunBcftools(argList):
    #bcftools mpileup -Ou -b AllBamFiles.txt -q 20 -Q 20 -C 50 -I -r NC_019458.2:1-1000000 -f ~/data/Shared/References/Oar_v4./GCF_000298735.2_Oar_v4.0_genomic.fna|bcftools call -mv -Ob -o NC_019458.2.1.bcf
    outDir=argList[1].replace(".txt","")
    fileGatkPath="./bcftools_"+outDir+"/"+argList[2]
    dest1=fileGatkPath+"/"+argList[2]+"_"+str(argList[-1])+".bcf"
    dest=open(str(argList[1])+"GlobalLog.txt","a")
    bash_command=('bcftools mpileup -Ou -b '+argList[1]+' -q 20 -Q 20 -C 50 -I -a "AD,ADF,ADR,DP,SP,INFO/AD,INFO/ADF,INFO/ADR"-p -r '+argList[2]+":"+str(argList[3])+"-"+str(argList[4])+" -f "+argList[0]+"| bcftools call -mv -Ob -o "+dest1)
    #print(bash_command)
    subprocess.call([bash_command],stderr=dest,shell=True)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="This python script will run the variant calling using samtools in parallel",pilog="author: Maulik Upadhyay (Upadhyaya.maulik@gmail.com)")
    parser.add_argument('-b',"--bFile",metavar="File",help="bam file",required=True)
    parser.add_argument('-r',"--rFile",metavar="File",help="reference",required=True)
    parser.add_argument("-i","--iFile",metavar="File",help="index file",required=True)
    parser.add_argument("-c","--nCore",metavar="Int",help="number of cores",required=True)
    parser.add_argument("-n","--nBins",metavar="Int",help="bins interval for parallel processes",required=True)
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        runParGatk=RunParallelBcftools(args.bFile,args.rFile,args.iFile,args.nCore,args.nBins)
        runParGatk.RunParallel()
