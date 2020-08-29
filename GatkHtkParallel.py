import sys
import os
import subprocess
import argparse
from collections import OrderedDict
from multiprocessing import Pool

class RunParallelHtk:

    def __init__(self,bam,reference,idx,cores,bins,sample):
        self.bam=bam
        self.reference=reference
        self.idx=idx
        self.cores=cores
        self.bins=bins
        self.sample=sample
        self.inputList=[]
        self.outputDict=OrderedDict()
        self.withinChrmDict=OrderedDict()
        self.acrossGenomeList=[]
        self.globalLog=open(self.sample+"GlobalLog.txt","w")


    def ReadIdx(self):
        newDirCommand=("mkdir ./Gatk_"+self.sample)
        subprocess.call([newDirCommand],shell=True)
        with open(self.idx) as source:
            for line in source:
                a=line.strip().split()
                self.chrm=a[0]
                newSubDirCommand=("mkdir ./Gatk_"+self.sample+"/"+a[0])
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
                    tmpList=[self.reference,self.bam,self.chrm,self.start,self.end,self.count,self.sample]
                    self.inputList.append(tmpList[:])
                    del tmpList[:]
        return self.inputList

    def MergeWithinChrm(self):
        #self.directory="./"
        for chrom in self.outputDict:
            self.withinChrmDict[chrom]=[]
            self.directory="./Gatk_"+self.sample+"/"+chrom
            for i in range(self.outputDict[chrom]):
                fileName=chrom+"_"+str(i+1)+".g.vcf"
                filePath=os.path.join(self.directory,fileName)
                if os.path.isfile(filePath) and os.access(filePath,os.R_OK):
                    self.globalLog.write(fileName+" exists")
                    self.withinChrmDict[chrom].append("-I "+self.directory+"/"+fileName)
                else:
                    self.globalLog.write(fileName+"does not exist")
        return self.withinChrmDict

    def MergeAcrossGenome(self):
        for chrom in self.outputDict:
            fileName=chrom+"_"+"mergedAllPoints.g.vcf"
            directoryG="./Gatk_"+self.sample+"/"+chrom
            filePath=os.path.join(directoryG,fileName)
            if os.path.isfile(filePath) and os.access(filePath,os.R_OK):
                self.globalLog.write(fileName+" exists")
                self.acrossGenomeList.append("-I "+directoryG+"/"+fileName)
            else:self.globalLog.write(fileName+"does not exist")
        return self.acrossGenomeList

    def RunParallel(self):
        self.inputList=self.ReadIdx()
        #print(self.inputList)
        self.cores=round(len(self.inputList)/int(self.cores))
        #print(self.cores)
        with Pool(processes=len(self.inputList)) as pool:
            pool.map(RunGatk,self.inputList,int(self.cores))
        self.withinChrmDict=self.MergeWithinChrm()
        for chrom in self.withinChrmDict:
            inputFilesCommand=" "+" ".join(self.withinChrmDict[chrom])+" "
            #print(inputFilesCommand)
            mergeVcfCommand=('/home/maulik/software/VariantCallingPipeline/SnpCalling/gatk-4.1.8.1/gatk --java-options "-Xmx30G" GatherVcfs -R '+\
            self.reference+inputFilesCommand+"-O ./Gatk_"+self.sample+"/"+chrom+"/"+chrom+"_mergedAllPoints.g.vcf")
            subprocess.call([mergeVcfCommand],stderr=self.globalLog,shell=True)
        self.acrossGenomeList=self.MergeAcrossGenome()
        acrossInFilesCommand=" "+" ".join(self.acrossGenomeList)+" "
        mergeChrmCommand=('/home/maulik/software/VariantCallingPipeline/SnpCalling/gatk-4.1.8.1/gatk --java-options "-Xmx30G" GatherVcfs -R '+\
        self.reference+acrossInFilesCommand+" -O ./Gatk_"+self.sample+"/"+self.sample+"_allChrm.g.vcf")
        print(mergeChrmCommand)
        subprocess.call([mergeChrmCommand],stderr=self.globalLog,shell=True)


def RunGatk(argList):
    #print(argList)
    fileGatkPath="./Gatk_"+argList[-1]+"/"+argList[2]
    dest1=fileGatkPath+"/"+argList[2]+"_"+str(argList[-2])+".g.vcf"
    dest=open(str(argList[-1])+"GlobalLog.txt","a")
    bash_command=('/home/maulik/software/VariantCallingPipeline/SnpCalling/gatk-4.1.8.1/gatk --java-options "-Xmx30G -Djava.library.path=/home/maulik/software/VariantCallingPipeline/SnpCalling/gatk-4.1.8.1/libs -XX:+UseParallelGC -XX:ParallelGCThreads=2" HaplotypeCaller -I '+argList[1]+" -O "+dest1+" -R "+argList[0]+" --sample-name "+argList[-1]+" --emit-ref-confidence GVCF -pairHMM FASTEST_AVAILABLE --native-pair-hmm-threads 2 -L "+argList[2]+":"+str(argList[3])+"-"+str(argList[4]))
    #print(bash_command)
    subprocess.call([bash_command],stderr=dest,shell=True)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="This python script will run the Gatk haplotypecaller in parallel",epilog="author: Maulik Upadhyay (Upadhyaya.maulik@gmail.com)")
    parser.add_argument('-b',"--bFile",metavar="File",help="bam file",required=True)
    parser.add_argument('-r',"--rFile",metavar="File",help="reference",required=True)
    parser.add_argument("-i","--iFile",metavar="File",help="index file",required=True)
    parser.add_argument("-c","--nCore",metavar="Int",help="number of cores",required=True)
    parser.add_argument("-n","--nBins",metavar="Int",help="bins interval for parallel processes",required=True)
    parser.add_argument("-s","--sample",metavar="Str",help="name of the sample (as mentioned in bam header",required=True)
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        runParGatk=RunParallelHtk(args.bFile,args.rFile,args.iFile,args.nCore,args.nBins,args.sample)
        runParGatk.RunParallel()
