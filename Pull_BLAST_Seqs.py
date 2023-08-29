import pandas as pd
from Bio import SeqIO
import sys

def pullseq(blastresults, multifasta, outfile):
    df=pd.read_csv(blastresults, sep='\t')
    df.columns=["Sample", "Start", "Stop"]
    seqlist=[]
    names=df["Sample"]
    namesset=list(set(names))
    file=open(multifasta, 'r')
    Lines=file.readlines()
    linecount=0
    SeqDict={rec.id : rec.seq for rec in SeqIO.parse(multifasta, "fasta")}
    # print(SeqDict)
    counter=0
    for i in range(len(names)):
        string=SeqDict[names[i]]
        start=df.iloc[i][1]
        stop=df.iloc[i][2]
        if start<stop:
            seq=string[start-1:stop-1]
            seqlist.append(seq)
        else:
            seq=string[stop-1:start-1]
            seq=seq.reverse_complement()
            seqlist.append(seq)


    df["Seq"]=seqlist

    realstart=[]
    realstop=[]
    finaldf=pd.DataFrame()
    for i in range(len(df)):
        start=df.iloc[i][1]
        stop=df.iloc[i][2]
        if start>stop:
            realstart.append(stop)
            realstop.append(start)
        else:
            realstart.append(start)
            realstop.append(stop)
    df["AdjustedStart"]=realstart
    df["AdjustedStop"]= realstop
    for j in range(len(namesset)):
        overlap=[]
        currentsamp=namesset[j]
        dfname=df[df.values==currentsamp]
        index=dfname.index.values.tolist()
        dfname.sort_values('AdjustedStart')
        firstset=set(list(range(int(dfname.iloc[0][4]),int(dfname.iloc[0][5]))))
        biggestset=firstset
        low=int(dfname.iloc[0][4])
        high=int(dfname.iloc[0][5])
        for k in range(len(dfname)):
            newstart = int(dfname.iloc[k][4])
            newstop = int(dfname.iloc[k][5])

            currset=set(list(range(newstart, newstop)))
            if(len(currset) >= len(biggestset) and newstart<=low and newstop>=high):
                low=newstart
                high=newstop
                biggestset=currset


    df.to_csv(outfile,sep=",", header=True)
def fastamaker(filein,outfasta):
    df=pd.read_csv(filein, sep=",", usecols=["Sample", "Seq"])
    sampnames=df["Sample"].tolist()
    seqs=df["Seq"].tolist()
    with open(outfasta, "w+") as f:
        for i in range(len(sampnames)):
            currentsamp=sampnames[i]
            currentseq=seqs[i]
            f.writelines(">"+currentsamp+"\n")
            f.writelines(currentseq+"\n")


if __name__ == '__main__':
    pullseq(sys.argv[1], sys.argv[2], sys.argv[3]) #1 is the blast output tsv file, 2 is the multifasta used in the database, 3 is the output csv file
    fastamaker(sys.argv[3], sys.argv[4]) # 3 is the output csv from above, and 4 is the output fasta file
