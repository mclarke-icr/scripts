ABC = [ line.rstrip() for line in open("ABC_samples.txt","r") ]
AB = [ line.rstrip() for line in open("AB_samples.txt","r") ]
A = [ line.rstrip() for line in open("A_samples.txt","r") ]

out = open("samples.txt","w")

for samp in ABC:
    if samp in A:
        out.write(samp+"\tA\n")
    elif samp in AB:
        out.write(samp+"\tB\n")
    else:
        out.write(samp+"\tC\n")
        
        
epk = [ line.split("\t") for line in open("Exome_Processing_Key_20150707.txt","r")]
sanger = [ line.split("\t") for line in open("TrioSangerAnnotation_20151202.txt","r") ]
samp = [ line.split("\t")[0] for line in open("../opex/COGprobands/samples.txt","r") ]
out=open("TSA_20151202.txt","w")
for sl in sanger :
    for line in epk: 
        if line[11] == sl[0]:
            if line[16] == "check": 
                if line[2] in samp:
                    out.write("\t".join([sl[0],line[2],sl[1],sl[2],sl[3]]))
                    
                    

sanger = [ line.split("\t") for line in open("/scratch/cancgene/mclarke/REF/AllSangerAnnotation_20151203_ER.txt","r") ]
samp = [ line.split("\t")[0] for line in open("/scratch/cancgene/mclarke/opex/COGprobands/samples.txt","r") ]
out=open("/scratch/cancgene/mclarke/REF/ProbandSangerAnnotation_20151203_MC.txt","w")
for sl in sanger :
    if sl[1] in samp:
        out.write("\t".join(sl))
        
        
import os 

for root, dirs, files in os.walk("/scratch/cancgene/mclarke/opex/opex-v1.0.0", topdown=True):
    for name in files:
        print(os.path.join(root, name))
    for name in dirs:
        print(os.path.join(root, name))
    