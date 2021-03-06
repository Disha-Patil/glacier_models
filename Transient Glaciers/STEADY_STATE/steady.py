import os
from os import system as trm

trm("rm -r ../output_steady")
trm("mkdir ../output_steady")

with open('tau.txt') as g: 
    ela_arr = [[j[0], float(j[1]), float(j[2]), float(j[3]), float(j[4])] for j in [k.split() for k in g.readlines()]]
i=0
for t in ela_arr:
    print(i+1)
    glid, tau, ela, beta, C = t
    if(tau==0):continue
    tmax = min([6000,1000+int(tau)])

    #COPY AND MAKE BED AND ELA FILES
    trm("cp ../bed/bed_%s bed.txt" %(glid))
    trm("cp ../bed/bed_%s smoothbed.txt" %(glid))
    trm("awk '{if($3>0){print $1,$2,%f}else{print $1,$2,0}}' bed.txt > ela.txt "%(ela))
    trm("python smoothen.py")

    #CALCULATE XMAX YMAX AND NNZ
    with open('bed.txt') as b: bed_arr = [[int(j[0]), int(j[1]), float(j[2])] for j in [k.split() for k in b.readlines()]]
    xmax , ymax = bed_arr[-1][0]+1, bed_arr[0][1]+1
    nnz = 0
    for c in bed_arr: 
        if(c[2]):nnz+=1
    #EDIT GLOBALS FILE
    with open('SIA_CODES/globals.h') as glb: lines = [l for l in glb.readlines()]
    temp_file = open('temp','w')
    for l in lines: 
        if(len(l)>10 and l[8:10]=="xm"):
            temp_file.write("#define xmax %d\n"%(xmax))
        elif(len(l)>12 and l[8:12]=="tmax"):
            temp_file.write("#define tmax %d\n"%(tmax))
        elif(len(l)>10 and l[8:10]=="ym"):
            temp_file.write("#define ymax %d\n"%(ymax))
        elif(len(l)>10 and l[8:10]=="n0"):
            temp_file.write("#define n0 %d\n"%(nnz))
        elif(len(l)>11 and l[8:11]=="ela"):
            temp_file.write("#define ela_change %d\n"%(0))
        elif(len(l)>11 and l[8:12]=="mbgr"):
            temp_file.write("#define mbgrad %.3f\n"%(beta))
        elif(len(l)>11 and l[8:11]=="C 0"):
            temp_file.write("#define C %.7f\n"%(C))
        elif(len(l)==0):continue
        else: temp_file.write(l)
    temp_file.close()
    trm("mv temp SIA_CODES/globals.h")
    trm("make")
    trm("./tune")

    #trm("mv SIA_CODES/ice.txt ice.txt")
    with open('replot.sh') as glb: lines = [l for l in glb.readlines()]
    temp_file = open('temp','w')
    for l in lines:
        if(l[0]== "p"):
            temp_file.write("p 'ice.txt' u 1:5 w l lw 1 ti '%s';\n"%(glid))
        elif(l[0:5]== "set x"):
            temp_file.write("set xrange[%d:];\n"%(tmax-100))
        else:temp_file.write(l)
    temp_file.close()
    trm("mv temp replot.sh")
    trm("chmod 777 replot.sh")
    trm("./replot.sh")
    
    #OUTPUT FILES
    trm("mv h_steady.txt ../output_steady/hinit_%s"%(glid))
    trm("mv ice.png ../output_steady/ice_%s.png"%(glid))
    i=i+1


trm("rm -rd bed.txt")
trm("rm -rd ice.txt")
trm("rm -rd ela.txt")
trm("rm -rd smoothbed.txt")
trm("rm -rf h_steady.txt")
trm("rm -rf ice.png")
