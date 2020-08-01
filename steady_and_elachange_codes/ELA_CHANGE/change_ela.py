import os
from os import system as trm

trm("rm -r ../output_change_ela")
trm("mkdir ../output_change_ela")

with open('tau.txt') as g: 
    ela_arr = [[j[0], float(j[1]), float(j[2]), float(j[3]), float(j[4])] for j in [k.split() for k in g.readlines()]]

for t in ela_arr:
    glid, tau, ela, beta, C = t
    if(tau==0):continue

    #COPY AND MAKE BED AND ELA FILES
    trm("cp ../bed/bed_%s bed.txt" %(glid))
    trm("cp ../bed/bed_%s smoothbed.txt" %(glid))
    trm("cp ../output_steady/hinit_%s h_init_file.txt" %(glid))
    trm("awk '{if($3>0){print $1,$2,%f}else{print $1,$2,0}}' bed.txt > ela.txt "%(ela))

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
            temp_file.write("#define tmax %d\n"%(1000))
        elif(len(l)>10 and l[8:10]=="ym"):
            temp_file.write("#define ymax %d\n"%(ymax))
        elif(len(l)>10 and l[8:10]=="n0"):
            temp_file.write("#define n0 %d\n"%(nnz))
#        elif(len(l)>11 and l[8:11]=="ela"):
#            temp_file.write("#define ela_rate %d\n"%(0.005))
        elif(len(l)>11 and l[8:12]=="mbgr"):
            temp_file.write("#define mbgrad %.3f\n"%(beta))
        elif(len(l)>11 and l[8:11]=="C 0"):
            temp_file.write("#define C %.7f\n"%(C))
        elif(len(l)==0):continue
        else: temp_file.write(l)
    temp_file.close()
    trm("python smoothen.py")
    trm("mv temp SIA_CODES/globals.h")
    trm("make")
    trm("./tune")

    '''with open('replot.sh') as glb: lines = [l for l in glb.readlines()]
    temp_file = open('temp','w')
    for l in lines:
        if(l[0]== "p"):
            temp_file.write("p 'ice.txt' u 1:5 w l lw 1 ti '%s';\n"%(glid))
        else:temp_file.write(l)
    temp_file.close()
    trm("mv temp replot.sh")
    trm("chmod 777 replot.sh")
    trm("./replot.sh")'''

    # CALCULATE SLOPE

    

    #OUTPUT FILES
    trm("cp ice.txt ../output_change_ela/consrv_%s"%(glid))
    trm("mv va ../output_change_ela/va_%s"%(glid))
    trm("mv h_steady.txt ../output_change_ela/hfinal_%s"%(glid))

trm("rm -rd bed.txt")
trm("rm -rd ice.txt")
trm("rm -rd ela.txt")
trm("rm -rd smoothbed.txt")
trm("rm -rf h_steady.txt")
trm("rm -rf ice.png")
trm("rm -rd h_init_file.txt")
