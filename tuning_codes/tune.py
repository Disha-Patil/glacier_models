import os
from os import system as trm

trm('rm -r ela')
trm('mkdir ela')
with open('guess_ela.txt') as g: ela_arr = [[j[0], float(j[1]), float(j[2]), float(j[3])] for j in [k.split() for k in g.readlines()]]

#make tau
tau_file = open('ela/tau.txt','w')
tau_file.close()

for t in ela_arr:
    glid = t[0]
    #if(glid!="15.05957"):continue
    guess_ela = t[1]
    beta = t[2]
    C = t[3]
    #COPY AND MAKE BED AND ELA FILES
    trm("cp bed/bed_%s bed.txt" %(glid))
    trm("cp bed/bed_%s smoothbed.txt" %(glid))
    trm("awk '{if($3>0){print $1,$2,1}else{print $1,$2,0}}' bed/ice_%s > ice_mask "%(glid))
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
    
    low = guess_ela
    high = guess_ela

    #BISECTION
    first = 1
    while True:
        print(guess_ela)
        trm("awk '{if($3>0){print $1,$2,%f}else{print $1,$2,0}}' bed.txt > ela.txt "%(guess_ela))
        trm("./tune")
        with open('steady_slope') as s: slope, af, tau = [[float(j[0]), float(j[1]), float(j[2])] for j in [k.split() for k in s.readlines()]][0]
        if(first):
            first = 0
            if(slope<=0.000001): case = 1
            elif(slope>0.000001): case = 2
            else: case = 3

        
        # ELA TOO HIGH
        if(case==1): 
            if(slope < 0.000001):
                high = low
                low -= 100
                guess_ela = low
            if(slope > 0.000001):
                print("    %f    %f"%(low,high)) 
                break

        # ELA TOO LOW
        if(case==2):
            if(slope > 0.000001):
                low = high 
                high += 100
                guess_ela = high
            if(slope < 0.000001):
                print("    %f    %f"%(low,high))
                break

        if(case==3):
            high = guess_ela+50
            low = guess_ela-50
            break
        

    while (high-low > 5):
        if(low == high):break
        ela = (high+low)/2
        #print(ela)
        trm("awk '{if($3>0){print $1,$2,%f}else{print $1,$2,0}}' bed.txt > ela.txt "%(ela))
        trm("./tune")
        with open('steady_slope') as s: slope, af, tau = [[float(j[0]), float(j[1]), float(j[2])] for j in [k.split() for k in s.readlines()]][0]
        if(slope >= 0.000001): low = ela
        else: high = ela
    ela = high

    #OUTPUT FILES
    #trm("mv h_steady.txt ela/hinit_%s"%(glid))
    #trm("mv ice.txt ela/ice_%s"%(glid))
    tau_file = open('ela/tau.txt','a+')
    tau_file.write("%s %f %f %.3f %.7f\n"%(glid, tau, ela, beta, C))
    tau_file.close()

        
    
    

