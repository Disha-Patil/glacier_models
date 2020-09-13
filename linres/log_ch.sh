#######################step change

#properties

#read from the file id_v_a_id_beta_e_bt.txt


#linres evolution:
#check formula with the manuscript
awk -F'[ ]' 'BEGIN{gamma=1.286;tmax=500}   { v0=$2; a0=$3; beta=$5; bt=$7;  h=1000*v0/a0;   tau0= -1.0/(bt/(gamma*h)+beta); alpha0= beta*tau0*50./(gamma*h); taua=2.56*tau0; Dv=-1.71*alpha0*v0;   tauv=taua*0.687; Da=a0*Dv/(1.93*v0);    vtot+=v0;atot+=a0;       for (t=1;t<tmax;t++){   De=50;       for (t1=t;t1<tmax;t1++) {v[t1]+=(Dv*De/(50*tauv))*exp((t-t1)/tauv); a[t1]+=(Da*De/(50*taua))*exp((t-t1)/taua);}  }         } END{v[0]=0;a[0]=0;for(t=0;t<tmax;t++) print t,vtot+v[t],atot+a[t];}' id_v_a_id_beta_e_bt.txt >./tva_linres_step.txt

#scaling evolution:

awk 'BEGIN{gamma=1.286;tmax=500; fnum=0;opath="./scaling_results/";}\
FILENAME=="id_v_a_id_beta_e_bt.txt" {v0[$1]=$2*1e9; a0[$1]=$3*1e6; beta[$1]=$5; e[$1]=$6; bt[$1]=$7;flag[$1]=1;}\
NR>FNR {if(FNR==1)fnum++; id[fnum]=substr(FILENAME,19); z[fnum,FNR]=$1; hyp[fnum,FNR]=$2*1e4; nband[fnum]=FNR;}\
END{\
for(i=1;i<=fnum;i++){ a=0;dv0=0;id0=id[i]; v=v0[id0];outf=opath"va_"id0"_step";\
   print id[i];
   print "after"
   for(j=1;j<=nband[i];j++){a+=hyp[i,j]; b=beta[id0]*(z[i,j]-e[id0]); dv0+=(b>1?1:b)*hyp[i,j];};hyp[i,1]+=(a0[id0]-a);  
   print v;
   #print outf,a,v, dv0,beta[id0],z[i,1],hyp[i,1],b;
   print 0,v*1e-9,a0[id0]*1e-6 >> outf;\
   for(t=1;t<tmax;t++){\
      de=50;a=0;dv=0;\
      for(j=1;j<=nband[i];j++){a+=hyp[i,j]; b=beta[id0]*(z[i,j]-e[id0]-de); dv+=(b>1?1:b)*hyp[i,j];}\
      dv-=dv0;da=a*dv/(gamma*v);v+=dv;\
      for(j=1;j<=nband[i];j++){if(hyp[i,j]==0)continue; hyp[i,j]+=da;\
         if(hyp[i,j]<0){da=-hyp[i,j];hyp[i,j]=0;}\
         else break;  };\
      if(v<0)v=0;print t,v*1e-9,a*1e-6 >> outf;\
    };\
};\
}' id_v_a_id_beta_e_bt.txt ./hypso_arr/hypso*


#sia total v a
awk -F'[ ]' 'NR==FNR{sel[$1]=1;} NR>FNR{split(FILENAME,f,"_"); if (sel[f[4]]+0){a[int($1)]+=$3;v[int($1)]+=$2;}}END{for (t=0;t<1000;t+=1) print t,v[t],a[t]}' id_v_a_id_beta_e_bt.txt ./output_change_ela/va* >tva_sia_step.txt

awk -F'[ ]' '{v[int($1)]+=$2; a[int($1)]+=$3;}END{for (t=0;t<1000;t+=1) print t,v[t],a[t]}' ./scaling_results/va*step >tva_scaling_step.txt

##plots
#v vs t
#se xti 100; se mxti 2; se yti 50; se myti 5;se yla "volume (km^3)"; se xla "time (yr)";p [:500] 'tva_linres_step.txt' u 1:2 w l lc 6 lw 2 ti "linres", 'tva_scaling_step.txt' u 1:2 w l lc 7 lw 2 ti "scaling", "tva_sia_step.txt" u 1:2 w l lc 5 lw 2 ti "SIA"

#se xti 100; se mxti 2; se yti 200; se myti 2;se yla "area (km^2)"; se xla "time (yr)";p [:500] 'tva_linres_step.txt' u 1:3 w l lc 6 lw 2 ti "linres",'tva_scaling_step.txt' u 1:3 w l lc 7 lw 2 ti "scaling", 'tva_sia_step.txt' u 1:3 w l lc 5 lw 2 ti "SIA"


#noisy linres
#for the band plot
for i in {1..30}
do
awk -v seed="$RANDOM" 'BEGIN{FS=OFS=" ";srand(seed); gamma=1.286;tmax=500;c1=2.56+2*0.04*(2*rand()-1);c2=1.71+2*0.03*(2*rand()-1);c3=0.687+2*0.004*(2*rand()-1);c4=1.93+2*0.02*(2*rand()-1);}   { v0=$2; a0=$3; beta=$5; bt=$7;  h=1000*v0/a0;   tau0= -1.0/(bt/(gamma*h)+beta); alpha0= beta*tau0*50./(gamma*h); taua=c1*tau0; Dv=-c2*alpha0*v0;   tauv=c3*taua; Da=a0*Dv/(c4*v0);    vtot+=v0;atot+=a0;       for (t=1;t<tmax;t++){   De=50;       for (t1=t;t1<tmax;t1++) {v[t1]+=(Dv*De/(50*tauv))*exp((t-t1)/tauv); a[t1]+=(Da*De/(50*taua))*exp((t-t1)/taua);}  }         } END{v[0]=0;a[0]=0;for(t=0;t<tmax;t++) print t,vtot+v[t],atot+a[t]; printf("\n")}' id_v_a_id_beta_e_bt.txt >tva_linres_step_noise.txt
done


