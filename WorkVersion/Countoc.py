#Calculate Orientation correlation up to half chain length
#Load in:
#P"+str(Nindex)+"step3.postforvtk
#Output
#f:Distance from center post (Bead level)
#f: P1OrientationCorr.dat
#Active: TAP1OrientationCorr.dat
#Efficiency: 2m48.106s 1 Sample 1000 timsteps
from numpy import *
from math import *
from dump import dump

Nlist=[160,80,40,20]
for Nindex in range(1,5,1):
    N=Nlist[Nindex-1];
    #print N
    #Load in Chain
    dC=dump("P"+str(Nindex)+"step3.postforvtk");
    dC.map(1,"id",2,"type",3,"x",4,"y",5,"z",6,"xu",7,"yu",8,"zu",9,"vx",10,"vy",11,"vz")
    t = dC.time();
    #f= open('P'+str(Nindex)+'OrientationCorr.dat', 'w')
    f2= open('TAP'+str(Nindex)+'OrientationCorr.dat', 'w')
    TC3D=[0]*(int(N/2));
    TC2D=[0]*(int(N/2));
    #print f
    for j in range(len(t)):
        xt,yt,zt=dC.vecs(j*(t[1]-t[0]),"xu","yu","zu")  ##Absolute Unwrapped Value
        C3D=[0]*int(len(xt)/2);
        C2D=[0]*int(len(xt)/2);
        for k in range(0,int(len(xt)/2),1): #k=i-j
            v1=[0]*3;
            v2=[0]*3;
            l=0;
            for i1 in range(len(xt)-1-k):
                i2=i1+k;
                v1[0]=(xt[i1+1]-xt[i1]);
                v2[0]=(xt[i2+1]-xt[i2]);
                v1[1]=(yt[i1+1]-yt[i1]);
                v2[1]=(yt[i2+1]-yt[i2]);
                v1[2]=(zt[i1+1]-zt[i1]);
                v2[2]=(zt[i2+1]-zt[i2]);
                C3D[k]+=((v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])/(math.sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])*math.sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2])));
                C2D[k]+=((v1[0]*v2[0]+v1[1]*v2[1])/(math.sqrt(v1[0]*v1[0]+v1[1]*v1[1])*math.sqrt(v2[0]*v2[0]+v2[1]*v2[1])));
                l+=1;
            C3D[k]=C3D[k]/l;
            C2D[k]=C2D[k]/l;
            TC3D[k]+=C3D[k];
            TC2D[k]+=C2D[k];
            #f.write('%d  ' % j)
            #f.write('%d  ' % k)
            #f.write('%r  ' % C3D[k])
            #f.write('%r\n  ' % C2D[k])
    for k  in range(1,int(len(xt)/2),1): #k=i-j
        TC3D[k]/=len(t);
        TC2D[k]/=len(t);
        f2.write('%d  ' % k)
        f2.write('%r  ' % TC3D[k])
        f2.write('%r\n  ' % TC2D[k])

    f2.close()
    #f.close()
