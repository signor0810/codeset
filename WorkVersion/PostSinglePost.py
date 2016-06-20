#This code use absolute coordinate to calculate polymer conformation near a single post
#Load in:
#step3.PE.dat
#P"+str(Nindex)+"step3.postforvtk
#Output:
#f:Distance from center post (Bead level) 
#Active:P1com.dat f3 COM close post index l, k (Polymer Level)
#Active:PCP1.dat f4 t RCOM ThetaCOM Rperp2,Rpara2, Therepara,Thetaee, Thetarange (Polymer Level)
#t RperpCOM ThetaparaCOM Rperp2 Rpara2 Show Ree Thetaee Rfar Rnear Thetarange (Polymer Level)
#Efficiency: 0m41.7s for 10 samples 1000 time steps


from numpy import *
from math import *
from operator import itemgetter
from ToUnitCell2 import ToUnitCell2
from dump import dump

import sys
#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)
if len(sys.argv) != 3:
    print 'exacute: python PostSinglePost.py  d R';
    sys.exit()
d=float(sys.argv[1]);
print 'd=',d
R=float(sys.argv[2]);
print 'R=',R

ascale=2*R+d;
a1=[ascale,0];
a2=[0.5*ascale,math.sqrt(3)/2*ascale];

for Nindex in range(1,5,1):
    #Load in Chain
    dC=dump("P"+str(Nindex)+"step3.postforvtk");
    dC.map(1,"id",2,"type",3,"x",4,"y",5,"z",6,"xu",7,"yu",8,"zu",9,"vx",10,"vy",11,"vz")

    #Load in C-PA Energy
    PEdata=loadtxt("step3.PE.dat",skiprows = 2);
    #Load in C Property
    t = dC.time();
    t0=t[0];
    tdiff=t[2]-t[1];
    f3= open('P'+str(Nindex)+'com.dat', 'w')
    #f = open('MP'+str(Nindex)+'Bdist2.dat', 'w')
    f4= open('PCP'+str(Nindex)+'.dat', 'w')
    #print f
    for j in range(len(t)):
        xcom=0.0;
        zcom=0.0;
        if abs(PEdata[j,Nindex])>0:
            xt,yt,zt=dC.vecs(t0+j*tdiff,"xu","yu","zu")  ##Absolute Unwrapped Value
            for i in range(len(xt)): #Count com
                xcom+=xt[i];
                zcom+=zt[i]

            xcom/=len(xt);
            zcom/=len(xt);
            #Get The index of the unit cell
            indexUnit=ToUnitCell2(xcom,zcom,R,d);
            #Determine which post is the nearest
            f3.write('%d  ' % t[j])
            f3.write('%r  ' % xcom)
            f3.write('%r  ' % zcom)
            f3.write('%d  ' % indexUnit.k)
            f3.write('%d  \n' % indexUnit.l)

            xPAU=[0,ascale,0.5*ascale,1.5*ascale]
            zPAU=[0,0,ascale*math.sqrt(3)/2,ascale*math.sqrt(3)/2]
            xfPAU=[0,ascale,0.5*ascale,1.5*ascale]
            zfPAU=[0,0,ascale*math.sqrt(3)/2,ascale*math.sqrt(3)/2]

            for i in range(0,len(xPAU)):
                xfPAU[i]=xPAU[i]+indexUnit.k*a1[0]+indexUnit.l*a2[0];
                zfPAU[i]=zPAU[i]+indexUnit.k*a1[1]+indexUnit.l*a2[1];

            ##f3.write('%r ' % xfPAU[0]);
            ##f3.write('%r ' % zfPAU[0]);
            ##f3.write('%r ' % xfPAU[1]);
            ##f3.write('%r ' % zfPAU[1]);
            ##f3.write('%r ' % xfPAU[2]);
            ##f3.write('%r ' % zfPAU[2]);
            ##f3.write('%r ' % xfPAU[3]);
            ##f3.write('%r ' % zfPAU[3]);


            ##Find the center post
            Cr2temp=[];
            for i in range(0,len(xPAU)):
                Cr2temp.append((xfPAU[i]-xcom)*(xfPAU[i]-xcom)+(zfPAU[i]-zcom)*(zfPAU[i]-zcom));
            pid=Cr2temp.index(min(Cr2temp))


            ##COM in post coordinate
            #Determinate Radical Center
            RperpCOM=0.0;
            ThetaparaCOM=0.0;
            for i in range(len(xt)): #writeout the distance with post pid
                xtemp=xt[i]-xfPAU[pid];
                ztemp=zt[i]-zfPAU[pid];
                rtemp=math.sqrt(xtemp*xtemp+ztemp*ztemp);
                RperpCOM+=rtemp-R;
                ThetaparaCOM+=math.atan2(xt[i]-xfPAU[pid],zt[i]-zfPAU[pid]);#-pi pi
            RperpCOM=RperpCOM/len(xt);
            ThetaparaCOM=ThetaparaCOM/len(xt);#-pi pi

            #RperpCOM=math.sqrt((xfPAU[pid]-xcom)*(xfPAU[pid]-xcom)+(zfPAU[pid]-zcom)*(zfPAU[pid]-zcom));
            #val, pid = min((i, pid) for (pid, val) in enumerate(Cr2temp))
            #f3.write('%d   \n' % pid);

            r2temp=0.0;
            drtemp=0.0;

            ##Gyration Radius:
            Rperp2=0.0;
            Thetapara2=0.0;

            #Stretch;
            Rfar2=0.0;
            Rnear2=100000;
            Rpara2=0.0;
            ##pata
            Thetamin=2*math.pi;
            Thetamax=-2*math.pi;

            for i in range(len(xt)): #writeout the distance with post pid
                xtemp=xt[i]-xfPAU[pid];
                ztemp=zt[i]-zfPAU[pid];
                r2temp=xtemp*xtemp+ztemp*ztemp;
                drtemp=(math.sqrt(r2temp)-R-RperpCOM);
                Thetatemp=math.atan2(xtemp,ztemp)-ThetaparaCOM;#Should in -pi pi
                if (Thetatemp)>math.pi:
                    Thetatemp=Thetatemp-2*math.pi;
                if (Thetatemp)<-1*math.pi:
                    Thetatemp=Thetatemp+2*math.pi;
                #f.write('%d  ' % t[j])
                #f.write('%d  ' % i)
                #f.write('%r  ' % Thetatemp)
                #f.write('%r\n  ' % r2temp)

                #Stretch
                if r2temp>Rfar2:
                    Rfar2=r2temp;
                if r2temp<Rnear2:
                    Rnear2=r2temp;
                if (Thetatemp)>Thetamax:
                    Thetamax=Thetatemp;
                if (Thetatemp)<Thetamin:
                    Thetamin=Thetatemp;
                #Gyration Radius
                Rperp2+=(drtemp)*(drtemp);
                Thetapara2+=Thetatemp*Thetatemp;
            ##Gyration Radius
            Rperp2=Rperp2/len(xt);
            Thetapara2=Thetapara2/len(xt);
            Show=math.sqrt(Thetapara2);
            Rpara2=(R+RperpCOM)*(R+RperpCOM)*Thetapara2;

            ##End to end
            Ree=math.sqrt((xt[0]-xfPAU[pid])*(xt[0]-xfPAU[pid])+(zt[0]-zfPAU[pid])*(zt[0]-zfPAU[pid]))-math.sqrt((xt[len(xt)-1]-xfPAU[pid])*(xt[len(xt)-1]-xfPAU[pid])+(zt[len(xt)-1]-zfPAU[pid])*(zt[len(xt)-1]-zfPAU[pid]));
            Thetaee=math.atan2(xt[0]-xfPAU[pid],zt[0]-zfPAU[pid])-math.atan2(xt[len(xt)-1]-xfPAU[pid],zt[len(xt)-1]-zfPAU[pid]);
            if (Thetaee)>math.pi:
                Thetaee=Thetaee-2*math.pi;
            if (Thetaee)<-1*math.pi:
                Thetaee=Thetaee+2*math.pi;
            Thetaee=Thetaee;#Should in -pi pi

            ##Stretch
            Rfar=math.sqrt(Rfar2)-R;
            Rnear=math.sqrt(Rnear2)-R;
            Thetarange=(Thetamax-Thetamin);


            f4.write('%d  ' % t[j])
            f4.write('%r  ' % RperpCOM)
            f4.write('%r  ' % ThetaparaCOM)
            f4.write('%r  ' % Rperp2)
            f4.write('%r  ' % Rpara2)
            f4.write('%r  ' % Show)
            f4.write('%r  ' % Ree)
            f4.write('%r  ' % Thetaee)
            f4.write('%r  ' % Rfar)
            f4.write('%r  ' % Rnear)
            f4.write('%r \n ' % Thetarange)

    #f.close()
    f3.close()
    f4.close()
