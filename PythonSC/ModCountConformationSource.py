#f:Distance from center post (Bead level)
#f3: COM close post index l, k (Polymer Level)
#f4: rCOM ThetaCOM Rperp2,Rpara2, Therepara,Thetaee, Thetarange (Polymer Level)
from numpy import *
from math import *
from operator import itemgetter
from ToUnitCell2 import ToUnitCell2

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
    Structuredata=loadtxt("P"+str(Nindex)+"step3.fixave.dat",skiprows = 2);
    tdiff=Structuredata[2,0]-Structuredata[1,0];

    t = dC.time();
    f3= open('P'+str(Nindex)+'com.dat', 'w')
    #f = open('MP'+str(Nindex)+'Bdist2.dat', 'w')
    f4= open('PCP'+str(Nindex)+'.dat', 'w')
    #print f
    for j in range(len(t)):
        xcom=0.0;
        zcom=0.0;
        if abs(PEdata[j,Nindex])>0:
            xt,yt,zt=dC.vecs(j*tdiff,"xu","yu","zu")  ##Absolute Unwrapped Value
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
            RperpCOM=0.0
            for i in range(len(xt)): #writeout the distance with post pid
                xtemp=xt[i]-xfPAU[pid];
                ztemp=zt[i]-zfPAU[pid];
                rtemp=math.sqrt(xtemp*xtemp+ztemp*ztemp);
                RperpCOM+=rtemp-R;
            RperpCOM=RperpCOM/len(xt);
            #RperpCOM=math.sqrt((xfPAU[pid]-xcom)*(xfPAU[pid]-xcom)+(zfPAU[pid]-zcom)*(zfPAU[pid]-zcom));
            ThetaparaCOM=math.atan2(xcom-xfPAU[pid],zcom-zfPAU[pid]);#-pi pi
            #val, pid = min((i, pid) for (pid, val) in enumerate(Cr2temp))
            #f3.write('%d   \n' % pid);
            r2temp=0.0;
            drtemp=0.0;
            Rperp2=0.0;
            Rpara2=0.0;
            Thetapara2=0.0;
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
                if (Thetatemp)>Thetamax:
                    Thetamax=Thetatemp;
                if (Thetatemp)<Thetamin:
                    Thetamin=Thetatemp;

                Rperp2+=(drtemp)*(drtemp);
                Thetapara2+=Thetatemp*Thetatemp;
                #f.write('%d  ' % t[j])
                #f.write('%d  ' % i)
                #f.write('%r  ' % Thetatemp)
                #f.write('%r\n  ' % r2temp)
            Rperp2=Rperp2/len(xt);
            Thetaee=math.atan2(xt[0]-xfPAU[pid],zt[0]-zfPAU[pid])-math.atan2(xt[len(xt)-1]-xfPAU[pid],zt[len(xt)-1]-zfPAU[pid]);
            #Should in -pi pi
            if (Thetaee)>math.pi:
                Thetaee=Thetaee-2*math.pi;
            if (Thetaee)<-1*math.pi:
                Thetaee=Thetaee+2*math.pi;
            Thetaee=Thetaee;
            Thetapara2=Thetapara2/len(xt);
            Thetarange=(Thetamax-Thetamin);
            Rpara2=(R+RperpCOM)*(R+RperpCOM)*Thetapara2;
            Show=math.sqrt(Thetapara2);
            f4.write('%d  ' % t[j])
            f4.write('%r  ' % RperpCOM)
            f4.write('%r  ' % ThetaparaCOM)
            f4.write('%r  ' % Rperp2)
            f4.write('%r  ' % Rpara2)
            f4.write('%r  ' % Show)
            f4.write('%r  ' % Thetaee)
            f4.write('%r \n ' % Thetarange)

    #f.close()
    f3.close()
    f4.close()
