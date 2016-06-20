#This code use coordinate in the unit rectangle cell to calculate polymer sub chain statistics
#Load in:
#"../Run/step2.dump.PA.lammpstrj"
#Output:
#Active: P1InteractionStat.dat #ofBead for P1 P2 P3 P4 (Polymerlevel)
#Active: P1SubChainIndexStat.dat Ntrain,Nloop,Ntail,Nloop2,trainCount,loopCount,tailCount,loop2Count (Polymerlevel)
#Efficiency: 1m30.310s for 10 samples 1000 timesteps

#Connect to MATLAB
#To get # of contact posts: Use P1InteractionStat.dat
#To get portion of forward transition and backward transition: Use # of contact, post and the closed lattice index for COM


#Content:
#P1InteractionStat.dat #ofBead for P1 P2 P3 P4 (Polymerlevel)
#P1Interactionlist.dat P1 P2 P3 P4 YesOrNo(Yes>0 No<0) (Beadlevelel)
#P1SubChainIndexStat.dat Ntrain,Nloop,Ntail,Nloop2,trainCount,loopCount,tailCount,loop2Count (Polymerlevel)
#P1SubChainIndex.dat   train=1, loop=2, tail=3, Crossloop=4 (Beadlevelel)

from numpy import *
from math import *
from operator import itemgetter
import sys
#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)
from dump import dump

if len(sys.argv) != 3:
    print 'exacute: python -f PostSubChainAnl.py  d R';
    sys.exit()
d=float(sys.argv[1]);
print 'd=',d
R=float(sys.argv[2]);
print 'R=',R

ascale=2*R+d;
InterCri=(R+3)*(R+3);
#Load in the Post Position
dP=dump("../Run/step2.dump.PA.lammpstrj");
Px,Py,Pz=dP.vecs(0,"x","y","z")

for Nindex in range(1,5,1):
    ##Load in Chain
    dC=dump("P"+str(Nindex)+"step3.postforvtk");
    dC.map(1,"id",2,"type",3,"x",4,"y",5,"z",6,"xu",7,"yu",8,"zu",9,"vx",10,"vy",11,"vz");
    time,box,atoms,bonds,tris,lines = dC.viz(0);
    dC.sort("id")
    t = dC.time();
    t0=t[0];
    tdiff=t[2]-t[1];
    #Ilist= open('P'+str(Nindex)+'Interactionlist.dat', 'w')
    Iliststat= open('NewP'+str(Nindex)+'InteractionStat.dat', 'w') #t, number of contact bead for each post
    Iliststatplus= open('NewP'+str(Nindex)+'InteractionStat.dat', 'w') #t, numbe of contacting post,

    #Ilist2= open('P'+str(Nindex)+'SubChainIndex.dat', 'w')
    Ilist2stat= open('NewP'+str(Nindex)+'SubChainIndexStat.dat', 'w')
    for j in range(len(t)):
        xt,yt,zt=dC.vecs(t0+j*tdiff,"x","y","z")  ##Use coordinate In rectangle unit cell
        interlist = [[0 for interlistindex in range(len(Px))] for interlistindex in range(len(xt))]  ##[Bead][PostNumber] = 1 contect ,0 no contact.
        interlist2 =zeros(len(xt)); ## [Time] = index of contact post , -1 if no contact
        for k in range(0,len(xt)):
            for i in range(0,len(Px)):
                dx=(Px[i]-xt[k]);
                dz=(Pz[i]-zt[k]);
                if (dx>0.5*box[3]):
                    dx=dx-box[3];
                if (dx<=-0.5*box[3]):
                    dx=dx+box[3];
                if (dz>0.5*box[5]):
                    dz=dz-box[5];
                if (dz<=-0.5*box[5]):
                    dz=dz+box[5];
                r2temp=dx*dx+dz*dz;
                if (r2temp-InterCri)<0:
                    interlist[k][i]=(1);
                    interlist2[k]=(i);
                    break
                else:
                    interlist[k][i]=(0);
                    interlist2[k]=(-1);#
            #Ilist.write('%d  ' % t[j])
            #Ilist.write('%d  ' % k)
            #    checkbead=0;
            #for i in range(0,len(Px)):
            #    Ilist.write('%d ' % interlist[k][i]);
            #Ilist.write('%d ' % interlist2[k]);
            #        checkbead+=interlist[k][i]
            #    if checkbead>1:
            #        print("something wrong 1 bead near 2 posts");
            #Ilist.write('\n')

        interliststat = zeros(len(Px));
        for k in range(0,len(xt)):
            for i in range(0,len(Px)):
                interliststat[i]+=interlist[k][i];
        Iliststat.write('%d  ' % t[j])
        checkpolymer=0;
        for i in range(0,len(Px)):
            Iliststat.write('%r '%interliststat[i]);
            checkpolymer+=interliststat[i];
        Iliststat.write('\n')
        if checkpolymer>len(xt):
            print("something wrong 1 polymer has excess beads %r " % checkpolymer);


        #For absorb case Mark the Beads
        #train=1, loop=2, tail=3, Crossloop=4
        #Mark the type for every bead:
        if checkpolymer>0:
            interlistSubChain = zeros(len(xt));
            #1. Categorize 1 vs 2,3,4
            for k in range(0,len(xt)):
                checkbead=0;
                for i in range(0,len(Px)):
                    checkbead+=interlist[k][i];
                if checkbead>0: #train
                    interlistSubChain[k]=1;
                else:
                    interlistSubChain[k]=2; #Assign to loop first

            #2. Identify 3 from 2,3,4
            #For head
            if interlistSubChain[0]==2:
                ii=0;
                stopcri=0;
                while stopcri!=1:
                    interlistSubChain[ii]=3;
                    ii+=1;
                    if interlistSubChain[ii]==1:
                        stopcri=1;
            #For tail
            if interlistSubChain[len(xt)-1]==2:
                ii=len(xt)-1;
                stopcri=0;
                while stopcri!=1:
                    interlistSubChain[ii]=3;
                    ii-=1;
                    if interlistSubChain[ii]==1:
                        stopcri=1;

            #3. Identify 4 from 2,4
            #For body (No 0 & N-1)
            loopstart=[]
            loopend=[]
            for k in range(1,len(xt)-1):
                if (interlistSubChain[k-1]==1) & (interlistSubChain[k]==2):
                    loopstart.append(k);
                if (interlistSubChain[k]==2) & (interlistSubChain[k+1]==1):
                    loopend.append(k);

            for sl in range(0,len(loopstart)):
                if interlist2[loopstart[sl]-1] != interlist2[loopend[sl]+1]:
                    for kii in range(loopstart[sl],loopend[sl]+1):
                        interlistSubChain[kii]=4;

             ##tempP1=interlist2[k-1];
             ##SL0=k;
             ##ii=0;
             ##while interlistSubChain[k+ii]==2:
             ##    ii+=1;
             ##SLE=k+ii-1;
             ##tempP2=interlist2[k+ii];
             ##if tempP1 != tempP2:  #Change loop type
             ##    for kii in range(SL0,SLE+1):
             ##        interlistSubChain[kii]=4;


            #4. Write Category on file
            #for k in range(0,len(xt)):
            #    Ilist2.write('%d  ' % t[j])
            #    Ilist2.write('%d  ' % k)
            #    Ilist2.write('%d  \n' % interlistSubChain[k])

            ##Structure Statistics
            Ntrain=0;
            Nloop=0;
            Nloop2=0;
            Ntail=0;
            trainCount=0;
            loopCount=0;
            loop2Count=0;
            #Checkloop2Count=0;
            #CheckloopCount=0;
            tailCount=0;
            #head (Cannot be loops)
            if interlistSubChain[0]==1:
                Ntrain+=1;
                trainCount+=1;
            elif  interlistSubChain[0]==3:
                Ntail+=1;
                tailCount+=1;
            #body
            for k in range(1,len(xt)):
                if interlistSubChain[k]==1:
                    trainCount+=1;
                    if interlistSubChain[k-1]!=1:
                        Ntrain+=1;
                elif  interlistSubChain[k]==2:
                    loopCount+=1;
                    if interlistSubChain[k-1]!=2:
                        Nloop+=1;
                elif  interlistSubChain[k]==3:
                    tailCount+=1;
                    if interlistSubChain[k-1]!=3:
                        Ntail+=1;
                elif  interlistSubChain[k]==4:
                    loop2Count+=1;
                    if interlistSubChain[k-1]!=4:
                        Nloop2+=1;
            #loopCount=len(xt)-tailCount-trainCount-loop2Count;
            #if loopCount!=CheckloopCount:
            #    print("something wrong CheckloopCount ! loopCount")

            Ilist2stat.write('%d  ' % t[j])
            Ilist2stat.write('%d %d %d %d %d %d %d %d\n' % (Ntrain,Nloop,Ntail,Nloop2,trainCount,loopCount,tailCount,loop2Count))

    #Ilist.close()
    #Ilist2.close()
    Iliststat.close()
    Ilist2stat.close()
