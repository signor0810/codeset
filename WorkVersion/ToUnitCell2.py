import math
def ToUnitCell2(X,Z,R,d):
    ##This function return the index of the unit cell (k,l)
    ascale=(2*R+d);
    a1=[ascale,0];
    a2=[0.5*ascale,math.sqrt(3)/2*ascale];
    l=math.floor(Z/a2[1]);
    k=math.floor((X-l*a2[0])/a1[0]);
    ##P1=[0,0];
    ##P2=[0,ascale];
    ##P3=[0.5*ascale,math.sqrt(3)/2*ascale];
    ##P4=[ascale,math.sqrt(3)/2*ascale];
    ##Pf1=[0,0];
    ##Pf2=[0,ascale];
    ##Pf3=[0.5*ascale,math.sqrt(3)/2*ascale];
    ##Pf4=[ascale,math.sqrt(3)/2*ascale];

    ##for i in range(0,1):
    ##    Pf1[i]=P1[i]+k*a1[i]+l*a2[i];
    ##    Pf2=P2[i]+k*a1[i]+l*a2[i];
    ##    Pf3=P3[i]+k*a1[i]+l*a2[i];
    ##    Pf4=P4[i]+k*a1[i]+l*a2[i];
    return ReturnValue(k,l)

class ReturnValue(object):
    def __init__(self,k,l):
        self.k = k
        self.l = l

