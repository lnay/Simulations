from numpy import matrix as mat
from math import sin as sin
from math import cos as cos
import pygame
import sys
import numpy as np
from time import time
pygame.init()

###############adjustable#area##############################

nB=7*1.72*10**-48
k=9.00*(10**9)
pi=3.1415926
e=1.6*(10**-19)
g=9.81*10**-7
proton_mass=1.67*10**-27
#ratio
ratio=10
scale=[ratio*10**-12]
#lengths in pm
sigma=300
SIGMA=[sigma/ratio]
l1=100
l2=70
Oxygen_radius=126
Hydrogen_radius=70
space=300
spacing=space/ratio

#charges in coulombs
Q1=[0.24*e]
Q2=[-0.24*e]
#angles in radians
Theta=104.5*pi/180
Alpha=109.5*pi/180
#window
Box_width=200
Box_height=500
Box_depth=300
Tab_height=100
Time_Interval=[(10**-24)/ratio]
################adjustable#area#############################

#position vectors of charges in H2O
V1=mat([[l1*cos(Theta/2)],[l1*sin(Theta/2)],[0.0]])/ratio
V2=mat([[l1*cos(Theta/2)],[-l1*sin(Theta/2)],[0.0]])/ratio
U1=mat([[-l2*cos(Alpha/2)],[0.0],[l2*sin(Alpha/2)]])/ratio
U2=mat([[-l2*cos(Alpha/2)],[0.0],[-l2*sin(Alpha/2)]])/ratio
#box
size=(Box_width,Box_height+Tab_height)
BOX=[Box_width,Box_height,Box_depth]
size = Box_width, Box_height+100
screen = pygame.display.set_mode(size)
pygame.display.set_caption("H2O Simunlation")

#masses
HM=[proton_mass]
OM=[16*proton_mass]
WM=[18*proton_mass]
#display
HR=Hydrogen_radius/ratio
OR=Oxygen_radius/ratio

red=0,255,0
grey=100,100,100
black=0,0,0
dark=20,20,20
#lists
H2O=[]
Atoms=[]
Particles=[]


###parent class
class particle:
    f=mat([[0.0],[0.0],[0.0]])
    Box=BOX
    def __init__(self,pos,rel,c):
        self.POS=pos
        self.rel=rel
        self.c=c

########

###classes
class lonepair(particle):
    m=[0]
    q=Q2
    
class hydrogen(particle):
    r=HR
    q=Q1
    m=HM
    def display(self):
        a=255*self.c[0]
        colour=(a,a,a)
        pos=self.POS+self.rel
        x=int((pos[0]-self.Box[0]/2)*self.c[0]+self.Box[0]/2)
        y=int((pos[1]-self.Box[1]/2)*self.c[0]+self.Box[1]/2)
        r=int(self.c[0]*self.r)
        pygame.draw.circle(screen,colour,(x,y),r)
        pygame.draw.circle(screen,(0,0,0),(x,y),r,1)

        
class oxygen(particle):
    r=OR
    m=OM
    def display(self):
        a=255*self.c[0]
        colour=(a,0,0)
        x=int((self.POS[0]-self.Box[0]/2)*self.c[0]+self.Box[0]/2)
        y=int((self.POS[1]-self.Box[1]/2)*self.c[0]+self.Box[1]/2)
        r=int(self.c[0]*self.r)
        pygame.draw.circle(screen,colour,(x,y),r)
        pygame.draw.circle(screen,(0,0,0),(x,y),r,1)
####
class water:
    Box=BOX
    t=Time_Interval
    r=SIGMA
    m=WM
    ang_v=mat([[0.0],[0.0],[0.0]])
    v=mat([[0.0],[0.0],[0.0]])
    def __init__(self,_x,_y,_z,Atoms):
        self.c=[0]
        self.pos=mat([[_x],[_y],[_z]])
        self.H1=hydrogen(self.pos,V1,self.c)
        self.H2=hydrogen(self.pos,V2,self.c)
        self.L1=lonepair(self.pos,U1,self.c)
        self.L2=lonepair(self.pos,U2,self.c)
        self.O=oxygen(self.pos,mat([[0],[0],[0]]),self.c)
        Atoms.extend([self.H1,self.H2,self.O])
        self.charges=[self.H1,self.H2,self.L1,self.L2]
        self.particles=[self.H1,self.H2,self.L1,self.L2,self.O]
        self.clouds=[self.O]
        Particles.extend([self.H1,self.H2,self.L1,self.L2,self.O])
    def spin(self):
        Torque_Inertia(self)
        rotate(self)

            
        

        
    
###

def Setup(H2O,Atoms):
    height=BOX[1]
    for x in range(5):
        for y in range(3):
            for z in range(3):
                x = float(x)
                y = float(y)
                z = float(z)
                H2O.append(water((x+1)*spacing,height/2-(y+1)*spacing,(z+1)*spacing,Atoms))

def ClearGlobal():
   pi=e=proton_mass=ratio=sigma=SIGMA=l1=l2=Oxygen_radius=Hydrogen_radiusspace=spacing=Q1=Q2=Theta=Alpha=Box_width=Box_height=Box_depth=Tab_height=Time_Interval=V1=V2=U1=U2=size=BOX=size=HM=OM=WM=HR=OR=None

def resetc(H2O):
    for body in H2O:
        body.c[0]=1-body.pos[2]/(2*body.Box[2])

def update(H2O):
    for body in H2O:
        forces(body)
        bounce(body)
        body.spin()
        
def forces(body):
    force=mat([[0.0],[0.0],[0.0]])  
    for i in body.particles:
        force+=i.f
    a=mat([[0],[g],[0]])+force/body.m[0]
    body.v+=a*body.t[0]
    body.pos+=body.v*body.t[0]

def Torque_Inertia(body):
    Torque=mat([[0.0],[0.0],[0.0]])
    a=b=c=d=e=f=0
    for i in body.charges:
        x=i.rel[0,0]
        y=i.rel[1,0]
        z=i.rel[2,0]
        x2=x**2
        y2=y**2
        z2=y**2
        Torque[0,0]+=i.rel[1,0]*i.f[2,0]-i.rel[2,0]*i.f[1,0]
        Torque[1,0]+=i.rel[2,0]*i.f[0,0]-i.rel[0,0]*i.f[2,0]
        Torque[2,0]+=i.rel[0,0]*i.f[1,0]-i.rel[1,0]*i.f[0,0]
        a+=(y2+z2)*i.m[0]
        b+=-x*y*i.m[0]
        c+=-x*z*i.m[0]
        d+=(x2+z2)*i.m[0]
        e+=-y*z*i.m[0]
        f+=x2+y2
    Inertia=mat([[a,b,c],[b,d,e],[c,e,f]])
    body.ang_v+=np.linalg.solve(Inertia,Torque)*body.t
    
def rotate(body):
    angle=(body.ang_v[0,0]**2+body.ang_v[1,0]**2+body.ang_v[2,0]**2)**0.5
    W=body.ang_v/angle
    angle=angle*body.t[0]
    a=W[0,0]
    b=W[1,0]
    c=W[2,0]
    d=(W[1,0]**2+W[2,0]**2)**0.5
    cod=c/d
    bod=b/d
    T=mat([[1,0,0],[0,cod,-bod],[0,bod,cod]])
    Ti=mat([[1,0,0],[0,cod,bod],[0,-bod,cod]])
    U=mat([[d,0,-a],[0,1,0],[a,0,d]])
    Ui=mat([[d,0,a],[0,1,0],[-a,0,d]])
    co=cos(angle)
    si=sin(angle)
    R=mat([[co,si,0],[-si,co,0],[0,0,1]])
    M=Ti*Ui*R*U*T
    for i in body.charges:
        i.rel=M*i.rel
        
def Screen(Atoms):
    Atoms.sort(key=lambda x: (x.POS[2,0]+x.rel[2,0]),reverse=True)
    resetc(H2O)
    for i in Atoms:
        i.display()

def bounce(body):
    for i in range(3):
        if body.pos[i,0]-body.r[0]<0:
            body.pos[i,0]=body.r[0]
            body.v[i,0]=-body.v[i,0]
        if body.pos[i,0]+body.r[0]>body.Box[i]:
            body.pos[i,0]=body.Box[i]-body.r[0]
            body.v[i,0]=-body.v[i,0]
            

    
def Interact(H2O):
    scl=scale[0]
    for i in Particles:
        i.f[0,0]=0.0
        i.f[1,0]=0.0
        i.f[2,0]=0.0
    for g in range(len(H2O)-1):
        for h in range(g+1,len(H2O)):
            for i in H2O[g].clouds:
                for j in H2O[h].clouds:
                    CloudInteract(i,j,scl)
            for i in H2O[g].charges:
                for j in H2O[h].charges:
                    ChargesInteract(i,j,scl)

def ChargesInteract(i,j,scl):
    d=i.POS+i.rel-j.POS-j.rel
    d=d*scl
    r2=d[0,0]**2+d[1,0]**2+d[2,0]**2
    if r2==0:
        d[0,0]=1
        r2=1
    electro=i.q[0]*j.q[0]*k/r2
    r=r2**0.5
    tot=electro/r*10**33
    f=d*tot
    i.f=mat(i.f+f)
    j.f=mat(j.f-f)


def CloudInteract(i,j,scl):
    d=i.POS+i.rel-j.POS-j.rel
    d=d*scl
    r2=d[0,0]**2+d[1,0]**2+d[2,0]**2
    repel=nB/r2              ###just for H2O
    r=r2**0.5
    tot=repel/r
    f=d*tot
    i.f=mat(i.f+f)
    j.f=mat(j.f-f)


    
def Background():
    screen.fill(grey)

    rect=pygame.Rect(Box_width/4, Box_height/4, Box_width/2, Box_height/2)
    pygame.draw.rect(screen,dark,rect)

    pygame.draw.line(screen,black,(0,0),(Box_width/4, Box_height/4))
    pygame.draw.line(screen,black,(3*Box_width/4, 3*Box_height/4),(Box_width, Box_height))
    pygame.draw.line(screen,black,(Box_width/4,3*Box_height/4),(0, Box_height))
    pygame.draw.line(screen,black,(3*Box_width/4,Box_height/4),(Box_width,0))
        
def IncrementSpeed(n,H2O):
    for water in H2O:
        water.v = water.v*n
        water.ang_v = water.ang_v*n
        
if __name__=="__main__":
    Setup(H2O,Atoms)
    ClearGlobal()
    Continue=True
    while Continue:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                break
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_UP:
                    IncrementSpeed(1.25,H2O)
                if event.key == pygame.K_DOWN:
                    IncrementSpeed(0.8,H2O)
                if event.key == pygame.K_RIGHT:
                    IncrementSpeed(0.0,H2O)
        time0=time()
        Background()
        Interact(H2O)
        
        update(H2O)
        Screen(Atoms)
        pygame.display.flip()
        time1=time()
        print(time1-time0)

        #ChargesInteract(Charges,k)
        #print(H2O[0].H1.force)
        



