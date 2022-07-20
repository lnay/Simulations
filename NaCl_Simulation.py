import sys, pygame

pygame.init()
from time import time


# testing pycharm git integration
Box_width=200
Box_height=500
Box_depth=300
size = Box_width, Box_height+100
DISPLAY_SET_MODE = pygame.display.set_mode(size)
screen = DISPLAY_SET_MODE
pygame.display.set_caption("Ionic simulations")

grey=100,100,100
black=0,0,0
dark=20,20,20

ratio=10
g=10**-8
G=6.67*(10**-11)

Actual_Na_Radius=120
Actual_Cl_Radius=170
Na_Radius=Actual_Na_Radius/ratio
Cl_Radius=Actual_Cl_Radius/ratio

Actual_Bond_length=240
Bond_length=Actual_Bond_length/ratio

n_Cl=9
n_Na=7
n=8
k=9.00*(10**9)
e=1.60*(10**(-19))


qNa_X_qCl=-e**2
qNa_X_qNa=e**2
qCl_X_qCl=e**2
qNa_X_qCl_k=qNa_X_qCl*k
qNa_X_qNa_k=qNa_X_qNa*k
qCl_X_qCl_k=qCl_X_qCl*k
B=qNa_X_qNa*k*(Bond_length**7)/8
print(B)
Bond_length=Bond_length*1.1
mass_Na=23*1.67 * 10**-27
mass_Cl=35.5*1.67 * 10**-27
mass_Cl_X_Na=mass_Na*mass_Cl
mass_Cl_X_Cl=mass_Cl*mass_Cl
mass_Na_X_Na=mass_Na*mass_Na




Cl_list=[]
Na_list=[]
LIST=[]

Frame_Interval=1000/ratio

class ion():
    def __init__(self,_x,_y,_z,_name):
        self.name=_name
        self.x=_x
        self.y=_y
        self.z=_z
        
        self.fx=0
        self.fy=0
        self.fz=0
        
        self.ax=0
        self.ay=0
        self.az=0
        
        self.vx=0
        self.vy=0
        self.vz=0

def ResetForcesOnIons(Cl_list,Na_list,g,mass_Na,mass_Cl):
    for ion in Na_list:
        ion.fx=0
        ion.fz=0
        ion.fy=9.81*mass_Na*g
    for ion in Cl_list:
        ion.fx=0
        ion.fz=0
        ion.fy=9.81*mass_Cl*g

    
def IonicInteraction1(cation,anion,n,B,G,mass_Cl_X_Na,qNa_X_qCl_k):
    dx=-(anion.x-cation.x)
    dy=-(anion.y-cation.y)
    dz=-(anion.z-cation.z)
    r2=dx**2+dy**2+dz**2
    if r2==0:
        dx=1
        r2=2
    r=r2**0.5
    
    Electrostatic=qNa_X_qCl_k/r2
    repell=n*B*(r**(-n-1))
    total=(Electrostatic+repell)/(r)
    cation.fx+=dx*total
    anion.fx+=-dx*total
    cation.fy+=dy*total
    anion.fy+=-dy*total
    cation.fz+=dz*total
    anion.fz+=-dz*total



def IonicInteraction2(cation1,cation2,n_Na,qNa_X_qNa_k,B,mass_Na_X_Na,G):
    dx=-(cation2.x-cation1.x)
    dy=-(cation2.y-cation1.y)
    dz=-(cation2.z-cation1.z)
    r2=dx**2+dy**2+dz**2
    if r2==0:
        dx=1
        r2=2
    r=r2**0.5
    Electrostatic=qNa_X_qNa_k/r2
    total=(Electrostatic)/(r)
    cation1.fx+=dx*total
    cation2.fx+=-dx*total
    cation1.fy+=dy*total
    cation2.fy+=-dy*total
    cation1.fz+=dz*total
    cation2.fz+=-dz*total

def IonicInteraction3(anion1,anion2,n_Cl,qCl_X_qCl_k,B,mass_Cl_X_Cl,G):
    dx=-(anion2.x-anion1.x)
    dy=-(anion2.y-anion1.y)
    dz=-(anion2.z-anion1.z)
    r2=dx**2+dy**2+dz**2
    if r2==0:
        dx=1
        r2=2
    Electrostatic=qCl_X_qCl_k/r2
    r=r2**0.5
    total=(Electrostatic)/(r)
    anion1.fx+=dx*total
    anion2.fx+=-dx*total
    anion1.fy+=dy*total
    anion2.fy+=-dy*total
    anion1.fz+=dz*total
    anion2.fz+=-dz*total
    


def IonicInteractions(Cl_list,Na_list,n,B,G,mass_Cl_X_Na,mass_Na_X_Na,mass_Cl_X_Cl,qNa_X_qNa_k,qNa_X_qCl_k,qCl_X_qCl_k):
    for Na in Na_list:
        for Cl in Cl_list:
            IonicInteraction1(Na,Cl,n,B,G,mass_Cl_X_Na,qNa_X_qCl_k)
    for i in range(len(Na_list)-1):
        for j in range(i+1,len(Na_list)):
            IonicInteraction2(Na_list[i],Na_list[j],n_Na,qNa_X_qNa_k,B,mass_Na_X_Na,G)
    for i in range(len(Cl_list)-1):
        for j in range(i+1,len(Cl_list)):
            IonicInteraction3(Cl_list[i],Cl_list[j],n_Cl,qCl_X_qCl_k,B,mass_Cl_X_Cl,G)

def SetupIonLists():
    global Na_list
    global Cl_list
    Na_list=[]
    Cl_list=[]
    LastIon="Cl"
    for _z in range(4):
        for _y in range(5):
            for _x in range(5):
                if LastIon=="Cl":
                    x=_x*(Bond_length*1.1)+50
                    y=Box_height-_y*(Bond_length*1.1)
                    z=_z*(Bond_length*1.1)+50
                    Na_list.append(ion(x,y,z,"Na"))
                    LIST.append(Na_list[len(Na_list)-1])
                    LastIon="Na"
                else:
                    x=_x*(Bond_length*1.1)+50
                    y=Box_height-_y*(Bond_length*1.1)
                    z=_z*(Bond_length*1.1)+50
                    Cl_list.append(ion(x,y,z,"Cl"))
                    LIST.append(Cl_list[len(Cl_list)-1])
                    LastIon="Cl"
"""
option7 = Button(window, text = 'Reset', command = SetupIonLists)
option7.pack()
"""
def Display(ion,Cl_Radius,Na_Radius,Box_width,Box_depth,Box_height):
    
    c=1-ion.z/(2*Box_depth)
    x=int((ion.x-Box_width/2)*c+Box_width/2)
    y=int((ion.y-Box_height/2)*c+Box_height/2)
    #x=int((Box_depth*(2*ion.x-Box_width)/(ion.z+Box_depth)+Box_width)/2)
    #y=int((Box_depth*(2*ion.y-Box_height)/(ion.z+Box_depth)+Box_height)/2)
    if ion.name=="Cl":
        colour=(0,255, int(255*ion.z/Box_depth))
        #r=int(Cl_Radius-Box_depth/1020)
        r=int(c*Cl_Radius)
        
    else:
        colour=(255,0, int(255*ion.z/Box_depth))
        #r=int(Na_Radius-Box_depth/1020)
        r=int(c*Na_Radius)
    if ion.z>Box_depth or ion.z<0:
        print("z is",ion.z)
    pygame.draw.circle(screen,colour,(x,y),r)
    pygame.draw.circle(screen,(0,0,0),(x,y),r,1)
    pygame.draw.circle(screen,(255,255,255),(int(x+r/3),int(y-r/3)),int(r/4))

def SortAndDisplay(Cl_list,Na_list):
    
    #Cl_list=Cl_list_original
    #Na_list=Na_list_original
    
    Cl_list2=[]
    Na_list2=[]
    LIST.sort(key=lambda x: x.z, reverse=True)
    for i  in LIST:
        Display(i,Cl_Radius,Na_Radius,Box_width,Box_depth,Box_height)

def IncrementSpeed(n):
   for ion in LIST:
      ion.vx *= n
      ion.vy *= n
      ion.vz *= n
def PauseSpeed():
   for ion in LIST:
      ion.vx = 0
      ion.vy = 0
      ion.vz = 0
def UpdateIons(Cl_list,Na_list,mass_Cl,mass_Na,ratio):
    for Cl in Cl_list:
        Cl.ax=Cl.fx/mass_Cl
        Cl.ay=Cl.fy/mass_Cl
        Cl.az=Cl.fz/mass_Cl
        Cl.vx+=Cl.ax*Frame_Interval
        Cl.vy+=Cl.ay*Frame_Interval
        Cl.vz+=Cl.az*Frame_Interval
        Cl.x+=Cl.vx*Frame_Interval
        Cl.y+=Cl.vy*Frame_Interval
        Cl.z+=Cl.vz*Frame_Interval
        if Cl.x>Box_width: #bounce x
            Cl.x=Box_width
            Cl.vx=-Cl.vx
        elif Cl.x<0:
            Cl.x=0
            Cl.vx=-Cl.vx
        if Cl.y>Box_height: #bounce y
            Cl.y=Box_height
            Cl.vy=-Cl.vy
        elif Cl.y<0:
            Cl.y=0
            Cl.vy=-Cl.vy
        if Cl.z>Box_depth: #bounce z
            Cl.z=Box_depth-Cl_Radius
            Cl.vz=-Cl.vz
            
        elif Cl.z<0:
            Cl.z=Cl_Radius
            Cl.vz=-Cl.vz
           
    for Na in Na_list:
        Na.ax=Na.fx/mass_Na
        Na.ay=Na.fy/mass_Na
        Na.az=Na.fz/mass_Na
        Na.vx+=Na.ax*Frame_Interval
        Na.vy+=Na.ay*Frame_Interval
        Na.vz+=Na.az*Frame_Interval
        Na.x+=Na.vx*Frame_Interval
        Na.y+=Na.vy*Frame_Interval
        Na.z+=Na.vz*Frame_Interval
        if Na.x>Box_width: #bounce x
            Na.x=Box_width
            Na.vx=-Na.vx
        elif Na.x<0:
            Na.x=0
            Na.vx=-Na.vx
        if Na.y>Box_height: #bounce y
            Na.y=Box_height
            Na.vy=-Na.vy
        elif Na.y<0:
            Na.y=0
            Na.vy=-Na.vy
        if Na.z>Box_depth: #bounce z
            Na.z=Box_depth-Na_Radius
            Na.vz=-Na.vz
            #print("bounc in z")
        elif Na.z<0:
            Na.z=Na_Radius
            Na.vz=-Na.vz
            #print("bounc in z")






SetupIonLists()
i=0
while True:
    #print(i)
    #i+=1
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            break
        elif event.type == pygame.KEYDOWN:
            if event.key == pygame.K_UP:
                IncrementSpeed(1.25)
            if event.key == pygame.K_DOWN:
                IncrementSpeed(0.8)
            if event.key == pygame.K_RIGHT:
                PauseSpeed()
    time0=time()
    screen.fill(grey)

    rect=pygame.Rect(Box_width/4, Box_height/4, Box_width/2, Box_height/2)
    pygame.draw.rect(screen,dark,rect)

    pygame.draw.line(screen,black,(0,0),(Box_width/4, Box_height/4))
    pygame.draw.line(screen,black,(3*Box_width/4, 3*Box_height/4),(Box_width, Box_height))
    pygame.draw.line(screen,black,(Box_width/4,3*Box_height/4),(0, Box_height))
    pygame.draw.line(screen,black,(3*Box_width/4,Box_height/4),(Box_width,0))
    
    ResetForcesOnIons(Cl_list,Na_list,g,mass_Na,mass_Cl)
    SortAndDisplay(Cl_list,Na_list)
    #update
    IonicInteractions(Cl_list,Na_list,n,B,G,mass_Cl_X_Na,mass_Na_X_Na,mass_Cl_X_Cl,qNa_X_qNa_k,qNa_X_qCl_k,qCl_X_qCl_k)
    
    UpdateIons(Cl_list,Na_list,mass_Cl,mass_Na,ratio)
    #print(Na_list[0].fx,Na_list[0].ax,Na_list[0].vx,Na_list[0].x)

    pygame.display.flip()
    time1=time()
    print(time1-time0)
    #canvas.after(0)

