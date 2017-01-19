import numpy as np
import matplotlib.pyplot as plt

T_step=1.

Beta=0.007
b_0=3900.
s=0.1
v=10.
E=2900.
alpha_m=3.0
length=20000.
g=9.81
rho=917.


def Bs(length):
    return (-1.0/2.0*Beta*s*length**2+Beta*((alpha_m/(1+v*s))*length**0.5+b_0-E)*length)

L=np.array([])
bed=np.array([])
for i in range(int(2000/T_step)):
    length=length+(2.0*(1.0+v*s)*(Bs(length))/(3.0*alpha_m))*length**(-0.5)*T_step
    L=np.append(L,length)
    if int(i*T_step)<500:
        E=2900.
    if 500<=int(i*T_step)<1000:
        E=2800.
    if 1000<=int(i*T_step)<1500:
        E=2900.
    if int(i*T_step)>=1500:
        E=3000.
    bed=np.append(bed,Bs(length)/length)

plt.plot(L)
plt.show()

plt.plot(bed)
plt.show()
length1=30000.
def tau(x):
    return (200.0/1.15)*((1.0-(abs(x/length1-0.5)/0.5)**3.0)+0.15)

def b(x):
    return 200.0*np.exp(-(x/8000.0)**2)

#def H(x):
#    return tau(x)/(rho*g)*(1.0/(h*(x+1)-h)
x=0
h=500.
H=5.
Height=np.array([])
bed=np.array([])
ice=np.array([])
for i in range(int(length1/T_step)):
   if (-b(x)-H+h)/T_step>0.01:
       f=-b(x)-H+h
   else:
       f=0.01
   H=tau(x)/(rho*g)*(1.0/(f))
   h=b(x)+H 
   Height=np.append(Height,H)
   bed=np.append(bed,b(x))
   ice=np.append(ice,h)
   x=x+1

   
   

plt.plot(ice)
plt.plot(bed)
plt.show()
#gufiasfa
#git pull origin master -m "hallo"
#git commit -m "hallo"
#git push origin master -m "hallo"

