import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import math as m

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
length1=30000.
tau0 = 150E3

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
def tau(x):
    delta = 0.15
    taum = 150000
    return (taum/(1+delta))*(1.0-(abs(x/(length1/3.)-0.5)/0.5)**3.0+delta)

def b(x):
    return 200.0*np.exp(-(x/8000.0)**2)

hoogtelijst = np.array([1930,1800,1700,1600,1500,1400,1300,1200,1100,1000,900,800,600,500,400])
lengtelijst = np.array([0,2000,4500,5500,6250,6500,6750,7000,7125,7525,7825,8025,8275,8525,8775])
plt.plot(lengtelijst, hoogtelijst)
plt.show()
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
   
def H(lengtearray,hoogtearray):
    Harray = np.array([])
    hoogte = 0.
    Harray = np.append(Harray, tau(lengtearray[0])/(rho*g)*(1.0/((-(hoogtearray[1]-hoogtearray[0]))/(lengtearray[1]-lengtearray[0]))))
    for i in range(0, len(lengtearray)-1):
        hoogte = tau(hoogtearray[i])/(rho*g)*(1.0/((-(hoogtearray[i+1]-hoogtearray[i]))/(lengtearray[i+1]-lengtearray[i])))
        Harray = np.append(Harray,hoogte)
    return Harray
    
def bed(lengtearray, hoogtearray):
    bedarray = np.array([])
    bed = 0.
    for i in range(0, len(lengtearray)):
        bed = hoogtearray[i]- H(lengtearray,hoogtearray)[i]
        bedarray = np.append(bedarray, bed)
    return bedarray
    
#print(bed(lengtelijst,hoogtelijst))    
print(H(lengtelijst, hoogtelijst))
plt.plot(lengtelijst, bed(lengtelijst,hoogtelijst))
plt.plot(lengtelijst, hoogtelijst, color = "green")
plt.show()

def Gaussian(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    
p0 = [1930., 0.,4000.]

coeff, var_matrix = curve_fit(Gaussian,lengtelijst, hoogtelijst, p0 = p0) 
fit = Gaussian(np.arange(0,10000, 10), *coeff)
plt.plot(lengtelijst, hoogtelijst)
plt.plot(np.arange(0,10000, 10), fit, color = 'red')
plt.show()
print('Fitted mean = ', coeff[1])
print('Fitted standard deviation = ', coeff[2])

def namaakbed(lengte):
    bedhoogte = 0.
    array = np.array([])
    if(type(lengte) == type(5.) or type(lengte) == type(5)):
        if (lengte < 1700.):
            bedhoogte = 1600.
        elif (1700. <= lengte < 3500.):
            bedhoogte = 1600. + Gaussian(lengte, -600., 2600., 300.)
        elif (3500. <=lengte < 6000.):
            bedhoogte = 1475.
        elif (6000. <= lengte < 7200.):
            bedhoogte = -lengte*(1475-900.)/(7400-6000)+3939.
        elif (lengte >= 7400):
            bedhoogte = 900 + Gaussian(lengte, 600, 7400, 400.)
        return bedhoogte
    elif('numpy' in str(type(lengte))):
        for el in range(0,len(lengte)):
            if (0 <= lengte[el] < 1700.):
                bedhoogte = 1600.
            elif (1700. <= lengte[el] < 3500.):
                bedhoogte = 1600. + Gaussian(lengte[el], -600., 2500., 600.)
            elif (3500. <=lengte[el] < 6000.):
                bedhoogte = 1475.
            elif (6000. <= lengte[el] < 7400.):
                bedhoogte = -lengte[el]*(1475-900.)/(7400-6000)+3939.
            elif (lengte[el] >= 7400):
                bedhoogte = 300 + Gaussian(lengte[el], 600, 7400, 400.)
            array = np.append(array, bedhoogte)
            bedhoogte = 0.
        return array
    else:
        print("Geef een array of int/float op.")

print(np.polyfit([6000,7400],[1475,900],1))
lijstlengtes = np.arange(0, 10000,10)
plt.plot(lijstlengtes,namaakbed(lijstlengtes))
plt.ylim(0, 1700)
plt.show()

#plt.plot(ice)
#plt.plot(bed)
#plt.show()



#git pull origin master -m "hallo"
#git commit -m "hallo"
#git push origin master -m "hallo"

