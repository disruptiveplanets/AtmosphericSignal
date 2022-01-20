#Written by Prajwal Niraula
#Program to calculate the atmospheric signal in TSM metric specified by Kempton 2018
#Python version 3.7
#This makes use of forecaster tool from Chen and Kipping (2017)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mr_forecast as mr


import matplotlib as mpl
mpl.rc('font',**{'sans-serif':['Helvetica'], 'size':15,'weight':'bold'})
mpl.rc('axes',**{'labelweight':'bold', 'linewidth':1.5})
mpl.rc('ytick',**{'major.pad':22, 'color':'k'})
mpl.rc('xtick',**{'major.pad':10,})
mpl.rc('mathtext',**{'default':'regular','fontset':'cm','bf':'monospace:bold'})
#mpl.rc('text', **{'usetex':True})
#mpl.rc('text.latex',preamble=r'\usepackage{cmbright},\usepackage{relsize},'+r'\usepackage{upgreek}, \usepackage{amsmath}')
mpl.rc('contour', **{'negative_linestyle':'solid'})


MJup = 1.898e30 #mass of jupiters in grams
RJup = 7.1492e9 #
RadRatio = 10.973 #Radius of Jupiter to that of earth
MassRatio = 317.83 #Mass ratio of Jupiter to that of the earth
GravConst=6.6260755e-8
au=1.4959789e13 #cm
RSun=6.9599e10 #cm
k=1.380658e-16 #CGS Unit
Albedo = 0.0 #albedo
mu = 20.0
EarthMass = 5.97219e27  # grams
EarthRadius = 6.371e8   # radius of earth in cm


#Read the data using pandas
Data = pd.read_csv("PS_2022.01.20_06.34.56.csv", skiprows=22)

#Assign the right variables
Name = Data["pl_name"]
Period = Data["pl_orbper"] #period day
a_Rs = Data["pl_ratdor"]
PlanetaryMass = Data["pl_bmasse"]       #Planetary mass in Earth mass is here   
PlanetaryRad = Data["pl_rade"]          #planet radius in terms of r_earth
TStar = Data["st_teff"]                 #Stellar temperature
RStar = Data["st_rad"]                     #stellar radius in terms of r_sun
JMag = Data["sy_jmag"]                  #J mag star


#Remove the first two 11 Com b & 11 UMi b
Name[0] = "TOI-4306 b"
Period[0] =  2.722                       #period day
a_Rs[0] = 26.21
PlanetaryMass[0] = 2.38                  #Planetary mass in Earth mass is here   
PlanetaryRad[0] = 1.339                  #planet radius in terms of r_earth
TStar[0] = 2850                          #Stellar temperature
RStar[0] = 0.1539                        #stellar radius in terms of r_sun
JMag[0] =  12.258                        #J mag star

Name[1] = "SPECULOOS-2 b"
Period[1] =  2.722                       #period day
a_Rs[1] = 55.69
PlanetaryMass[1] = 2.31                  #Planetary mass in Earth mass is here   
PlanetaryRad[1] = 1.353                  #planet radius in terms of r_earth
TStar[1] = 2850                          #Stellar temperature
RStar[1] = 0.1539                        #stellar radius in terms of r_sun
JMag[1] =  12.258                        #J mag star



RemoveIndex = np.logical_or(PlanetaryRad>1.8, ~np.isfinite(JMag))


#Remove the right index
Name = np.ascontiguousarray(Name[~RemoveIndex])
Period = np.ascontiguousarray(Period[~RemoveIndex]) #period day
a_Rs = np.ascontiguousarray(a_Rs[~RemoveIndex])
PlanetaryMass = np.ascontiguousarray(PlanetaryMass[~RemoveIndex])       #Planetary mass in Earth mass is here   
PlanetaryRad = np.ascontiguousarray(PlanetaryRad[~RemoveIndex])         #planet radius in terms of r_earth
TStar = np.ascontiguousarray(TStar[~RemoveIndex])                                             #Stellar temperature
RStar = np.ascontiguousarray(RStar[~RemoveIndex])                                             #stellar radius in terms of r_sun
JMag = np.ascontiguousarray(JMag[~RemoveIndex])                                               #J mag star



#Use Chen and kipping to estimate the mass
for i in range(len(Name)):
    print("The current target is::", i)
    #Generate mass from the radius    
    if np.isfinite(PlanetaryMass[i]) and not(np.isfinite(PlanetaryRad[i])):
        pass
        #PlanetaryRad[i], _, _ = mr.Mstat2R(PlanetaryMass[i], 1e-9, sample_size=100, unit='Earth')
        #print("The planetary radius is :", PlanetaryRad[i], " and the mass is ", PlanetaryMass[i])  
    elif np.isfinite(PlanetaryRad[i]) and not(np.isfinite(PlanetaryMass[i])):
        #print("The value of given by ::", Name[i])
        PlanetaryMass[i], _, _ = mr.Rstat2M(PlanetaryRad[i], 1e-9, sample_size=1000, grid_size=1000,unit='Earth')  
        #print("The planetary mass is :", PlanetaryMass[i], " for the radius of the mass ", PlanetaryRad[i])
    else:
        print("Skipping this for now")




#Calculate the effective temperature in the planet.
T_Eff_Planet = (1-Albedo)**0.25*TStar*(1./(2.0*a_Rs))**0.5

#Calculate the gravity of the planet
Grav_pl = GravConst*PlanetaryMass*EarthMass/(PlanetaryRad*EarthRadius)**2




#Calculate the transmission spectroscopy signal
h_eff = 7*k*T_Eff_Planet/(mu*1.67e-24*Grav_pl)
S = 2*PlanetaryRad*EarthRadius*h_eff/(RStar*RSun)**2*1e6 #In terms of ppm


print("The maximum value of S is given by ")


SNR = S*10**(-JMag)

#Remove the 
RemoveIndex = ~np.isfinite(SNR)


Name = np.ascontiguousarray(Name[~RemoveIndex])
T_Eff_Planet = np.ascontiguousarray(T_Eff_Planet[~RemoveIndex])
S = np.ascontiguousarray(S[~RemoveIndex])
SNR = np.ascontiguousarray(SNR[~RemoveIndex])
PlanetaryRad = np.ascontiguousarray(PlanetaryRad[~RemoveIndex])
Grav_pl = Grav_pl[~RemoveIndex]

SelectIndex = np.where(Name=="TRAPPIST-1 b")
SNR = SNR/SNR[SelectIndex]

#Arrange by the SNR
#ArrangeIndex = np.argsort(SNR)[::-1]


#Name = np.ascontiguousarray(Name[ArrangeIndex])
#T_Eff_Planet = np.ascontiguousarray(T_Eff_Planet[ArrangeIndex])
#S = np.ascontiguousarray(S[ArrangeIndex])
#SNR = np.ascontiguousarray(SNR[ArrangeIndex])
#PlanetaryRad = np.ascontiguousarray(PlanetaryRad[ArrangeIndex])
#Grav_pl = Grav_pl[ArrangeIndex]


#Now the get the size of the signal for JWST
print("The gravity of the planet is given by:: ", Grav_pl[SelectIndex])


print("The first case:", Name[0])
print("The second case:", Name[1])









#The color bar is given by SNR
plt.figure(figsize=(12,8))

plt.scatter(T_Eff_Planet[0], S[0], s=(PlanetaryRad[0]*100), c="None", edgecolor='navy', lw=2, zorder=20)
plt.scatter(T_Eff_Planet[1], S[1], s=(PlanetaryRad[1]*100), c="None", edgecolor='navy', lw=2, zorder=20)


plt.scatter(T_Eff_Planet[SelectIndex], S[SelectIndex], s=(PlanetaryRad[SelectIndex]*10000)**0.5, c="None", edgecolor='black', lw=2, zorder=20)
#plt.text(T_Eff_Planet[SelectIndex]+30, S[SelectIndex], "TRAPPIST-1 b", va="center", color="black")

#for 1.5 vs 1.0 Earth radii
plt.scatter(975, 230, s=(1.5*100), c="k", edgecolor='None', lw=2, zorder=20)
plt.text(1000,  230, "1.5 $R_{\oplus}$", va="center")

plt.scatter(975, 220, s=(100), c="k",  edgecolor='None', lw=2, zorder=20)
plt.text(1000,  220, "1.0 $R_{\oplus}$", va="center")

RightHandSide = np.array(["SPECULOOS-2 b", "TRAPPIST-1 f", "TRAPPIST-1 e", "TRAPPIST-1 h", "K2-315 b",  "K2-72 b", "LHS 1140 c"])

for i in range(len(Name)):
    if S[i]>75.0 and T_Eff_Planet[i]<1000:
        if not(np.sum(Name[i]==RightHandSide)==1):
            plt.text(T_Eff_Planet[i]+20, S[i], Name[i], va="center")
        else:
            plt.text(T_Eff_Planet[i]-20, S[i], Name[i], ha="right", va="center")

    print(i, ":", Name[i])

Scatter = plt.scatter(T_Eff_Planet, S, c=SNR, s=(PlanetaryRad*100), edgecolor='None', cmap='rainbow', vmax=1.05, vmin=0.001)
cbar = plt.colorbar()
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.set_ylabel("Relative SNR", fontsize=20, rotation=270)
plt.ylabel("Expected Transmission Signal (ppm)", fontsize=20)
plt.xlabel("Equilibrium Temperature (K)", fontsize=20)
plt.ylim(5, 250)
plt.xlim(10, 1100)
#plt.yscale('log')
plt.yticks([12.5, 25., 50., 100., 200.,])
plt.tight_layout()
plt.savefig("SPECULOOS_2b.png")
plt.savefig("SPECULOOS_2b.pdf")
plt.show()