

import numpy as np
import pandas as pd
# import campbell
import sys
import filehandling
import time
import mms_funcs
import copy


#====================================================
# User input
#====================================================
calc_type = 'jit'
variable_dt = False

#====================================================
# Load and ready rainfall data
#====================================================
# Set minimum timestep
dt_min = 0.1 #min

# load rain
rain = filehandling.load_rain('Regn_viby_main.xlsx') #2022
# rain = filehandling.load_rain('Input_main2021.xlsx') #2021

# load Epot
E_pot = filehandling.load_rain2('Epotmain.xlsx') #2022

# load VWC målinger
VWC_10 = filehandling.load_rain('Jord_main_10.xlsx') #2022
VWC_10_value = VWC_10.VWC.values #2022
VWC_10_date = VWC_10.dato #2022
# VWC_10_value = rain.VWC.values #2021
# VWC_10_date = rain.VWC_date #2021


VWC_20 = filehandling.load_rain('Jord_main_20.xlsx')
VWC_20_value = VWC_20.VWC.values 
VWC_20_date = VWC_20.dato 


#Load dato
dato_ramme = filehandling.load_rain('Date_main.xlsx') #2022
dato_date = dato_ramme.main #2022
# dato_date = rain.Dato #2021


# number of simulations
n_sim = len(rain)
rainfall = np.zeros(n_sim)
rainfall = rain.Intensity.values #2022
# rainfall = rain.regn.values #2021
ET_p = E_pot.E_pot.values #2022
# ET_p = rain.epot.values #2021
ETa = np.zeros(n_sim)  

#====================================================
# Load and initialize model
#====================================================
# Load model from excel file
model = filehandling.load_model('box_builder_viby.xlsx')


# Find number of layers in loaded model
n_layers = np.int64(model.shape[0])

# Setup result arrays
theta = np.zeros((n_layers, n_sim)) # Volumetric water content cm3 / cm3
K     = np.zeros((n_layers, n_sim)) # Current hydrualic conductivity cm / min
psi   = np.zeros((n_layers, n_sim)) # cmH2O
alpha_L = np.zeros((n_layers, n_sim)) 
mmsK_L = np.zeros((n_layers, n_sim))
v     = np.zeros((n_layers-1, n_sim)) # Velocity between each box


# Set initial values
theta[0,0] = 0.205
theta[1,0] = 0.206
theta[2,0] = 0.308
theta[3,0] = 0.308
theta[4,0] = 0.309
theta[5,0] = 0.309
theta[6,0] = 0.309
theta[7,0] = 0.309
theta[8,0] = 0.310
theta[9,0] = 0.310
theta[10,0] = 0.310
theta[11,0] = 0.310
theta[12,0] = 0.310
theta[13,0] = 0.310
theta[14,0] = 0.310
theta[15,0] = 0.310
theta[16,0] = 0.310
theta[17,0] = 0.310
theta[18,0] = 0.310
theta[19,0] = 0.310

# Create values for model parameters

if calc_type=='vec':
    psi[:,0]   = model.psi_e * (theta[:, 0] / model.theta_s) ** (-model.b)
    K[:,0] = model.Ks * (theta[:,0] / model.theta_s) ** (2 * model.b + 3)
else:
    for i in range(n_layers):
        psi[i,0]   = model.psi_e[i] * (theta[i, 0] / model.theta_s[i]) ** (-model.b[i])
        K[i,0] = model.Ks[i] * (theta[i,0] / model.theta_s[i]) ** (2 * model.b[i] + 3)




if calc_type=='vec': 
    alpha_L[:,0] = (2+3/model.b.values) / -psi[:,0]
    mmsK_L[:,0] = K[:,0] * np.exp(-alpha_L[:,0] * psi[:,0])
else:
    for i in range(n_layers):
        alpha_L[i,0] = (2+3/model.b.values[i]) / -psi[i,0]
        mmsK_L[i,0] = K[i,0] * np.exp(-alpha_L[i,0] * psi[i,0])
        
    

alpha_n = np.zeros((n_layers-1, n_sim))
K_n = np.zeros((n_layers-1, n_sim))
v_model = np.zeros((n_layers+1, n_sim))
dt_mult = np.zeros((n_layers-1, n_sim))
dt_masse = np.zeros((n_layers, n_sim)) 
theta_dt = np.zeros(n_layers)
theta_dt[0] = 0.205
theta_dt[1] = 0.206
theta_dt[2] = 0.308
theta_dt[3] = 0.308
theta_dt[4] = 0.309
theta_dt[5] = 0.309
theta_dt[6] = 0.309
theta_dt[7] = 0.309
theta_dt[8] = 0.310
theta_dt[9] = 0.310
theta_dt[10] = 0.310
theta_dt[11] = 0.310
theta_dt[12] = 0.310
theta_dt[13] = 0.310
theta_dt[14] = 0.310
theta_dt[15] = 0.310
theta_dt[16] = 0.310
theta_dt[17] = 0.310
theta_dt[18] = 0.310
theta_dt[19] = 0.310


# Evapotranspiration
theta_mktop = 0.21 #field capacity cm3 / cm3
theta_vgtop = 0.04 #Wilting point cm3 / cm3
WDm = theta_mktop-theta_vgtop
WDa = np.zeros(n_sim) 


OverfladeA_i = np.zeros(n_sim)
OverfladeA_ii = np.zeros(n_sim)# Surface runoff 


# Subsurface runoff 
theta_psiCW_top = 0.43*(-2.4/-7.5)**(1/5.2) # Vandindhold ved pF -7.5, kriterie for runoff
theta_psiCW_kom = 0.36*(-15/-50)**(1/12.69) # Vandindhold ved pF -50, kriterie for runoff
UnderjordiskA_i = np.zeros(n_sim)
UnderjordiskA_ii = np.zeros(n_sim)
UnderjordiskA_iii = np.zeros(n_sim)


#====================================================
# Run model!
#====================================================
t_start = time.time()

# counter variables - needed because of non-consistant timestep
ctr = 1
n = 1
print('Running Model!')

while n<n_sim:
    # Update box values
    if calc_type=='jit':
        psi[:, ctr], K[:, ctr], alpha_L[:, ctr], mmsK_L[:, ctr] = mms_funcs.box_calc_vec(psi[:, ctr-1], model.psi_e.values, theta[:, ctr-1], model.theta_s.values, model.b.values, K[:, ctr-1], model.Ks.values, alpha_L[:, ctr-1], mmsK_L[:,ctr-1])
    else:
        psi[:, ctr], K[:, ctr], alpha_L[:, ctr], mmsK_L[:, ctr] = mms_funcs.box_calc(psi[:, ctr-1], model.psi_e.values, theta[:, ctr-1], model.theta_s.values, model.b.values, K[:, ctr-1], model.Ks.values, alpha_L[:, ctr-1], mmsK_L[:,ctr-1])
    
    # Update values between
    if calc_type=='jit':
        v[:,ctr] = mms_funcs.darcy_v_jit(n_layers, alpha_n[:, ctr-1], alpha_L[:, ctr], K_n[:, ctr-1], mmsK_L[:, ctr], psi[:, ctr], model.dz.values, v[:,ctr])    
    else:
        v[:,ctr] = mms_funcs.darcy_v(n_layers, alpha_n[:, ctr-1], alpha_L[:, ctr], K_n[:, ctr-1], mmsK_L[:, ctr], psi[:, ctr], model.dz.values, v[:,ctr])    
        
    # Set added water to each box
    v_model[0, ctr] = rainfall[n];
    v_model[1:-1, ctr] = v[:,ctr]
    v_model[-1, ctr]= v[-1, ctr]
      
    
    ## Surface Runoff ##
    # Surface runoff 2: theta > theta_s
    if theta_dt[0]>model.theta_s[0]:
        OverfladeA_ii[n-1] = v_model[0,n-1]*dt_min/model.dz-(model.theta_s[0]-theta_dt[0])
        v_model[0] = 0
    else:
        # Surface runoff 1: v>1.5*Ks
        if v_model[0,n-1]>K[0,n-1]:
              OverfladeA_i[n-1]=v_model[0,n-1]-K[0,n-1]
              v_model[0,n] = K[0,n-1]
        else:
              OverfladeA_i[n-1]=0
        OverfladeA_ii[n-1] = 0
    
    # Subsurface Runoff ##
    # Subsurface Runoff 3: theta > theta_s -> første kompakt (boks 9)
    if theta_dt[2]>model.theta_s[2]:
        UnderjordiskA_iii[n-1] = v_model[2,n-1]*dt_min/2.5-(model.theta_s[2]-theta_dt[2])
        v_model[2] = 0
    else:
         if theta_dt[1]>theta_psiCW_top and theta_dt[2]>theta_psiCW_kom:
             UnderjordiskA_i[n-1] = theta_dt[1]-theta_psiCW_top
         else:
          # Subsurface Runoff 2: v>K, boks 3 -> første kompakt 
              if v_model[2,n-1]>K[2,n-1]:
                  UnderjordiskA_ii[n-1] = v_model[2,n-1]-K[2,n-1]
                  v_model[2] = K[2,n-1]
              else:    
                  UnderjordiskA_ii[n-1] = 0
              UnderjordiskA_i[n-1] = 0
         UnderjordiskA_iii[n-1] = 0


 
    # Evapotranspiration (boks 1-8) ##  
    if sum(theta[0:2,n-1])/2>theta_mktop:
        WDa[n-1] = 0
    else:
        WDa[n-1] = theta_mktop-sum(theta[0:2,n-1])/2

    if WDa[n-1] <=((theta_mktop-theta_vgtop)*0.5):
        ETa[n-1]=(ET_p[n-1]/2)/2.5
    else:
        ETa[n-1]=(2*(1-(WDa[n-1]/WDm)))*ET_p[n-1]/2/2.5 

   

    # Update water content in each box
    if calc_type=='jit':
        theta_dt = mms_funcs.update_jit(theta_dt, theta[:, ctr-1], v_model[:,ctr-1], dt_min, model.dz.values, n_layers)
    else:
        theta_dt = mms_funcs.update(theta_dt, theta[:, ctr-1], v_model[:,ctr-1], dt_min, model.dz.values, n_layers)
    
    
   
    if any(theta_dt<=0):
        print(v_model[:, -1])
        print(theta_dt)
        sys.exit('Found values below zero')
    
    #  Update water content with runoff and evapotranspiration 
    theta[0:2, ctr] = theta_dt[0:2]-ETa[n-1]
    theta[2:20, ctr] = theta_dt[2:20]
    
    ctr+=1
    
    n += 1
    #break
        


#print(result_c)
t_total = time.time() - t_start
print('Sim. tid:',t_total)


print('Fordampning:', sum(ETa)*2*25/10) #cm fordampning på simulationstiden
print('v>k overflade:',sum(OverfladeA_i*0.1)) #cm overfladeafstrømning v>K
print('Theta>theta_s overflade:',sum(OverfladeA_ii)) #cm overfladeafstrømning theta>theta_s 
print('underjordisk kritisk VWc:',sum(UnderjordiskA_i*2.5)) #cm Underjordiskafstrømning psi = -7,5
print('underjordisk v>k:',sum(UnderjordiskA_ii*0.1)) #cm Underjordiskafstrømning v>K
print('underjordisk theta>theta_s:',sum(UnderjordiskA_iii)) #cm Underjordiskafstrømning theta>theta_s 

afstromning = (OverfladeA_i*0.1)+(OverfladeA_ii)+(UnderjordiskA_i*2.5)+(UnderjordiskA_ii*0.1)+(UnderjordiskA_iii) 
print('Samlet afstrømning:',sum(afstromning)) #samlet afstrømning cm
print('Samlet nedbør:',sum(rainfall)) #



import matplotlib.pyplot as plt
# VWC vs simuleringer
plt.plot(dato_date,(theta[3]+theta[4])/2)
plt.plot(VWC_10_date,VWC_10_value)
# plt.plot(VWC_20_date,VWC_20_value)
plt.ylabel('Vandindhold [cm^3 /cm^3]')
plt.xlabel('Dato')
plt.xticks(rotation=45)
plt.title('11/4 2021')
plt.legend(['Modelleret','10cm','20cm'],loc = "center left")
plt.show()

plt.plot(dato_date,(psi[3]+psi[4])/2)
plt.plot(dato_date, psi[1])
# plt.plot(dato_date, psi[2])
plt.ylabel('Psi [cmH2O]')
plt.xlabel('Dato')
plt.xticks(rotation=45)
plt.title('11/4 2021')
plt.legend(['Modelleret 10cm','3,75cm','6,25'],loc = "lower right")
plt.show()

plt.plot(dato_date,np.cumsum(afstromning*10))
plt.ylabel('Akk. afs')
plt.xlabel('Dato')
plt.xticks(rotation=45)
plt.title('Afstrømning  akkumuleret')
plt.show()






