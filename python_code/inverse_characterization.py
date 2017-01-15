# inversse characterization
# -*- coding: utf-8 -*-
from scipy.optimize import differential_evolution
import numpy as np
import scipy.io
#import math
import matplotlib.pyplot as plt
'''
#best1exp
#f = 0.935
#c_r = 1
#tamanho_populacao = 40
'''
#%%
def vector_frequencies(fmax,fmin,deltaf):

    f = np.arange(fmin,fmax+deltaf,deltaf)
    
    omega=2*np.pi*f
    return [f,omega]
        
#%%
def objetivo_allard_limp_funcao_dupla(x):
    params_JC_simples = m_Johnson_Champoux3_deltaf(fmax,fmin,deltaf,espessura1,x[0],x[1],x[2],x[3],x[4],x[5],structure='rigid')
    A1 = params_JC_simples['alpha']
    params_JC_duplo = m_Johnson_Champoux3_deltaf(fmax,fmin,deltaf,espessura2,x[0],x[1],x[2],x[3],x[4],x[5],structure='rigid')
    A2 = params_JC_duplo['alpha']


    F_obj_1 = np.absolute(A_referencia_esp1[600:5000+1]-A1[600:5000+1])**2
#    F_obj_1 = np.absolute(Z_referencia_esp1[500:1650+1]/413-Zs1[500:1650+1]/413)**2
    F_obj_1 = F_obj_1.sum()


    F_obj_2 = np.absolute(A_referencia_esp2[600:5000+1]-A2[600:5000+1])**2
#    F_obj_2 = np.absolute(Z_referencia_esp2[500:1650+1]/413-Zs2[500:1650+1]/413)**2
    F_obj_2 = F_obj_2.sum()

    F_obj = F_obj_1 + F_obj_2
 
    return F_obj

#%%
def m_Johnson_Champoux3_deltaf(fmax,fmin,deltaf,espessura,sigma,phi,alpha_inf,lambda_v,fator,rho,structure='rigid'):
           
    rho_ef_jc = rho0*alpha_inf*(1+((phi*sigma)/(1j*rho0*alpha_inf*omega))*((1+1j*((4*eta0*rho0*alpha_inf**2*omega)/(sigma**2*phi**2*lambda_v**2)))**(1./2.)))
    rho_eq_jc=rho_ef_jc

    if structure=='limp':
        rho_eq_jc = (rho_eq_jc*rho-rho0**2)/(rho+rho_eq_jc-2*rho0)

    K_ef_jc = gamma*P0/(gamma-((gamma-1)*(1+((8*eta0)/(1j*omega*Pr*(fator*lambda_v)**2*rho0))*(1+1j*((omega*Pr*rho0*(fator*lambda_v)**2)/(16*eta0)))**(1./2.))**(-1)))
    
    
    k=omega*(rho_eq_jc/K_ef_jc)**(1./2.)
    Zc=(rho_eq_jc*K_ef_jc)**(1./2.)
    
    Zs=-(1j/phi)*Zc*(1/np.tan(k*espessura))
    R=(Zs-rho0*c0)/(Zs+rho0*c0)
    Abs=1-np.absolute(R)**2

    params_JC = {'Zs': Zs, 'alpha': Abs}
    
    print fator

    return params_JC
    

if __name__ == "__main__":

	#Air properties
	P0 = 101325.0                   #static pressure [Pa]
	R_spec = 287.058                #specific gas constant [J/(kg·K)]
	T = 20.0                        #temperature [°C]
	rho0 = P0/(R_spec*(T+273.15))      #air density [kg/m^3] (ideal gas law)
	c0 = 331.3+0.606*T                  #speed of sound in dry air [m/s]
	eta0 = ((1.458e-6)*(T+273.15)**(3./2.))/((T+273.15)+110.4)    #dynamic viscosity [kg/(m*s)] (Sutherland's law)
	k_f = ((2.334e-3)*(T+273.15)**(3./2.))/((T+273.15)+164.54)      #thermal conductivity [W/(m*K)]
	gamma = 1.4                         #heat capacity ratio (ratio of specific heats)
	c_p = 1.0035e3                      #specific heat capacity at constant pressure [J/(kg*K)]
	Pr = eta0*c_p/k_f                   #Prandtl number

	# Setup parameters
	fmax=5500
	fmin=0
	deltaf=1
	[f,omega]=vector_frequencies(fmax,fmin,deltaf)
	espessura1=0.024
	espessura2=2*0.024

	simples = scipy.io.loadmat('../data/espuma_preta/espumapreta_24mm_simples.mat')
	duplo = scipy.io.loadmat('../data/espuma_preta/espumapreta_48mm_duplo.mat')
	A_referencia_esp1 = simples['A_A1']
	A_referencia_esp1=np.ravel(A_referencia_esp1)
	A_referencia_esp2 = duplo['A_A1']
	A_referencia_esp2=np.ravel(A_referencia_esp2)


	bounds = [ (5000, 300000), (0.80, 0.99), (1, 4), (10e-6, 500e-6), (1,6), (64, 64)]
	result = differential_evolution(objetivo_allard_limp_funcao_dupla, bounds, strategy='best1exp', disp=True, popsize=20, mutation=0.935, recombination=1)

	print result.x
	print result.fun

	params_JC_otim1=m_Johnson_Champoux3_deltaf(fmax,fmin,deltaf,espessura1,result.x[0],result.x[1],result.x[2],result.x[3],result.x[4],result.x[5],structure='rigid')
	params_JC_otim2=m_Johnson_Champoux3_deltaf(fmax,fmin,deltaf,espessura2,result.x[0],result.x[1],result.x[2],result.x[3],result.x[4],result.x[5],structure='rigid')
	alpha_otim1=params_JC_otim1['alpha']
	alpha_otim2=params_JC_otim2['alpha']


	plt.figure(1)
	plt.plot(f[200:5500+1],A_referencia_esp1[200:5500+1],f[200:5500+1],A_referencia_esp2[200:5500+1],f[200:5500+1],alpha_otim1[200:5500+1],f[200:5500+1],alpha_otim2[200:5500+1])
	plt.xlabel('frequency')
	plt.ylabel('alpha')
	plt.grid()
	plt.show()
