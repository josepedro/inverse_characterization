# -*- coding: utf-8 -*-
from scipy.optimize import differential_evolution
import numpy
import scipy.io
import math

'''
#best1exp
#f = 0.935
#c_r = 1
#tamanho_populacao = 100
'''

def ackley(x):
	arg1 = -0.2 * numpy.sqrt(0.5 * (x[0] ** 2 + x[1] ** 2))
	arg2 = 0.5 * (numpy.cos(2. * numpy.pi * x[0]) + numpy.cos(2. * numpy.pi * x[1]))
	return -20. * numpy.exp(arg1) - numpy.exp(arg2) + 20. + numpy.e

def objetivo_allard_limp_funcao_dupla(vector):
    params_allard_rigido_simples = allard_rigido(0.025,2000,vector[0],vector[1],vector[2],vector[3],vector[4])
    Zs1 = params_allard_rigido_simples['Zs']
    params_allard_rigido_duplo = allard_rigido(2*0.025,2000,vector[0],vector[1],vector[2],vector[3], vector[4])
    Zs2 = params_allard_rigido_duplo['Zs']

    mat = scipy.io.loadmat('Z_exp.mat')
    
    Z_referencia_esp1 = mat['Z_3']
    F_obj_1 = (numpy.absolute(numpy.subtract((Z_referencia_esp1[500-1:1650-1]/413),(Zs1[500-1:1650-1]/413) ) )**2)
    F_obj_1 = F_obj_1.sum()

    Z_referencia_esp2 = mat['Z_6']
    F_obj_2 = (numpy.absolute(numpy.subtract((Z_referencia_esp2[500-1:1650-1]/413),(Zs2[500-1:1650-1]/413) ) )**2)
    F_obj_2 = F_obj_2.sum()

    print F_obj_1
    print F_obj_2

    F_obj = 3*F_obj_1 + F_obj_2

    #print 'F_obj'
    #print F_obj

    return F_obj

def allard_rigido(L, fmax, sigma, phi, alpha_inf, Lambda_material, fator):
    '''
    ## Modelo de material poroso de estrutura rigido
    # Calcula densidade efetiva dinamica e compressibilidade efetiva para
    # modelos de fluido equivalente de Johnson-Lafarge. 
    # [Zs alpha] = Lafarge(L,freq,sigma,phi,alpha_inf,Lambda_material,Lambda_l_material)
    # Paramentros de entrada:
    # L = comprimento da amostra;
    # freq = frequencia maxima;
    # sigma = resistividade ao fluxo [Ns/m^2];
    # phi = porosidade [%]
    # alpha_inf = tortuosidade (normalmente varia de 1 a 4)
    # Lambda_material = Comprimento caracteristico viscoso
    # Lambda_l_material = Comprimento caracteristico termico

    # Como saida tem-se:
    # Zs = impedancia de superficie
    # Alpha = Coeficiente de absorcao
    '''

    L_material = numpy.complex256(L); # espessura do material poroso
    rho_0 = numpy.complex256(1.204);      # Densidade [kg/m^3]
    T = numpy.complex256(20);               # Temperatura
    c_0 = numpy.complex256(331.2+0.6*T);# velocidade de propagacao no meio [m/s]
    eta = numpy.complex256(1.84e-5);      # Coeficiente de viscosidade, viscosidade absoluta ou viscosidade dinamica
    gamma = numpy.complex256(1.4);        # Razão de calores específicos (gamma=c_p/c_v)
    k_f = numpy.complex256(0.026);        # Condutibilidade térmica do fluido [W/mK]
    P_0 = numpy.complex256(101320);       # Pressão estática do meio [Pa]
    c_p = numpy.complex256(1.0035e3);     # Coeficiente de Calor específico a pressão constante (c_p), (c_v é coef. calor esp. a VOLUME cte)
    Pr = (eta*c_p)/k_f;   # Número de Prandtl (???)--- c_p calor específico, k_f é

    # Parâmetros Macro
    sigma_material = sigma;      # Resistividade ao fluxo [kN/s]
    phi_material =  phi;         # Porosidade 60#
    alpha_inf_material = alpha_inf;    # Tortuosidade
    Lambda = Lambda_material;
    #fator = 1
    Lambda_l = Lambda_material*fator 
    #Lambda_l = Lambda_l_material;
    q_0 = eta/sigma;
    # q_l_0 = q_0*Lambda_l/Lambda;
    q_l_0 = (phi_material*Lambda_l*Lambda_l)/8.;

    f = numpy.arange(0, fmax, 1)

    w = 2*math.pi*f ## MODELO DE Johnson - Lafarge

    D = (eta*(phi_material**2)*(Lambda**2));
    C = (4*w*rho_0*(q_0**2)*(alpha_inf_material**2));
    #numpy.divide(a,b,dtype=float)
    G = numpy.nan_to_num((1j*w*rho_0*alpha_inf_material*q_0))
    B = numpy.nan_to_num(numpy.divide((phi_material*eta),G))
    #B = ((phi_material*eta)/(1j*w*rho_0*alpha_inf_material*q_0));
    E = ( 1 + 1j*numpy.divide(C, D));
    E = numpy.power(E,1/2)
    A = 1 + numpy.multiply(B, E)
    rho_rigido = rho_0*alpha_inf_material * A; # densidade dinâmica efetiva

    print 'rho_rigido'
    print numpy.mean(rho_rigido)

    F = numpy.complex256(eta*(phi_material**2)*(Lambda_l**2))
    E = numpy.complex256(numpy.multiply(4*rho_0*Pr*(q_l_0**2), w))
    D =  1 + (1j* numpy.divide(E, F))
    G = (1j*w*rho_0*Pr*q_l_0)
    C = numpy.divide(numpy.complex256(phi_material*eta), numpy.complex256(G))
    C = numpy.nan_to_num(C)
    # C e D tao errados
    B = 1 + numpy.multiply(C, numpy.power(D,(1/2)))
    A = numpy.divide(gamma-(gamma-1), B)
    K_ef = numpy.nan_to_num(numpy.divide(gamma*P_0, A)) # módulo de compressibilidade 
    print 'A'
    print numpy.mean(A)

    #print 'compressibilidade'
    #print numpy.mean(K_ef)

    Zc = numpy.sqrt(numpy.multiply(rho_rigido, K_ef))
    A = numpy.divide(rho_rigido, K_ef)
    B = numpy.power(A,(1/2))
    kc_L = numpy.multiply(w, B)
    A = numpy.multiply(kc_L, L_material)
    H = numpy.nan_to_num(numpy.multiply(phi_material,numpy.arctan(A)))
    Zs = numpy.nan_to_num(-1j*numpy.divide(Zc, H))
    
    alpha = 1 - numpy.power((numpy.absolute(numpy.divide((Zs - rho_0*c_0),(Zs + rho_0*c_0)))), 2);

    params_allard_rigido = {'Zs': Zs, 'alpha': alpha}

    return params_allard_rigido

#bounds = [(-5, 5), (-5, 5)]
# Lambda_material = Comprimento característico viscoso
# Lambda_l_material = Comprimento característico térmico
# Lambda_l = Lambda_l_material;
# Lambda_l = Lambda_material*fator
# Comprimento característico térmico = (Comprimento característico viscoso)*fator
'''
paramDefCell = {'parameter1', [5000 300000], 10     # RESISITIVDADE AO FLUXO
            'parameter2', [0.90   0.99], 0.01   # POROSIDADE
            'parameter3', [1 4], 0.01           # TORTUOSIDADE
            'parameter4', [10e-6 500e-6], 1e-6  # COMRPIMENTO CARACTERÍSTICO VISCOSO
            'parameter5', [10e-6 500e-6], 1e-6  # COMPRIMENTO CARACTERÍSTICO TÉRMICO ou FATOR (1 3)
            'parameter6', (64, 64), 0.01 };      # DENSIDADE
'''

#bounds = [ (5000, 300000), (0.90, 0.99), (1, 4), (10e-6, 500e-6), (1,3), (64, 64)]
#result = differential_evolution(ackley, bounds, strategy='best1exp', disp=True, popsize=100, mutation=0.935, recombination=1)

#print result.x 
#print result.fun

vector = numpy.array([24930, 0.91, 1, 4.1000e-05, 2.847222222222222, 64])
f_obj = objetivo_allard_limp_funcao_dupla(vector)

print f_obj

'''
compressibilidade da errado

compressibilidade
(5.34825643657e+12-116.617165994e+04j)

1.1110e+05 + 1.1826e+04i

'''