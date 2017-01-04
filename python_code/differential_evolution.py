# -*- coding: utf-8 -*-
"""
Spyder Editor

Este é um arquivo de script temporário.
"""
import numpy
import random
import math
import scipy.io

class differential_evolution_optimizer(object):
  """
This is a python implementation of differential evolution
It assumes an evaluator class is passed in that has the following
functionality
data members:
 n              :: The number of parameters
 domain         :: a  list [(low,high)]*n
                   with approximate upper and lower limits for each parameter
 x              :: a place holder for a final solution

 also a function called 'target' is needed.
 This function should take a parameter vector as input and return a the function to be minimized.

 The code below was implemented on the basis of the following sources of information:
 1. http://www.icsi.berkeley.edu/~storn/code.html
 2. http://www.daimi.au.dk/~krink/fec05/articles/JV_ComparativeStudy_CEC04.pdf
 3. http://ocw.mit.edu/NR/rdonlyres/Sloan-School-of-Management/15-099Fall2003/A40397B9-E8FB-4B45-A41B-D1F69218901F/0/ses2_storn_price.pdf


 The developers of the differential evolution method have this advice:
 (taken from ref. 1)

If you are going to optimize your own objective function with DE, you may try the
following classical settings for the input file first: Choose method e.g. DE/rand/1/bin,
set the number of parents NP to 10 times the number of parameters, select weighting
factor F=0.8, and crossover constant CR=0.9. It has been found recently that selecting
F from the interval [0.5, 1.0] randomly for each generation or for each difference
vector, a technique called dither, improves convergence behaviour significantly,
especially for noisy objective functions. It has also been found that setting CR to a
low value, e.g. CR=0.2 helps optimizing separable functions since it fosters the search
along the coordinate axes. On the contrary this choice is not effective if parameter
dependence is encountered, something which is frequently occuring in real-world optimization
problems rather than artificial test functions. So for parameter dependence the choice of
CR=0.9 is more appropriate. Another interesting empirical finding is that rasing NP above,
say, 40 does not substantially improve the convergence, independent of the number of
parameters. It is worthwhile to experiment with these suggestions. Make sure that you
initialize your parameter vectors by exploiting their full numerical range, i.e. if a
parameter is allowed to exhibit values in the range [-100, 100] it's a good idea to pick
the initial values from this range instead of unnecessarily restricting diversity.

Keep in mind that different problems often require different settings for NP, F and CR
(have a look into the different papers to get a feeling for the settings). If you still
get misconvergence you might want to try a different method. We mostly use DE/rand/1/... or DE/best/1/... .
The crossover method is not so important although Ken Price claims that binomial is never
worse than exponential. In case of misconvergence also check your choice of objective
function. There might be a better one to describe your problem. Any knowledge that you
have about the problem should be worked into the objective function. A good objective
function can make all the difference.

Note: NP is called population size in the routine below.)
Note: [0.5,1.0] dither is the default behavior unless f is set to a value other then None.

  """

  def __init__(self,
               evaluator,
               population_size=100,
               f=0.935,
               cr=1,
               eps=1e-2,
               n_cross=1,
               max_iter=150,
               monitor_cycle=200,
               out=None,
               show_progress=False,
               show_progress_nth_cycle=1,
               insert_solution_vector=None,
               dither_constant=0.4):
    print 'quatro'
    self.dither=dither_constant
    self.show_progress=show_progress
    self.show_progress_nth_cycle=show_progress_nth_cycle
    self.evaluator = evaluator
    self.population_size = population_size
    self.f = f
    self.cr = cr
    self.n_cross = n_cross
    self.max_iter = max_iter
    self.monitor_cycle = monitor_cycle
    self.vector_length = evaluator.n
    self.eps = eps
    self.population = []
    self.seeded = False
    if insert_solution_vector is not None:
      assert len( insert_solution_vector )==self.vector_length
      self.seeded = insert_solution_vector
    for ii in xrange(self.population_size):
      self.population.append( numpy.zeros(self.vector_length) )


    self.scores = 1000*numpy.ones(self.population_size)
    self.optimize()
    self.best_score = numpy.min( self.scores )
    self.best_vector = self.population[ numpy.nanargmin( self.scores ) ]
    self.evaluator.x = self.best_vector

    if self.show_progress:
      self.evaluator.print_status(
            numpy.min(self.scores),
            numpy.mean(self.scores),
            self.population[ numpy.nanargmin( self.scores ) ],
            'Final')

  # ESSA FUNCAO EH O CERNE DE TODO CODIGO
  def optimize(self):
    print 'cinco'
    # initialise the population please
    self.make_random_population()
    # score the population please
    self.score_population()
    converged = False
    monitor_score = numpy.min( self.scores )
    self.count = 0
    # ESSE WHILE EH IMPORTANTE PRA CARALHO POIS ELE QUE VAI EVOLUIR AS POPULACOES ATE CONVERGIR
    while not converged:
      self.evolve()
      location = numpy.nanargmin( self.scores )
      if self.show_progress:
        if self.count%self.show_progress_nth_cycle==0:
          # make here a call to a custom print_status function in the evaluator function
          # the function signature should be (min_target, mean_target, best vector)
          self.evaluator.print_status(
            numpy.min(self.scores),
            numpy.mean(self.scores),
            self.population[ numpy.nanargmin( self.scores ) ],
            self.count)

      self.count += 1
      if self.count%self.monitor_cycle==0:
        if (monitor_score - numpy.min(self.scores) ) < self.eps:
          converged = True
        else:
         monitor_score = numpy.min(self.scores)
      rd = (numpy.mean(self.scores) - numpy.min(self.scores) )
      rd = rd*rd/(numpy.min(self.scores)*numpy.min(self.scores) + self.eps )
      if ( rd < self.eps ):
        converged = True


      if self.count>=self.max_iter:
        converged =True

  def make_random_population(self):
    print 'seis'
    for ii in xrange(self.vector_length):
      delta  = self.evaluator.domain[ii][1]-self.evaluator.domain[ii][0]
      offset = self.evaluator.domain[ii][0]
      random_values = numpy.random.uniform(size=self.population_size)
      random_values = random_values*delta+offset
      # now please place these values ni the proper places in the
      # vectors of the population we generated
      for vector, item in zip(self.population,random_values):
        vector[ii] = item
    if self.seeded is not False:
      self.population[0] = self.seeded

  def score_population(self):
    print 'sete'
    for vector,ii in zip(self.population,xrange(self.population_size)):
      tmp_score = self.evaluator.target(vector)
      print '------------------'
      print 'tamamho desse array de scores'
      if type(tmp_score) is numpy.ndarray:
        print len(tmp_score)
      print '------------------'
      self.scores[ii]=tmp_score

  # FUNCAO EVOLUIR MUITO IMPORTANTE
  def evolve(self):
    print 'evoluindo a galera'
    for ii in xrange(self.population_size):
      rnd = numpy.random.uniform(size=self.population_size-1)
      permut = numpy.argsort(rnd)
      # make parent indices
      i1=permut[0]
      if (i1>=ii):
        i1+=1
      i2=permut[1]
      if (i2>=ii):
        i2+=1
      i3=permut[2]
      if (i3>=ii):
        i3+=1
      #
      x1 = self.population[ i1 ]
      x2 = self.population[ i2 ]
      x3 = self.population[ i3 ]

      if self.f is None:
        use_f = random.random()/2.0 + 0.5
      else:
        use_f = self.f

      vi = x1 + use_f*(x2-x3)
      # prepare the offspring vector please
      rnd = numpy.random.uniform(size=self.vector_length)
      permut = numpy.argsort(rnd)
      test_vector = self.population[ii].copy()
      # first the parameters that sure cross over
      for jj in xrange( self.vector_length  ):
        if (jj<self.n_cross):
          test_vector[ permut[jj] ] = vi[ permut[jj] ]
        else:
          if (rnd[jj]>self.cr):
            test_vector[ permut[jj] ] = vi[ permut[jj] ]
      # get the score please
      test_score = self.evaluator.target( test_vector )
      # check if the score if lower
      if test_score < self.scores[ii] :
        self.scores[ii] = test_score
        self.population[ii] = test_vector


  def show_population(self):
    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    for vec in self.population:
      print list(vec)
    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"


class test_function(object):
  def __init__(self):
    self.x = None
    self.n = 9
    self.domain = [ (-100,100) ]*self.n
    self.optimizer =  differential_evolution_optimizer(self,population_size=100,n_cross=5)
    assert numpy.sum(self.x*self.x)<1e-5

  def target(self, vector):
    tmp = vector.copy()
    result = (numpy.sum(numpy.cos(tmp*10))+self.n+1)*numpy.sum( (tmp)*(tmp) )
    return result


class test_rosenbrock_function(object):
  def __init__(self, dim=5):
    print "tres"
    self.x = None
    self.n = 2*dim
    self.dim = dim
    # ESSE CARA AQUI DOMINIO
    self.domain = [ (1,3), (2,3) ]
    print '----------------------------------'
    print self.domain
    print '----------------------------------'
    # ESSE CARA AQUI PEGA O MINIMO PARA SETAR O TAMANHO DA POPULACAO population_size=min(self.n*10,40)
    self.optimizer =  differential_evolution_optimizer(self,population_size=min(self.n*10,40),n_cross=self.n,cr=0.9, eps=1e-8, show_progress=True)
    print list(self.x)
    for x in self.x:
      assert abs(x-1.0)<1e-2


  def target(self, vector):
    print 'passei na funcao principal'
    tmp = vector.copy()
    # AQUI NESSA BAGAÇA EH SOH UM VALOR DE X
    x_vec = vector[0:self.dim]
    # AQUI NESSA BAGAÇA EH SOH UM VALOR DE Y
    y_vec = vector[self.dim:]
    print 'printando os valores de x'
    print x_vec
    print 'printando os valores de y'
    print y_vec
    result=0
    for x,y in zip(x_vec,y_vec):
      # EIS QUE AQUI EU COLOCO A PARADA DA RESTRICAO!!! ACHEI!!!
      if 2*x > y: # <= RESTRICAO
        result+=100.0*((y-x*x)**2.0) + (1-x)**2.0
      else:
        result+=9999999999999999999999999999
    #print list(x_vec), list(y_vec), result
    return result

  def print_status(self, mins,means,vector,txt):
    print txt,mins, means, list(vector)
    #a = 3

class parametros_inversos_dupla_funcao(object):
  def __init__(self, dim=5):
    self.x = None
    self.n = 6*dim
    self.dim = dim
    # ESSE CARA AQUI DOMINIO
    '''
    paramDefCell = {'parameter1', [5000 300000], 10     # RESISITIVDADE AO FLUXO
                'parameter2', [0.90   0.99], 0.01   # POROSIDADE
                'parameter3', [1 4], 0.01           # TORTUOSIDADE
                'parameter4', [10e-6 500e-6], 1e-6  # COMRPIMENTO CARACTERÍSTICO VISCOSO
                'parameter5', [10e-6 500e-6], 1e-6  # COMPRIMENTO CARACTERÍSTICO TÉRMICO
                'parameter6', (64, 64), 0.01 };      # DENSIDADE
    '''

    self.domain = [ (5000, 300000), (0.90, 0.99), (1, 4), (10e-6, 500e-6), (10e-6, 500e-6), (64, 64)]
    print '----------------------------------'
    print self.domain
    print '----------------------------------'
    # ESSE CARA AQUI PEGA O MINIMO PARA SETAR O TAMANHO DA POPULACAO population_size=min(self.n*10,40)
    self.optimizer =  differential_evolution_optimizer(self,population_size=100,n_cross=self.n,cr=1, eps=1e-6, show_progress=True)
    print list(self.x)
    for x in self.x:
      assert abs(x-1.0)<1e-2


  def target(self, vector):
    print 'passei na funcao principal'
    print 'printando os valores'
    print vector

    result=0
    # EIS QUE AQUI EU COLOCO A PARADA DA RESTRICAO!!! ACHEI!!!
    '''
    COMPRIMENTO CARACTERÍSTICO VISCOSO
    COMPRIMENTO CARACTERÍSTICO TÉRMICO
    '''
    #if vector[3] - vector[4] < 0:
    result+=self.objetivo_allard_limp_funcao_dupla(vector)
    #else:
     # result+=9999999999
    #print list(x_vec), list(y_vec), result
    return result

  def print_status(self, mins,means,vector,txt):
    print txt,mins, means, list(vector)
    #a = 3

  def objetivo_allard_limp_funcao_dupla(self, vector):
    params_allard_rigido_simples = self.allard_rigido(0.025,2000,vector[0],vector[1],vector[2],vector[3],vector[4])
    Zs1 = params_allard_rigido_simples['Zs']
    params_allard_rigido_duplo = self.allard_rigido(2*0.025,2000,vector[0],vector[1],vector[2],vector[3],vector[4])
    Zs2 = params_allard_rigido_simples['Zs']

    mat = scipy.io.loadmat('Z_exp.mat')
    
    Z_referencia_esp1 = mat['Z_3']
    F_obj_1 = (numpy.absolute(numpy.subtract((Z_referencia_esp1[500-1:1650-1]/413),(Zs1[500-1:1650-1]/413) ) )**2)
    F_obj_1 = F_obj_1.sum()

    Z_referencia_esp2 = mat['Z_6']
    F_obj_2 = (numpy.absolute(numpy.subtract((Z_referencia_esp2[500-1:1650-1]/413),(Zs2[500-1:1650-1]/413) ) )**2)
    F_obj_2 = F_obj_2.sum()

    F_obj = 3*F_obj_1 + F_obj_2
    return F_obj

  def allard_rigido(self, L, fmax, sigma, phi, alpha_inf, Lambda_material, Lambda_l_material):
    '''
    ## Modelo de material poroso de estrutura rigido
    # Cálcula densidade efetiva dinâmica e compressibilidade efetiva para
    # modelos de fluído equivalente de Johnson-Lafarge. 
    # [Zs alpha] = Lafarge(L,freq,sigma,phi,alpha_inf,Lambda_material,Lambda_l_material)
    # Parâmentros de entrada:
    # L = comprimento da amostra;
    # freq = frequência máxima;
    # sigma = resistividade ao fluxo [Ns/m^2];
    # phi = porosidade [%]
    # alpha_inf = tortuosidade (normalmente varia de 1 à 4)
    # Lambda_material = Comprimento característico viscoso
    # Lambda_l_material = Comprimento característico térmico

    # Como saída tem-se:
    # Zs = impedancia de superfície
    # Alpha = Coeficiente de absorção
    '''

    L_material = L; # espessura do material poroso
    rho_0 = 1.204;      # Densidade [kg/m^3]
    T = 20;               # Temperatura
    c_0 = (331.2+0.6*T);# velocidade de propagação no meio [m/s]
    eta = 1.84e-5;      # Coeficiente de viscosidade, viscosidade absoluta ou viscosidade dinâmica
    gamma = 1.4;        # Razão de calores específicos (gamma=c_p/c_v)
    k_f = 0.026;        # Condutibilidade térmica do fluido [W/mK]
    P_0 = 101320;       # Pressão estática do meio [Pa]
    c_p = 1.0035e3;     # Coeficiente de Calor específico a pressão constante (c_p), (c_v é coef. calor esp. a VOLUME cte)
    Pr = (eta*c_p)/k_f;   # Número de Prandtl (???)--- c_p calor específico, k_f é

    # Parâmetros Macro
    sigma_material = sigma;      # Resistividade ao fluxo [kN/s]
    phi_material =  phi;         # Porosidade 60#
    alpha_inf_material = alpha_inf;    # Tortuosidade
    Lambda = Lambda_material;
    Lambda_l = Lambda_l_material;
    q_0 = eta/sigma;
    # q_l_0 = q_0*Lambda_l/Lambda;
    q_l_0 = (phi_material*Lambda_l*Lambda_l)/8.;

    f = numpy.arange(0, fmax, 1)

    w = 2*math.pi*f ## MODELO DE Johnson - Lafarge

    D = (eta*(phi_material**2)*(Lambda**2));
    C = (4*w*rho_0*(q_0**2)*(alpha_inf_material**2));
    #numpy.divide(a,b,dtype=float)
    G = (1j*w*rho_0*alpha_inf_material*q_0)
    G[G==0] = numpy.finfo(float).eps
    B = numpy.divide((phi_material*eta),G)
    #B = ((phi_material*eta)/(1j*w*rho_0*alpha_inf_material*q_0));
    E = ( 1 + 1j*numpy.divide(C, D));
    E = numpy.power(E,1/2)
    A = 1 + numpy.multiply(B, E)
    rho_rigido = rho_0*alpha_inf_material * A; # densidade dinâmica efetiva

    F = eta*(phi_material**2)*(Lambda_l**2)
    E = numpy.multiply(4*rho_0*Pr*(q_l_0**2), w)
    D =  1 + (1j* numpy.divide(E, F))
    G = (1j*w*rho_0*Pr*q_l_0)
    G[G==0] = numpy.finfo(float).eps
    C = numpy.divide((phi_material*eta), G)
    B = 1 + numpy.multiply(C, numpy.power(D,(1/2)))
    A = numpy.divide(gamma-(gamma-1), B)
    K_ef = numpy.divide(gamma*P_0, A) # módulo de compressibilidade 

    Zc = numpy.sqrt(numpy.multiply(rho_rigido, K_ef))
    A = numpy.divide(rho_rigido, K_ef)
    B = numpy.power(A,(1/2))
    kc_L = numpy.multiply(w, B)
    A = numpy.multiply(kc_L, L_material)
    H = numpy.multiply(phi_material,numpy.arctan(A))
    H[H==0] = numpy.finfo(float).eps
    Zs = -1j*numpy.divide(Zc, H)
    
    alpha = 1 - numpy.power((numpy.absolute(numpy.divide((Zs - rho_0*c_0),(Zs + rho_0*c_0)))), 2);

    params_allard_rigido = {'Zs': Zs, 'alpha': alpha}

    return params_allard_rigido

def run():
  random.seed(0)
  numpy.random.seed(0)
  parametros_inversos_dupla_funcao(1)
  print "OK"


if __name__ == "__main__":
  run()
