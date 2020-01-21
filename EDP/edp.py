import numpy as np 
import matplotlib.pyplot as plt

#variaveis de inicializacao
ageFim  = 50
ageH   = 1 #passo no age
ageNpts = int(ageFim/ageH )
a = 0

tempoFim = 90
tempoH   = 1 #passo no tempo
tempoNpts = int(tempoFim/tempoH)
t = 0

#parametros
s     = 1
d     = 1
betta = 1
delta = 1
alpha = 1
mu    = 1
rho   = 1
c     = 1


#condicoes iniciais
T0 = 100
V0 = 100
I0_t0 = betta*T0*V0
R0_t0 = 1


#variaveis 

I = np.zeros((ageNpts, tempoNpts))
R = np.zeros((ageNpts, tempoNpts))
V = np.zeros((ageNpts, tempoNpts))
T = np.zeros((ageNpts, tempoNpts))

T[0][0] = T0
V[0][0] = V0
R[0][0] = R0_t0
I[0][0] = I0_t0


#condicoes de contorno
I0_t = np.zeros(tempoNpts)
I0_t[t] = betta*T[a][t-1]*V[a][t-1]
R0_t = np.zeros(tempoNpts)
R0_t = 1


print(I)