import numpy as np 
import matplotlib.pyplot as plt


#Inicializacao

#variaveis de inicializacao
ageFim  = 50
deltaA   = 0.1 #passo no age
ageNpts = int(ageFim/deltaA ) + 1
agePt = np.linspace(0, ageFim, ageNpts)
ageCont = 0


tempoFim = 9
deltaT   = 0.01 #passo no tempo
tempoNpts = int(tempoFim/deltaT) + 1 
tempoPt = np.linspace(0, tempoFim, tempoNpts)
tempoCont = 0


#parametros
s     = 1.3*10**5
d     = 0.01
betta = 5*10**-8
c     = 22.30
delta = 0.62
rho   = 8.180
alpha = 30.0
r     = 1.61


k     = 0.80 # coeficiente da funcao exponencial de atraso na exportacao de RNA positivo
tau   = 0.50 # tempo de atraso para a exportacao de RNA positivo
n     = 1.00 # atraso de delta
Rmax  = 50.0 # numero maximo de RNA negativo / complexo de replicacao (Rn)
sigma = 1.30 # taxa de formacao de complexos de replicacao
mu_t  = 0.89 # decaimento natural de Rt
theta = 1.20 # taxa de disponibilidade para traducao
mu_c  = 2.39 # decaimento natural de Rc e Rn


epsilon_s     = 0.998 # efetividade da terapia em diminuir ou bloquear a exportacao de RNA positivo
epsilon_alpha = 0.924 # efetividade da terapia em diminuir ou bloquear a replicacao de RNA positivo
epsilon_r     = 0.290 # efetividade da terapia em diminuir ou bloquear a replicacao de RNA negativo
kappa_t       = 1.000 # fator para aumentar a degradacao de RNA positivo disponivel para traducao
kappa_c       = 1.000 # fator para aumentar a degradacao de RNA positivo e negativo no complexo de replicacao

tol  = 10**-4 # tolerancia

#variaveis 

I = np.zeros((tempoNpts, ageNpts))
Rp = np.zeros((tempoNpts, ageNpts))
Rt = np.zeros((tempoNpts, ageNpts))
Rn = np.zeros((tempoNpts, ageNpts))
V = np.zeros(tempoNpts)
T = np.zeros(tempoNpts)

#condicoes iniciais
T0 = 1.3*10**5
V0 = 5.45*10**5
I0_t0 = betta*T0*V0
Rt0_t0 = 1
Rp0_t0 = 0
Rn0_t0 = 0


T[0] = T0
V[0] = V0

Rt[0][0] = Rt0_t0
Rp[0][0] = Rp0_t0
Rn[0][0] = Rn0_t0
I[0][0] = I0_t0

rho1 = 0.00

# Condicao inicial para as EDPs
for ageCont in range(1, ageNpts):
    I[0][ageCont] = (betta*T[0]*V[0]*np.exp(-deltaA*ageCont))*deltaA + I[0][ageCont-1]
    
    if ageCont*deltaA < tau:
        rho1 = 0
    else:
        rho1 = (1 - np.exp(-k*((ageCont*deltaA) - tau)))*rho
    
    Rn[0][ageCont] = (r*Rp[0][ageCont-1] - r*Rp[0][ageCont-1]*(Rn[0][ageCont-1]/Rmax) 
    - mu_c*Rn[0][ageCont-1])*deltaA + Rn[0][ageCont-1]

    Rp[0][ageCont] = (alpha*Rn[0][ageCont-1] + sigma*Rt[0][ageCont-1] 
    - theta*Rp[0][ageCont-1]- rho1*Rp[0][ageCont-1]
    - mu_c*Rp[0][ageCont-1])*deltaA + Rp[0][ageCont-1]

    Rt[0][ageCont] = (theta*Rp[0][ageCont-1] - sigma*Rt[0][ageCont-1] 
    - rho1*Rt[0][ageCont-1] - mu_t*Rt[0][ageCont-1])*deltaA + Rt[0][ageCont-1]
    
    

#Fim da inicializacao

def calcIntegral(I,Rp,Rt):
    soma = 0.0
    for a in range(0, ageNpts):
        soma = soma + I[a]*Rp[a]*Rt[a]
    return soma/float(ageNpts)
#duvida: a calcIntegral2 so eh usada para calcular N, mas estava comentado no HCV_model.cpp
    #tiro essa funcao ou deixo??
def calcIntegral2(v1, v2):
    soma = 0.0
    for a in range(0, ageNpts):
        soma = soma + v1[a]*v2[a]*np.exp(-delta*a)
    return soma/float(ageNpts)    


#Solve
#duvida: Eh para comecar em t = 1 ou t = 0
for t in range(1, tempoNpts - 1):
    
    T[t] = (s - d*T[t-1] - betta*V[t-1]*T[t-1])*deltaT + T[t-1]
    # duvida: Na integral os parametros sao do tempo anterior(I[t-1]...) ou I[t] ??
    V[t] =  deltaT*((1 - epsilon_s)*rho*calcIntegral(I[t-1],Rp[t-1],Rt[t-1]) 
    - c*V[t-1]) + V[t-1]
    
    rho1 = 0.0
    
    Rp[t][0] = 0
    Rn[t][0] = 0
    Rt[t][0] = 1
    I[t][0]  = betta*V[t]*T[t] #duvida: O lado da direita eh V[t] ou V[t-1]?
    #duvida: por que a < AGE-1?? a ultima posicao, Rp[AGE-1], nao eh utilizada nunca?
    #equivalente no HCV.cpp: linha 259
    for a in range(1, ageNpts - 1):
        #duvida: rho fica variavel?
        if(float(a*deltaA) < tau):
            rho1 = 0
        else:
            rho1 = rho*(1 - np.exp(-k*(float(a*deltaA) - tau)))
        #duvida: utilizando a equacao com delta fixo //equivalente: linha 280
        I[t][a] = (-delta*I[t-1][a] - (I[t-1][a] - I[t-1][a-1])/deltaA)*deltaT + I[t-1][a]
        
        Rn[t][a] = ((1 - epsilon_r)*r*Rp[t-1][a]*(1 - (Rn[t-1][a]/Rmax)) - kappa_c*mu_c*Rn[t-1][a]
                     - (Rn[t-1][a] - Rn[t-1][a-1])/(deltaA ))*deltaT + Rn[t-1][a]
        
        Rp[t][a] = ((1 - epsilon_alpha)*alpha*Rn[t-1][a] + sigma*Rt[t-1][a] - theta*Rp[t-1][a]
		            - (1- epsilon_s)*rho1*Rp[t-1][a] - kappa_c*mu_c*Rp[t-1][a]
                    - (Rp[t-1][a] - Rp[t-1][a-1])/(deltaA))*deltaT + Rp[t-1][a]
        
        Rt[t][a] = (theta*Rp[t-1][a] - sigma*Rt[t-1][a] - (1 - epsilon_s)*rho1*Rt[t-1][a] - kappa_t*mu_t*Rt[t-1][a]
                    - (Rt[t-1][a] - Rt[t-1][a-1])/(deltaA))*deltaT + Rt[t-1][a]
        

#fazer os plots dos graficos
















    
    
    