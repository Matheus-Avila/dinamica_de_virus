import scipy.integrate as integrate
import numpy as np
#definindo as variaveis // VALORES DO HCV_MODEL.cpp

alpha = 30 #taxa de producao de vRNA
beta = 5*10**-8 #taxa de infecção
rho = 1  #taxa de montagem de vRNA
mu = 0.89 #taxa de degradacao de vRNA // TROCAR ESSE VALOR mu_t
delta = 0.62 #taxa de morte de T
s = 130000 #taxa de produção de T
d = 0.01 #taxa de morte de T
c = 22.30 #taxa de morte de V

#Assumindo t = 0

t0 = 0
R0_t = 1
R0_a = 0

#Definindo as EDOs

T_t = np.zeros([20])
I_t = np.zeros([20])
I_a = np.zeros([20])
R_t = np.zeros([20])
R_a = np.zeros([20])
V_t = np.zeros([20])

#Target 
T_t[i] = s - beta*V_t*T_t[i-1] - d*T_t[i-1]

#Infectadas
I_t = -I_a - delta*I_t[i-1]

I_a#????

#vRNA intracelular
R_t#????

R_a = -R_t + alpha - (rho - mu)*R_a[i-1]

#Virus
V_t = integrate.quad(rho*R_a*I_t) - c*V_t[i-1]

#2 opcao: HCV_MODEL
T = (s - d*T - beta*V*T)*deltaT + T
V = ((1-epsilon_s)*rho*calcIntegral(I,Rp,Rt) - c*V)*deltaT + V

