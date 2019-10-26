import numpy as np
import matplotlib.pyplot as plt
import time

from scipy.integrate import quad


def ode_FE(f, U_0, dt, T):
    N_t = int(round(float(T)/dt))
    # Ensure that any list/tuple returned from f_ is wrapped as array
    f_ = lambda u, t: np.asarray(f(u, t))
    u = np.zeros((4,N_t+1))
    t = np.linspace(0, N_t * dt, N_t+1) # conferir N+1
    u = U_0[:,:22]
    print("t = ", t)
    print('u = ', u);
    for n in range(N_t):
        u[:,n] = u[:,n-1] + dt*f_(u[:,n-1], t[n-1])# passa uma coluna por iteraçao
    return u, t

def rhs(u, t):
    N = len(u) - 1
    rhs = np.zeros((4,N+1))
    print(len(u))
    print(N)
    print("rhs ",len(rhs))
    for i in range(1, N): # conferir os dois pontos
        rhs[0,i] = -(1/2*dx)*(u[i+1] - 2*u[i] + u[i-1]) - delta*u[i] # I
        rhs[1,i] = -(1/2*dx)*(u[i+1] - 2*u[i] + u[i-1]) + alpha - rho*u[i] -mu*u[i] # R
        
    rhs[2,:] = rho*rhs[0,:]*rhs[1,:] - c*rhs[2,:]# equacao dos virus
    rhs[3,:] = s - d*rhs[3,:] - beta*rhs[2,:]*rhs[3,:]# equacao das celulas alvo
    return rhs

def demo():
    global beta, V, Ta, dx, L, x  # needed in rhs
    # Ta e V nao sao usados
    global rho, mu, delta, s, d, c, omega, alpha
    L = 50 # age final
    N = 201 # numero de pontos no age
    x = np.linspace(0, L, N) # Mudei para N (estava N+1)linha 63 dava erro
    dx = x[1] - x[0]
    dt = 0.1
    # parametros
    alpha = 30  # taxa de producao de vRNA
    beta = 5 * 10 ** -8  # taxa de infecção
    rho = 8.18  # taxa de montagem de vRNA
    mu = 0.89  # taxa de degradacao de vRNA // TROCAR ESSE VALOR mu_t
    delta = 0.62  # taxa de morte de T
    s = 130000  # taxa de produção de T
    d = 0.01  # taxa de morte de T
    c = 22.30  # taxa de morte de V


    # condicoes iniciais e condicoes de contorno
    U_0 = np.zeros((4,N+1))
    
    U_0[0,0] = 1 #condicao inicial do RNA
    def RNA0(x):
        return (alpha/(rho+mu))+(1-(alpha/(rho+mu)))*np.exp(-(rho+mu)*x) # Arrumei a equacao
    U_0[0,1:] = RNA0(x) # condicao de contorno do RNA (pegar do artigo Eq. (9))
    # print("U_0[0,:]",U_0[0,:])
    #criar funcao pra jogar lixo fora
    [NN,lixo] = quad(lambda x: -delta*dx, 0, 50)
    plt.plot(x, U_0[0,1:])
    plt.show()
    omega = np.exp(NN) # omega deve ser pequeno
    [NNN,lixo] = quad(lambda x: RNA0(x), 0, 50)# Integral correta wolfram alpha
    print("NNN: ",NNN)
    print("omega: ", omega, " //Devia ser bem menor")
    print("NN: ", NN, "//Esperado: -31") # Diferente do wolfram alpha(-31)
    N_v = float(rho)*float(NNN)*float(omega) # Isso ta certo??
    print("N_v: ",N_v)
    V0 = (beta * N_v * s - d * c) / (beta * c) # beta muito pequeno
    print('V0 = ', V0)
    U_0[1,0] = V0 # condicao inicial do virus (V)
    Ta0 = c / (beta * N_v)
    U_0[2,0] = Ta0 # condicao inicial das celulas alvo (Target)
    U_0[3,0] = beta*V0*Ta0 # condicao inicial da infectada
    U_0[3,1:] = beta*V0*Ta0*np.exp(-delta*x) # condicao de contorno da infectada

    t0 = time.perf_counter()
    #from ode_system_FE import ode_FE
    u, t = ode_FE(rhs, U_0, dt, T=2)
    t1 = time.clock()
    print('CPU time: %.1fs' % (t1 - t0))

    # Make movie
    #import os
    #os.system('rm tmp_*.png')
    #plt.ion()
    y = u[0,:]
    plt.plot(x, y)
    plt.axis([x[0], x[-1]])
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.show()
    counter = 0
    # Plot each of the first 100 frames, then increase speed by 10x
    #change_speed = 100
    #for i in range(0, u.shape[0]):
    #    print(t[i])
    #    plot = True if i <= change_speed else i % 10 == 0
    #    lines[0].set_ydata(u[i,:])
    #    if i > change_speed:
    #        plt.legend(['t=%.0f 10x' % t[i]])
    #    else:
    #        plt.legend(['t=%.0f' % t[i]])
    #    plt.draw()
    #    if plot:
    #       plt.savefig('tmp_%04d.png' % counter)
    #        counter += 1
        #time.sleep(0.2)


if __name__ == '__main__':
    demo()