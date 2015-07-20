from numpy import *
import matplotlib.pyplot as plt

def rho(x,y):
	return sin(pi*x)*sin(pi*y)

N = 2
a, index, epsilon = 1./N,0,10e-5
phi, phi2 = zeros([N+1,N+1],float), zeros([N+1,N+1],float)
phi[:,0],phi[:,N] = 0,1

for i in range(1,N+1):
    phi[0,i],phi[N,i] = phi[0,i-1]+a,phi[N,i-1]+a
while (a > 1./1024):
    x, y = linspace(0,1,N), linspace(0,1,N)
    x2, y2 = ix_(x, y)
    rho2 = rho(x2,y2)
    delt = 1

    while (delt > epsilon):
        phi2[:,0],phi2[:,N] = phi[:,0],phi[:,N]
        phi2[0,:],phi2[N,:] = phi[0,:],phi[N,:]
        for i in range(1,N):
            for j in range(1,N):
                phi2[i,j] = (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])/4 + a**2*pi*rho2[i,j]         
        delt = amax(abs(phi-phi2))
        phi,phi2 = phi2,phi

    fig, ax = plt.subplots()
    im = ax.imshow(phi, vmin=abs(phi).min(), vmax=abs(phi).max(), extent=[0, 1, 0, 1])
    im.set_interpolation('bilinear')
    cb = fig.colorbar(im, ax=ax)
    plt.xlabel('x')
    plt.ylabel('y')
    title = 'Charge Distribution for dx = '+ str(a)
    plt.title(title)
    plt.show()

    a /= 2; N *= 2
    phi2, phi3 = zeros([N+1,N+1],float), zeros([N+1,N+1],float)
    phi3[:,0],phi3[:,N] = 0,1
    for i in range(1,N+1):
        phi3[0,i],phi3[N,i] = phi3[0,i-1]+a,phi3[N,i-1]+a
    
    for i in range(1,N):
        for j in range(1,N):
            phi3[i,j] = phi[i/2,j/2]        
    for i in range(1,N):
        for j in range(1,N):
            if phi3[i,j]==0:
                phi3[i,j] = (phi3[i,j-1]+phi3[i,j+1])/2

    phi,phi3 = phi3,phi
