from Question1 import omg, M, ka_mat, ks_mat
import numpy as np
import matplotlib.pyplot as plt
import math

b   = 6                 # halfspan
c   = 2                 # chord
xf  = 0.5 * c           # flexural axis
e   = (xf/c)-0.25       # eccentricity
cla = 2 * math.pi       # lift-curve slope
rho = 1.225             # density of air

w   = (omg[0][0] + omg[0][2])/2
k1  = (w*c)/(2*10)
k2  = (w*c)/(2*80)
i   = 1j

komg = []
kzet = []
kvel = []

for k in np.arange(k2, k1, 0.1):

    mt = -((5)/(2+(5*k)))

    # Aerodynamic Damping (C)
    c11     = c*b/10*cla
    c12     = 0
    c21     = -c**2*b/8*e*cla
    c22     = -c**3*b/24*mt
    C       = np.matrix([[c11, c12], [c21, c22]])

    A   = M
    B   = i*rho*((c/2)/k)*C
    C   = rho*(((c/2)/k)**2)*ka_mat
    
    F = A - B - C
    Q = np.linalg.inv(ks_mat)*F

    eigv = np.linalg.eigvals(Q)

    tmpkOmg = [] #temporary arrays, ignore pls
    tmpkZet = []
    for x in range(len(eigv)): 
        omgeq = 1/np.sqrt(eigv[x].real)
        zeteq = ((eigv[x].imag/eigv[x].real)/2)
        tmpkOmg.append(omgeq)
        tmpkZet.append(zeteq)
    
    komg.append(tmpkOmg)
    kzet.append(tmpkZet)

    w = [(1/np.sqrt(i.real))for i in eigv]
    kvel.append([(j*c)/(2*k) for j in w])
    
# Graph Plot
plt.subplot(2,1,1)
plt.plot(kvel, komg) 
plt.xlabel("V in m/s")
plt.ylabel("ω in rad/s")
plt.grid()
plt.title("ω vs V")

plt.subplot(2,1,2)
plt.plot(kvel, kzet)
plt.xlabel("V in m/s")
plt.ylabel("ζ")
plt.grid()
plt.title("ζ vs V")
plt.show()
