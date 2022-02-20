import numpy as np
import math
import matplotlib.pyplot as plt

# --------- Parameters ---------
b   = 6                 # halfspan
c   = 2                 # chord
xf  = 0.5 * c           # flexural axis
e   = (xf/c)-0.25       # eccentricity
xa  = 0.5 * c           # mass axis
ma  = 100               # mass/area
mt  = -1.2              # non-dimensional pitch damping derivative
EI  = 1 * (10 ** 7)     # bending rigidity
GJ  = 1 * (10 ** 6)     # torsional rigidity
cla = 2 * math.pi       # lift-curve slope
rho = 1.225             # density of air

# -------- Arrays ----------
omg = []
zet = []
vel = []

for v in range(0, 90):
    
    vel.append(v)
    
    # Structural Mass Matrix (M)
    m11     = (b*c) / 5
    m12     = (b/4) * (((c**2)/2)-(c*xf))
    m21     =  m12
    m22     = (b/3) * (((c**3)/3)-((c**2)*xf)+((xf**2)*c))
    m_mat   = np.matrix([[m11, m12], [m21, m22]]) # Mass
    M       = m_mat*ma
    
    # Aerodynamic Damping (C)
    c11     = c*b/10*cla
    c12     = 0
    c21     = -c**2*b/8*e*cla
    c22     = -c**3*b/24*mt
    c_mat   = np.matrix([[c11, c12], [c21, c22]]) # Damping
    C       = c_mat*rho*v
    
    # Aerodynamic Stiffness Matrix (ka)
    ka11    = 0
    ka12    = (c*b*cla)/8
    ka21    = 0
    ka22    = -((c**2)*b*cla*e)/6
    ka_mat  = np.matrix([[ka11, ka12], [ka21, ka22]])

    # Structural Stiffness Matrix (ks)
    ks11    = (4*EI)/(b**3)
    ks12    = 0
    ks21    = 0
    ks22    = GJ/b
    ks_mat  = np.matrix([[ks11, ks12], [ks21, ks22]])

    K       = (ka_mat*rho*(v**2)) + ks_mat # Stiffness 

    # Q matrices
    q11     = np.matrix([[0.0,0.0],[0.0,0.0]]) # 2x2 zero matrix
    q12     = np.matrix([[1.0,0.0],[0.0,1.0]]) # identity matrix
    q21     = -np.linalg.inv(M)*K # Stiffness element
    q22     = -np.linalg.inv(M)*C # Damping element

    # tolist() returns the first order matrices above, into a NESTED LIST
    q11     = q11.tolist()
    q12     = q12.tolist()
    q21     = q21.tolist() 
    q22     = q22.tolist()

    for i in range (4): 
        if i < 2:
            q11.append(q21[i])
            for j in range(2):
                q11[i].append(q12[i][j])
        else:
            for j in range(2):
                q11[i].append(q22[i-2][j-2])

    q11 = np.matrix(q11)

    # Eigenvalues
    eigv    = np.linalg.eigvals(q11) # linalg.eigvals() returns matrix q11 into eigenvalues
    n       = len(eigv)
    
    # Omega and Zeta calculations
    tmpOmg = [] #temporary arrays, ignore pls
    tmpZet = []
    for x in range(n): 
        calc = np.sqrt((eigv[x].real**2)+(eigv[x].imag**2))
        tmpOmg.append(calc)
        tmpZet.append((-eigv[x].real)/calc)
    
    omg.append(tmpOmg)
    zet.append(tmpZet)
    

# -----------CALCULATION ENDS HERE---------------

# Graph Plot
plt.subplot(2,1,1)
plt.plot(vel, omg) 
plt.xlabel("V in m/s")
plt.ylabel("ω in rad/s")
plt.grid()
plt.title("ω vs V")

plt.subplot(2,1,2)
plt.plot(vel, zet)
plt.xlabel("V in m/s")
plt.ylabel("ζ")
plt.grid()
plt.title("ζ vs V")
plt.show()
