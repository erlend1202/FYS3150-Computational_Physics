import numpy as np 
import matplotlib.pyplot as plt 


k_e = 1.38935333 * 10**5
T = 9.64852558 * 10 
V = 9.64852558 * 10**7

B_0 = T
V_0 = 10*V 
d = 10**4
V_0_div_d2 = 9.64852558 


q = 2
m = 40.07

x_0 = 5
z_0 = 10
v_0 = 2

#V_0 = 0.0025*V
#d = 5*10**2

w_0 = q*B_0/m 
w_z2 = 2*q*V_0/(m*d*d)

w_m = (w_0 - np.sqrt(w_0**2 - 2*w_z2))/2
w_p = (w_0 + np.sqrt(w_0**2 - 2*w_z2))/2

print(w_m,w_p)
print(w_z2)

"""
w_m = -2.7
w_p = 0.6


A_p = (v_0 + x_0*w_m)/(w_m - w_p)
A_m = -(v_0 + x_0*w_p)/(w_m - w_p)


def x(t):
    return A_p*np.cos(w_p * t) + A_m*np.cos(w_m * t)

def y(t):
    return -(A_p*np.sin(w_p * t) + A_m*np.sin(w_m * t))

def z(t):
    return z_0*np.cos(np.sqrt(w_z2)*t)

R_m = np.abs(A_p - A_m)
R_p = A_p + A_m



time = np.linspace(0,60,10000)    
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(x(time), y(time), z(time), label='part_1')

plt.show()
circle1 = plt.Circle((0,0), R_m, color='red', fill=False, label="R_-")
circle2 = plt.Circle((0,0), R_p, color='purple', fill=False, label="R_-")
fig, ax = plt.subplots()

ax.add_patch(circle1)
ax.add_patch(circle2)

ax.plot(x(time),y(time), color="lime")
plt.legend()
plt.show()


timesteps = 10000

dt = 0.001

t = np.linspace(0,60,timesteps)
x = np.zeros(timesteps)
y = np.zeros(timesteps)
z = np.zeros(timesteps)

x[0] = 0.2
z[0] = 0.2

r = np.array([x,y,z])

xv = np.zeros(timesteps)
yv = np.zeros(timesteps)
zv = np.zeros(timesteps)

yv[0] = 0.5

v = np.array([xv,yv,zv])

def E(x,y,z):
    return (V_0_div_d2*2*x, V_0_div_d2*2*y, -V_0_div_d2*4*z) 

B = (0,0,B_0)

def cross_prod(A,B):
    return (A[1]*B[2] - A[2]*B[1], A[0]*B[2] - A[2]*B[0], A[0]*B[1] - A[1]*B[0])

def Force(x,y,z,t):
    Electric = E(x,y,z)
    Magnetic = cross_prod((q*v[0][t], q*v[1][t], q*v[2][t]), B)
    return (q*Electric[0] + Magnetic[0], q*Electric[1] + Magnetic[1], q*Electric[2] + Magnetic[2])


def velocity(r,v,t):
    total_force = Force(r[0][t], r[1][t], r[2][t], t)
    v[0][t+1] = total_force[0]/m*dt
    v[1][t+1] = total_force[1]/m*dt
    v[2][t+1] = total_force[2]/m*dt


def posisjon(r,v,t):
    velocity(r,v,t)
    r[0][t+1] = v[0][t+1]*dt
    r[1][t+1] = v[1][t+1]*dt
    r[2][t+1] = v[2][t+1]*dt


for i in range(timesteps-1):
    posisjon(r,v,i)


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(r[0], r[1], r[2], label='test')

plt.show()
"""