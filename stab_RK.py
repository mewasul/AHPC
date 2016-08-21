
import matplotlib.pyplot as plt
import numpy as np

# stability region

# set up "complex plane"
x = np.linspace(-4,4,100)
y = np.linspace(-4,4,100)

[X, Y] = np.meshgrid(x, y)
Z = X + 1j*Y

# calculate region: 
U = abs(1 + Z + 0.5*Z**2 + (1./6)*Z**3 + (1./24)*Z**4)

# and plot
plt.contour(X, Y, U, [1])             # contour for U = 1
plt.contourf(X, Y, U, [0, 1])         # colour area in interval [0,1]

plt.xlabel("$\lambda_{Re} \Delta t$")
plt.ylabel("$\lambda_{Im} \Delta t$")

plt.show()

# calculate amplification error + phase error
a = lambda z: (1 + 0.5*z**2 + (1./24)*z**4)**2 + (z - (1./6)*z**3)**2
theta = lambda z: np.tan((z - (1./6)*z**3)/(1 + 0.5*z**2 + (1./24)*z**4)) - z

z = np.linspace(-3,3,100)             # stable around approx +-2.8

plt.plot(z, a(z), z, theta(z))
plt.xlabel("$\lambda_{Im} \Delta t$")
plt.ylabel("Error")
plt.show()
