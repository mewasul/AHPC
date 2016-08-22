
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class ShallowWater:

    def __init__(self,xmin, xmax, ymin, ymax, tmin, tmax, dx, dy, dt):
        """
            Doc
        """
        self.dx, self.dy, self.dt = dx, dy, dt

        Nx = int((xmax-xmin)/dx) + 1
        Ny = int((ymax-ymin)/dy) + 1
        Nt = int((tmax-tmin)/dt) + 1

        x = np.linspace(xmin,xmax,Nx)
        y = np.linspace(xmin,xmax,Nx)

        [X, Y] = np.meshgrid(x,y)

        vx = np.zeros((Nt, Nx, Ny))
        vy = np.zeros((Nt, Nx, Ny))
        h  = np.zeros((Nt, Nx, Ny))

        self.Nx, self.Ny, self.Nt, self.X, self.Y, self.vx, self.vy, self.h = \
             Nx, Ny, Nt, X, Y, vx, vy, h


    def du(self,u,dim,h):
        """ 
            Gives spacial derivative in one dimension using 7 pt stencil
            assuming periodic boundary conditions (using modulo operation)

            Args:
                Array u of spacial points (of iteration t)
                dim = 0 or 1; in x or y direction
                h = spacial step (typical dx or dy)
            Returns:
                Array u' of spacial points of first derivative approximation
        """
      
        du = np.zeros(np.shape(u))
        coeff = (-1./60, 3./20, -3./4, 0, 3./4, -3./20, 1./60)

        # use shift operator to calculate derivative
        for i in range(7):
            du += (1./h)*coeff[i]*np.roll(u, i-3,axis=dim)

        return du
                         
    def calculate(self,ic_vx, ic_vy, ic_h):
        """
            Runs calculations using initializations from init; time derivatives are
            calculated using Runge-Kutta of 4th order, spacial derivatives the 7 pt stencil
            
            Args:
                ic - initial function, giving position of all points at time 0
            Returns:
                Calculated value of all points at all time steps (of size Nx x Ny x Nt)
        """
        Nx, Ny, Nt, dx, dy, X, Y, vx, vy, h = \
             self.Nx, self.Ny, self.Nt, self.dx, self.dy, self.X, self.Y, self.vx, self.vy, self.h

        vx[0] = ic_vx(X, Y)       
        vy[0] = ic_vy(X, Y)
        h[0] = ic_h(X, Y)

        # save memory by not saving all of the k's (save as we go)
        k_vx = np.zeros(np.shape(vx[0]))
        k_vy = np.zeros(np.shape(vy[0]))
        k_h  = np.zeros(np.shape(h[0]))

        vx_tmp = np.zeros(np.shape(vx[0]))
        vy_tmp = np.zeros(np.shape(vy[0]))
        h_tmp  = np.zeros(np.shape(h[0]))

        coeff1 = (dt/6, dt/3, dt/3, dt/6)
        coeff2 = (0, dt/2, dt/2, dt)

        Nt2 = Nt
        for m in range(1,Nt2):     
             print m 
             # first step of RK4: u(n+1) = u(n) + ...
             vx[m] = vx[m-1].copy()
             vy[m] = vx[m-1].copy()
             h[m]  = h[m-1].copy()

             # second step: ... + h/6*(k1 + 2k2 + 2k3 + k4)
             for c in range(4):       # calculate k 4 times

                 # calculate y to use in f(u,t); this is values
                 # u, u + h/2*k1; u + h/2*k2; u + h*k3 
                 vx_tmp = vx[m-1] + coeff2[c]*k_vx
                 vy_tmp = vy[m-1] + coeff2[c]*k_vy
                 h_tmp  = h[m-1] + coeff2[c]*k_h

                 dvxdx = self.du(vx_tmp,0,dx)
                 dvxdy = self.du(vx_tmp,1,dy)
                 dhdx  = self.du(h_tmp,0,dx)

                # print self.maximum(vx_tmp*dvxdx)
                # print self.maximum(vy_tmp*dvxdy)
                # print self.maximum(g*dhdx)
                 
                 k_vx = -(vx_tmp*dvxdx + vy_tmp*dvxdy + g*dhdx)

                 dvydx = self.du(vy_tmp,0,dx)
                 dvydy = self.du(vy_tmp,1,dy)
                 dhdy  = self.du(h_tmp,1,dy)

                 k_vy = -(vx_tmp*dvydx + vy_tmp*dvydy + g*dhdy)
                
                 k_h = -(dvxdx*h_tmp + vx_tmp*dhdx + dvydy*h_tmp + vy_tmp*dhdy)
                 
                 # add values to next time step from this k (iteration c)
                 vx[m] += coeff1[c]*k_vx
                 vy[m] += coeff1[c]*k_vy
                 h[m]  += coeff1[c]*k_h

             print coeff1
             print self.maximum(h[m])
             if m>20:
                 self.plot(vx[:m])
                 #self.plot(vy[:m])
                 self.plot(h[:m])

        return vx[:Nt2], vy[:Nt2], h[:Nt2]

    def maximum(self,y):
        return max([max([y[i][j] for j in range(len(y[i]))]) for i in range(len(y))])

    def plot(self,y):
        """
            Plots y over domain given by [X, Y] for all given time steps
        """
        X, Y, dt = self.X, self.Y, self.dt

        Nt = len(y)       # first dim. = time
            
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        surf = ax.plot_surface(X, Y, y[-1])
        plt.show()

            
xmin = ymin = tmin = 0
xmax = ymax = 100
tmax = 10
dx = dy = 0.5
dt = 0.5*dx*dx
g = 9.81

# ic_s is not really needed but it makes the code setup consistent 
ic_h = lambda x, y : np.exp(-(1./25)*((x-30)**2+(y-30)**2))
ic_s = lambda x, y : 0                                          

sw = ShallowWater(xmin,xmax,ymin,ymax,tmin,tmax,dx,dy,dt)
vx, vy, h = sw.calculate(ic_s, ic_s, ic_h)

sw.plot(h)
