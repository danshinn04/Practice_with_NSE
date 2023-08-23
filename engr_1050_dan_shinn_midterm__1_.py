# -*- coding: utf-8 -*-
"""ENGR 1050 - Dan Shinn - Midterm

# Imperfect NSE Simulator
"""

#ReadME
#Main is in seond section but run this so we initialize the navier stokes simulation "class": we didn't go over this in class but I guess this is my new addition to the project.
#Reference 1: http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf
#Reference 2: https://mikeash.com/pyblog/fluid-simulation-for-dummies.html 
#Reference 3: (Goes down in order of importance): https://en.wikipedia.org/wiki/Immersed_boundary_method)
#Reference 4: https://www.youtube.com/watch?v=JBmS--3L2eQ&ab_channel=MachineLearning%26Simulation, https://www.youtube.com/watch?v=BQLvNLgMTQE&ab_channel=MachineLearning%26Simulation 
#Reference 5: https://www.youtube.com/watch?v=alhpH6ECFvQ&ab_channel=TheCodingTrain 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
N = 50
iter = 2

"""
Solves the incompressible Navier Stokes equations for a 2D simulation
using small dt approximation, and some ideas come from Chorin's Projection.
NSE IBM: ∂u/∂t + (u ⋅ ∇) u = − 1/ρ ∇p + ν ∇²u + f
Incompressibility:  ∇ ⋅ u = 0 => Unbounded and since I do a dt simulation, I basically run some linear algebra functions referenced from 2 and iterate it multiple times to approach 
∇ ⋅ u = 0
"""

def IX(x, y): #Finds the indexes for the density array when you enter position x, y of the grid.
    x = np.clip(x,0,N-1)
    y = np.clip(y,0,N-1) #Basically I kept getting out of bound error because for linear algebra solving function, since it considers density points near the point of calculation, 
    #when you were doing a lin alg computation at the edge of the simulation field, simulation goes out of bound. So, another source of error, minor density propagation errors at edges.
    return x+(y*N) 

class Fluid: #Initializes the class of density array N by N to calculate how density will propagate over time.
    def __init__(self, dt, diffusion, viscosity):
        self.size = N
        self.dt = dt
        self.diff = diffusion
        self.visc = viscosity 
        self.s = np.zeros(N*N)
        self.density = np.zeros(N*N)
        self.Vx = np.zeros(N*N)
        self.Vy = np.zeros(N*N)
        self.Vx0 = np.zeros(N*N)
        self.Vy0 = np.zeros(N*N)
    def step(self): # Step as in Time Step, to simulate after 0.1sec, what happend to the density gradient.
        visc     = self.visc
        diff     = self.diff
        dt       = self.dt
        Vx      = self.Vx
        Vy      = self.Vy
        Vx0     = self.Vx0
        Vy0     = self.Vy0
        s       = self.s
        density = self.density
        self.diffuse(1, Vx0, Vx, visc, dt) #Computation of "diffusion" of velocity of Vx. Compares between Vxprev and Vx because we cannot create more velocity out of thin air. 
        #Can be thought of as transfer of momentum based on random movement also considering the density
        self.diffuse(2, Vy0, Vy, visc, dt)
        #An analogy is pepperflakes and soap dispersion effect.
        self.project(Vx0, Vy0, Vx, Vy)
        #You cannot create fluid out of thin air so this function makes sure that happens.
        self.advect(1, Vx, Vx0, Vx0, Vy0, dt) #Kind of confusing because it seems similar to diffuse term but advection of velocity is bulk movement of fluid so it is not like due to random
        #movement of individual velocities leading to diffusion of velocity. 
        self.advect(2, Vy, Vy0, Vx0, Vy0, dt)
        self.project(Vx, Vy, Vx0, Vy0)
        self.diffuse(0, s, density, diff, dt)
        self.advect(0, density, s, Vx, Vy, dt) 
    def addDensity(self, x, y, amount): #Adds density to the certain point of density array.
        index = IX(x,y)
        self.density[index] += amount
    def addVelocity(self, x, y, amountX, amountY): #Add velocity to the density at certain point.
        index = IX(x,y)
        self.Vx[index] += amountX
        self.Vy[index] += amountY
    def renderD(self): #Basically we plot the Density matrix and figure out how it changed after small dt.
        dmap = np.zeros(N*N)
        for i in range(N):
            for j in range(N):
                dmap[IX(i,j)] = self.density[IX(i,j)]
        dmap = dmap.reshape(N,N)
        fig, ax = plt.subplots()
        im = ax.imshow(dmap, cmap='jet')
        fig.colorbar(im)

        # Show the plot
        plt.show()
    def diffuse(self, b, x, x0, diff, dt):
        a = dt*diff*(N-2)*(N-2)
        self.lin_Solve(b,x,x0,a,1+6*a ) #These equations (Like the 1+6*a lde) were referenced from paper 2 and video 4.

    def lin_Solve(self, b, x, x0, a, c):
        cRecip = 1.0/c
        for k in range(iter): #Runs the LDE solver multiple times, the higher the iter, the greater the accuracy but 2 should be accurate enough
            for j in range(1,N-1):
                for i in range(1,N-1):
                    x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i+1,j)] + x[IX(i-1,j)]+x[IX(i,j+1)]+x[IX(i,j-1)]))*cRecip
            self.set_bnd(b,x) #The method is pretty popular, and lin solve equations models are written by this famous University of Toronto guy that 
            # A lot of video games and fluid simulators refer to his equations to write their simulations so it is reputable and accurate.
    def project(self, velocX, velocY, p, div): #Explained above but in a nutshell, makes sure conservation of mass is held so like random 1.0 * 10^9 m/s doesn't happen by
        for j in range(1, N - 1):
            for i in range(1, N - 1):
                div[IX(i, j)] = -0.5*(  #Comparing the two indexes above (Yup, Ydown) and (Xright, Xleft) and averaging its velocity as a temp variable associated to that index
                         velocX[IX(i+1, j  )]
                        -velocX[IX(i-1, j ) ]
                        +velocY[IX(i  , j+1)]
                        -velocY[IX(i  , j-1)]
                    )/N;
                p[IX(i, j)] = 0
        self.set_bnd(0, div) #And then doing a boundary solve to make sure it's not next to wall
        self.set_bnd(0, p)
        self.lin_Solve(0, p, div, 1, 6) #And then linearly compute the velocities.
        for j in range(1, N - 1):
            for i in range(1, N - 1):
                velocX[IX(i, j)] -= 0.5 * (  p[IX(i+1, j)]
                                                -p[IX(i-1, j)]) * N #double check
                velocY[IX(i, j)] -= 0.5 * (  p[IX(i, j+1)]
                                                -p[IX(i, j-1)]) * N
        self.set_bnd(1, velocX)
        self.set_bnd(2, velocY)

    def advect(self, b, d, d0, velocX, velocY, dt):
        i0, i1, j0, j1 = 0, 0, 0, 0
        
        dtx = dt * (N - 2) #Basically except for the two edges, we are considering the bulk movement.
        dty = dt * (N - 2)

        s0, s1, t0, t1 = 0, 0, 0, 0
        tmp1, tmp2, x, y = 0, 0, 0, 0
        Nfloat = N
        ifloat, jfloat = 0, 0
        i, j = 0, 0
        
        for j in range(1, N - 1):
            jfloat = float(j)
            for i in range(1, N - 1):
                ifloat = float(i)
                
                tmp1 = dtx * velocX[IX(i, j)]
                tmp2 = dty * velocY[IX(i, j)] 

                x = ifloat - tmp1
                y = jfloat - tmp2

                if x < 0.5:
                    x = 0.5
                if x > Nfloat + 0.5:
                    x = Nfloat + 0.5
                    
                i0 = math.floor(x)
                i1 = i0 + 1.0
                if y < 0.5:
                    y = 0.5
                if y > Nfloat + 0.5:
                    y = Nfloat + 0.5
                    
                j0 = math.floor(y)
                j1 = j0 + 1.0

                s1 = x - i0
                s0 = 1.0 - s1
                t1 = y - j0
                t0 = 1.0 - t1

                i0i = int(i0)
                i1i = int(i1)
                j0i = int(j0)
                j1i = int(j1)#referred to reference2 code for linear algebra equation and converted the C code snippets into python and changed some things so it runs.

                d[IX(i, j)] = (
                    s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)])
                    + s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)])
                )
        
        self.set_bnd(b, d)

    def set_bnd(self, b, x): 
#It's basically like a mirror existing at the wall to make sure the velocity at the edge is reflected back into the N*N array so that nothing is lost from the simulation.
        for k in range(1, N-1):
            for i in range(1, N-1):
                x[IX(i, 0)] = -x[IX(i, 1)] if b == 2 else x[IX(i, 1)]
                x[IX(i, N-1)] = -x[IX(i, N-2)] if b == 2 else x[IX(i, N-2)]
        
        for k in range(1, N-1):
            for j in range(1, N-1):
                x[IX(0, j)] = -x[IX(1, j)] if b == 1 else x[IX(1, j)]
                x[IX(N-1, j)] = -x[IX(N-2, j)] if b == 1 else x[IX(N-2, j)]
        
        x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)])
        x[IX(0, N-1)] = 0.5 * (x[IX(1, N-1)] + x[IX(0, N-2)])
        x[IX(N-1, 0)] = 0.5 * (x[IX(N-2, 0)] + x[IX(N-1, 1)])
        x[IX(N-1, N-1)] = 0.5 * (x[IX(N-2, N-1)] + x[IX(N-1, N-2)])

"""# 2D Boundary Generator"""

# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 01:35:56 2023

@author: Dan
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import math

xlist = np.random.randint(50, size = 1000)
ylist = np.random.randint(50, size = 1000) #Since I don't have training data in 2D for this, I will create one!
vxlist = np.random.randint(50, size = 1000) 
vylist = np.random.randint(50, size = 1000)

def findAngle(reference, vector):
    x1, y1 = reference
    x2, y2 = vector
    dot = x1*x2 + y1*y2      # dot product
    det = x1*y2 - y1*x2      # determinant
    angle = math.atan2(det, dot)  
    angle = math.degrees(angle) #Like we did in lecture, returns angle but I used determinants because it looks cleaner.
    
    return angle

def ListWithAngleGenerator(points, centroid): #Updates list to have angle as well now.
    referencevector = [1, 0]
    listwithangle = []
    for point in points:
        xy = [point[0], point[1]]
        vectorx = point[0]-centroid[0]
        vectory = point[1]-centroid[1]
        vector = [vectorx, vectory]
        if vectorx == 0 and vectory == 0: #For all the points it needs an angle appended to it to work properly.
            angle = 0
        else:
            angle = findAngle(referencevector, vector)
        xy.append(angle+180)
        listwithangle.append(xy)
    return listwithangle    
    
def findCentroid(points): #Finds the Centroid of the data set at which we will compare angles to find boundary
    sumx = 0
    sumy = 0
    N = len(points)
    for point in points:
        sumx = sumx + point[0]
        sumy = sumy + point[1]
    Centroid = (sumx/N, sumy/N) #basically the mean point
    return Centroid

def distancecalculator(pointA, pointB): #Uses Pythagorean Theorem to calculate distance between likely reference point and another point.
    x1, y1 = pointA
    x2, y2 = pointB
    D = ((x2-x1)**2+(y2-y1)**2)**(1/2)
    return D

def indexesOfPointsBoundary(listwithangle, centroid, anglelist, incrementanglecount):
    #Basically each list points have an angle associated to it from the centroid point. This algorithm compares the distances between all the points of similar angle and 
    #spits out the one with the longest length at that range. This is how boundary is calculated.
    longestpointindexes = []
    for angle in anglelist:
        length = 0
        lengthlongest = 0
        index2 = 0
        indexlongest = 0 #basically updates index as you compare and find something bigger.
        for point in listwithangle:
            theta = point[2]
            if theta < angle + incrementanglecount and theta > angle:
                length = distancecalculator(centroid, (point[0], point[1]))
                if lengthlongest < length: 
                    lengthlongest = length
                    indexlongest = index2
            index2 = index2 + 1
        longestpointindexes.append((indexlongest))
    return longestpointindexes
        
def plotter(xList, yList, Indexes_to_Highlight): #Plots the initial figure and the boundaries are highlighted in red.
    index = 0
    RedscatterlistX = []
    RedscatterlistY = []
    for x in xList:
        y = yList[index]
        for index1 in Indexes_to_Highlight:
            if index == index1:
                RedscatterlistX.append(x)
                RedscatterlistY.append(y)
        index = index + 1
    plt.scatter(xList, yList, c = 'b')
    plt.scatter(RedscatterlistX, RedscatterlistY, c = 'r')
    plt.show()

def calculate_boundary(x, y, sensitivity): #Basically hub where every other calculation before indexing happens after it is referenced to by the main() function.
    angleList = [0]
    incrementanglecount = (360/sensitivity)
    for angle in range(1, sensitivity):
        angleList.append(angle*incrementanglecount)
    points = np.column_stack((xlist, ylist)) #Organizes them into a array of lists that has individual x y coordinates for point.
    Centroid = findCentroid(points)
    print(Centroid)
    ListWithAngle = ListWithAngleGenerator(points, Centroid)
    print(ListWithAngle)
    Indexes_to_Highlight = indexesOfPointsBoundary(ListWithAngle, Centroid, angleList, incrementanglecount)
    print(angleList)
    print("")
    plotter(x, y, Indexes_to_Highlight)
    return Indexes_to_Highlight

def InsertFluidClass(fluidthing, Indexes, xList, yList, vxList, vyList):# Initializes the fluid class with the Indexes of boundaries and their velocities.
    for index in Indexes:
      x = xList[index]
      y = yList[index]
      vx = vxList[index]
      vy = vyList[index]
      fluidthing.addDensity(x, y, 10000)
      fluidthing.addVelocity(x, y, vx, vy)


def main():
    sensitivity = 50
    #sensitivity = int(input("Enter an integer number > 3:"))
    Indexes = calculate_boundary(xlist, ylist, sensitivity)
      
    fluid = Fluid(0.1, 0.0000001, 0.00000000001)
    print(Indexes)
    InsertFluidClass(fluid, Indexes, xlist, ylist, vxlist, vylist)
    AmountOfTimestep = 50
    for i in range(0, AmountOfTimestep): #Runs the timestep 50 times.
      fluid.step()
      fluid.renderD()

if __name__ == "__main__":
    main()