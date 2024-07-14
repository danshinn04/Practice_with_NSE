**Imperfect NSE Simulator - ENGR 1050 Midterm Project**

This project simulates the incompressible Navier-Stokes equations (NSE) for a 2D fluid using the immersed boundary method (IBM). The simulator approximates the flow of fluid through a grid and calculates the density and velocity fields over time.

**References**

Reference 1: http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf

Reference 2: https://mikeash.com/pyblog/fluid-simulation-for-dummies.html 

Reference 3: https://en.wikipedia.org/wiki/Immersed_boundary_method)

Reference 4: https://www.youtube.com/watch?v=JBmS--3L2eQ&ab_channel=MachineLearning%26Simulation, https://www.youtube.com/watch?v=BQLvNLgMTQE&ab_channel=MachineLearning%26Simulation 

Reference 5: https://www.youtube.com/watch?v=alhpH6ECFvQ&ab_channel=TheCodingTrain 

The Fluid class initializes a grid of size 
{𝑁
×
𝑁} to simulate the propagation of density and velocity over time using the NSE and IBM.

**Key Functions:**

**step():** Advances the simulation by one time step.

**addDensity(x, y, amount):** Adds density at a specific point.

**addVelocity(x, y, amountX, amountY):** Adds velocity at a specific point.

**renderD():** Plots the current density field.

**Boundary Generator**

The boundary generator creates a 2D dataset with randomly positioned points and velocities. It calculates the boundary of the dataset using the centroid and angle method.

**calculate_boundary(x, y, sensitivity):** Calculates the boundary points.

**InsertFluidClass(fluidthing, Indexes, xList, yList, vxList, vyList):** Initializes the Fluid class with boundary points and their velocities.

**Main Simulation**

The main function runs the simulation for a specified number of time steps and visualizes the density field.
	
**Results** (You get velocity fields as a series of time)

![image](https://github.com/user-attachments/assets/89489ba5-992f-4145-9b19-b8ddbaa9dade)
![image](https://github.com/user-attachments/assets/36e8b938-836c-40a9-8099-40cf3dc8c7c6)
![image](https://github.com/user-attachments/assets/8179aedd-f776-450c-b4ee-508ab0c151fa)


