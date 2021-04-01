# Title: Computational Fluid Dynamics Simulation of the Vorticity

# Simulator developer: Md. Nazmus Sanib Chowdhury

"""Reference:  A Finite Difference Code for the Navier-Stokes
Equations in Vorticity/Stream Function
Formulation

Instructor: Hong G. Im
University of Michigan
Fall 2001

University of Michigan."""



"""Suggested initial values:
hlen = 1
vlen = 1
xGrid = 32
yGrid = 32
dt = 0.005
timeStep = 100
Re = 10"""



""" Instruction to run the simulation:
Step 1: Pre-processing:

Define a object in whcih call the simulator class and pass the argument as the initial condition of the simulation, e.g.,
>>> simulator = preProcessor(1,1,16,16,0.005,100,10)

Step 2: Solver:

Run the simulation by 
>>> simulator.solver()

Step 3: Post processing

See the graphical output of the simulation,
for stream-plot:
>>> simulator.postProcess(1)

for vorticity plot:
>>> simulator.postProcess(2)
"""




class simulator:
    def __init__(self,hlen,vlen,xGrid,yGrid,dt,timeStep,Re,precision=200,tm=0.0):
        
        import numpy as np
        
        self.hlen = hlen
        self.vlen = vlen
        self.xGrid = xGrid

        # Setting the condition that the grid number along the x-axis can't be upto 16 or less than 9
        if self.xGrid>16:
            self.xGrid = 16
        elif self.xGrid<9:
            self.xGrid = 9

            
        self.yGrid = yGrid

        # Setting the condition that the grid number along the x-axis can't be upto 16 or less than 9
        if self.yGrid>16:
            self.yGrid = 16
        elif self.yGrid<9:
            self.yGrid = 9

        # Setting the condition that the grid number along the x-axis and y-axis can't be unequal    
        if self.xGrid!=self.yGrid:
            self.yGrid = self.xGrid

            
        self.stream = np.zeros((self.xGrid,self.yGrid),dtype = np.longdouble)
        #self.stream = np.ones((self.xGrid,self.yGrid),dtype = np.longdouble)

        self.dt = dt
        self.dx = hlen/(self.xGrid-1)
        
        self.timeStep = timeStep

        
        # Setting the condition that the timeStep can't be exceeded 30
        if self.timeStep!=100:
            self.timeStep = 100
            
        self.dy = vlen/(self.yGrid-1)
        
        self.vor = self.stream.copy()
        self.Re = Re
        
        self.precision = precision

        # Setting the condition that the precision have to set to 200
        if self.precision!=200:
            self.precision = 200
            
        
        self.tm = tm

        # Setting the condition that the initial time can't be unlike than 0.0
        if self.tm!=0.0:
            self.tm = 0.0

        self.vor[:,:] = .113
        #self.vor[:,:] = 0.314
        #self.vor[:,:] = 1.314

        
        
              
        


    def solver(self):
        import numpy as np
        
        def p(aq):
            import numpy as np
            #return float(np.format_float_positional(aq,unique =  False, precision = self.precision))
            return np.round(aq,decimals = self.precision)



        
        for n in range(self.timeStep):
            # Variables
            #self.stream[:,:] = 1.315
            
            
            maxErr = 0.001
            sf = self.stream
            vt = self.vor
            h = self.dx
            Re = self.Re
            sigma = 1
            dt = self.dt*sigma
            nx = self.xGrid
            ny = self.yGrid
            beta = 1.5
            w = sf.copy()
            
            maxIt = 100
           
           


            # Upadating the Stream function in space along with Navier-Stokes (N-S) equation
            for time in range(maxIt):
                for i in range(1,nx-1):
                    for j in range(1,ny-1):
                        if abs(self.stream[j,i])>=1e4:
                            break
                        self.stream[j,i]= 0.25*beta*(p(sf[j,i+1])+p(sf[j,i-1])+p(sf[j+1,i])+p(sf[j-1,i])+p(vt[j,i])*(h**2))+(1-beta)*p(sf[j,i])
                        #print(self.stream)
                        
                        
                        #print("Sub-loop: {time}".format(time = time))
                err = 0.0    # Stop iteration if converged
                for i in range(nx-1):
                    for j in range(ny-1):
                        err = err+abs(w[j,i]-sf[j,i])
                        if err.any()<=maxErr:
                            break





            # Boundary conditions

            
            vt[0,1:nx-2]=-2.0*sf[1,1:nx-2]/(h**2)
            vt[ny-1,1:nx-2]=-2.0*sf[ny-2,1:nx-2]/(h**2)-2.0/h
            vt[1:ny-2,0]=-2.0*sf[1:ny-2,1]/(h**2)
            vt[1:ny-2,nx-1]=-2.0*sf[1:ny-2,nx-2]/(h**2)
            



            # Updating the Vorticity in space along with N-S eq.
            for i in range(1,nx-1):
                for j in range(1,ny-1):
                    if abs(self.vor[j,i])>=1e4:
                        break
                    w[j,i] = 2000*(0.25*(1/h**2)*(((p(sf[j+1,i])-p(sf[j-1,i])))*((p(vt[j,i+1])-p(vt[j,i-1])))-((p(sf[j,i+1])-p(sf[j,i-1])))*((p(vt[j+1,i])-p(vt[j-1,i]))))+(1/Re)*((p(vt[j,i+1])+p(vt[j,i-1])+p(vt[j+1,i])+p(vt[j-1,i])-4*p(vt[j,i]))/h**2))
                    

           
            
            self.vor[2:ny-2,2:nx-2]=2000*vt[2:ny-2,2:nx-2]+self.dt*w[2:ny-2,2:nx-2]
            self.vor = self.vor/2000
            
            #print(self.vor)
            self.tm += self.dt

            print("\n Loop: {Loop}".format(Loop = n))
            
            print("\n Time: {Time}s".format(Time = (np.round(self.tm,decimals = 3))))
            
            

                   

        return("Stream:",sf,"""//##########################//""","Vorticity:",vt)           



    def postProcess(self,num):
        #print("Stream:",self.stream)
        #print("Vorticity:",self.vor)
        import matplotlib.pyplot as plt
        import seaborn as sn

        if num==1:
            sn.heatmap(self.stream,annot = False)
            plt.title("Heatmap Plot of the Stream function")
            plt.show()
            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            cax = ax.matshow(self.stream,cmap = 'hsv',interpolation = 'nearest')
            fig.colorbar(cax)
            plt.title("Stream function Plot")
            plt.show()
        else:
            sn.heatmap(self.vor,annot = False)
            plt.title("Heatmap Plot of the Vortcity")
            plt.show()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            cax = ax.matshow(self.stream,cmap = 'hsv', interpolation = 'nearest')
            fig.colorbar(cax)
            plt.title("Vorticity Plot")
            plt.show()
