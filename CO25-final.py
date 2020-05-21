#This is computing project CO25 - Laplace's equation
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# key formula for the numerical method for solving
# where R(i,j)=psi(i,j+1)+psi(i,j-1)+psi(i-1,j)+psi(i+1,j)-4psi(i,j) for each iteration


n =7   #size of the grid 
init_psi = np.zeros( (n, n))
init_psi[:,0] = init_psi[0,:]= 0  
init_psi[n-1,:]=[np.sinh(1)*np.sin(i*1/(n-1)) for i in range (n)] #sin x sinh 1
init_psi[:,n-1]=[np.sinh(i*1/(n-1))*np.sin(1) for i in range (n)] #sin 1 sinh y

m = 1000  #number of iteration 
np.set_printoptions(threshold=200)

#defining list psi and hist-historical values of elements in psi
psi = []
hist=[]

#over-relaxation constant and the optimum value for it
a = 1.35
a_opt = 2/(1+np.sin(np.pi/(n)))


#criteria for convergence
conv = 1e-12

'''
The flowing functions are for graph plotting that help illustrating datas
'''

'''
plotgraph(matrixZ,size) gives the contour plot of matrixZ
'''
def plotgraph(matrixZ,size):
    plt.figure(1)
    x = np.linspace(0, 1, size)
    y = np.linspace(0, 1, size)
    X, Y = np.meshgrid(x, y)

    # Set colour interpolation and colour map
    colorinterpolation = 50
    colourMap = plt.cm.jet

    #producing the contour plot
    CS = plt.contourf(X, Y, matrixZ, colorinterpolation, cmap=colourMap)

    # Set Colorbar
    plt.colorbar()
    #CS = plt.contour(X, Y, matrixZ)
    plt.title('contour plot of psi')

    #turn off interactive mode
    #plt.ioff()
    plt.show()



    '''
    The function plotimag(matrixZ,size,fig) illustrates the dynamical change of matrixZ in iteration on figure fig
    '''
def plotimag(matrixZ,size,fig):
    fig.clf()
    x = np.linspace(0, 1, size)
    y = np.linspace(0, 1, size)
    X, Y = np.meshgrid(x, y)
    
    #shows the value of matrixZ by colour image with gaussian bluring and colour range (5,-5)
    plt.title("Image for psi changing in iterations")
    plt.imshow(matrixZ, interpolation='gaussian',vmax = 1,vmin= -1,origin="lower")
    plt.colorbar()
    
    #this determins the pause between each new frame for showing the updated value of psi
    plt.pause(0.01) 
    
    
'''
this function is for ploting the 3D image of psi
'''
def plot3D(matrixZ,size):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Make data.
    X = np.linspace(0, 1, size)
    Y = np.linspace(0, 1, size)
    X, Y = np.meshgrid(X, Y)
    # Plot the surface.
    surf = ax.plot_surface(X, Y, matrixZ, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()



#this function is to plot data relate to x column of the Matrix
def hisgraphplot(Matrix,x):                     
    plt.figure(3)

    #converting Matrix to array
    arr = np.array(Matrix)

    #ploting the x coulumn of the Matrix
    plt.plot(arr[:,x])
    if x==0:
        str='upper'
    if x==1:
        str='middle'
    if x==2:
        str='lower'
    plt.ylabel('hist_values_'+str)
    plt.xlabel('the n th iteration')
    plt.show()






'''
function that solve psi and output dynamical plot for change of psi 
input: number of iteration, initial psi  and over-relaxation parameter a
output: array containing info of psi and historical values of elements in psi
'''
def solvlaplace(num_iter,init_psi,a):
    psi = np.array(init_psi)                      #input initial value of psi
    fig = plt.figure()
    size = np.size(psi,1)

    #initializing the residue
    R = np.zeros((size,size))        

    #turn on interactive mode
    plt.ion()
    for k in range (num_iter):
        idcator = 0
        for i in range (1,n-1):
            for j in range(1,n-1):
                R[i,j] = psi[i,j+1]+psi[i,j-1]+psi[i-1,j]+psi[i+1,j]-4*psi[i,j]


                #identifying if the convergence is reached 
                if abs((a/4)*R[i,j]) > conv:
                    idcator=1

                #updating the value of psi_ij
                psi[i,j] =psi[i,j]+(a/4)*R[i,j]


        #ploting the coloured image for psi at current step
        plotimag(psi,size,fig)

        #storing the historical values of the psi in 3 different locations (upper half, middle, lower half)
        hist.append([psi[int(n*3/4),int(n*3/4)],psi[int(n/2),int(n/2)],psi[int(n/4),int(n/4)]])

        #terminating the iteration
        if idcator == 0:
            break


    hist_out = np.array(hist)                       #converting data to array
    np.savetxt("hist_value.txt",hist_out)           #saving hist_out to file hist_value.txt
    psi_hist = np.array([psi,hist])                 #creating array made of matrix psi and the historical value of elements in psi
    print ("end")
    
    
    plt.ioff()
    plt.show()
    return psi_hist

'''
main function solving psi
input: number of iteration, initial psi and over-relaxation parameter a
output: array containing info of psi and historical values of elements in psi 
plot: producing contour plot for psi
'''
def solvlaplace2(num_iter,init_psi,a):         #function (number of iteration, initial psi func, size of grid, a is over-relaxation constant)
                                                    #->output( numerical solution psi and historical value of psi at 3 points;plot for contour line and histroical values) 
    psi = np.array(init_psi)                                  #input initial value of psi
    hist.append([psi[int(n*3/4),int(n*3/4)],psi[int(n/2),int(n/2)],psi[int(n*1/4),int(n*1/4)]])     #appending the initial value to the record
    size = np.size(psi,1)
    R = np.zeros((size,size))        
   
    for k in range (num_iter):
        idcator = 0
        for i in reversed(range (1,n-1)):
            for j in reversed(range(1,n-1)):

                #updating the residual
                R[i,j] = psi[i,j+1]+psi[i,j-1]+psi[i-1,j]+psi[i+1,j]-4*psi[i,j]

                #identifying if the convergence is reached 
                if abs((a/4)*R[i,j]) > conv:
                    idcator=1

                #updating the value of psi_ij
                psi[i,j] =psi[i,j]+(a/4)*R[i,j]


        #storing the historical values of the psi in 3 different locations (upper half, middle, lower half)
        hist.append([psi[int(3*n/4),int(3*n/4)],psi[int(n/2),int(n/2)],psi[int(n/4),int(n/4)]])

        #terminating the iteration
        if idcator == 0:
            break
          
    hist_out = np.array(hist)
    psi_hist = np.array([psi,hist])         #creating array made of matrix psi and the historical value of elements in psi

    #plot the contour graph of psi
    plotgraph(psi,size)
    return psi_hist




'''
this part try to improve the efficiency for solving the laplace equation by using the mean value property of solutions of laplace equation
we can choose init-psi closer to that of real solution
'''
def init_psi_mean_value(Matrix):          #by applying the mean value property of harmonic function to choose closer initial value
    sum = 0                               #vairable that represent the total sum of boundary value
    for i in init_psi[:,0]:
        sum += i
    for i in init_psi[:,n-1]:
        sum += i
    for i in init_psi[0,:]:
        sum += i
    for i in init_psi[n-1,:]:
        sum += i
    aver = sum/(4*n)                    #calculating the approximate average of the boundary value and thus the mean value for the interior
    
    for j in range(1,n-1):              #take values of init_psi to be the average value of the boundary in interior of the region
        for k in range(1,n-1):
            Matrix[j,k]=aver
    return Matrix
  
    

'''
for the example with specific boundary on script has analytical solution
psi=sinxsinhy, 
therefore this function is for comparing numerical result and true value
input:size of grid
output:contour graph of analytical solution
'''

def analyticalsol(size):
    psi_true=np.zeros( (size, size))
    for i in range(size):
        for j in range(size):
            psi_true[i][j]=np.sin(j/(size-1))*np.sinh(i/(size-1))
    plotgraph(psi_true,size)
    print(psi_true)
    np.savetxt("psi_true"+"_size_"+str(size)+".txt",psi_true)
    return psi_true






solvlaplace2(m,init_psi,a)
#np.savetxt("hist_non_zero_start_a="+str(a)+"_size_"+str(n)+"iter"+str(m)+".txt",solvlaplace2(m,init_psi,n,a)[1])
#np.savetxt("psi_non_zero_start_a="+str(a)+"_size_"+str(n)+"iter"+str(m)+".txt",psi)

