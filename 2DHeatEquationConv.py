#!/usr/local/bin/python
from math import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import os
import errno

## 2DAdvDiffEq.py
## Remi Carmigniani
## Solves the 2D Diffusion Equation :
## pd_t u[x,y,t]  - 1/Pe(pd_xx u[x,y,t]+pd_yy u[x,y,t])  = 0 (x,y) in [0 2pi]^2
## BC periodic
## IC u[x,y,0]=cos(x)+cos(y) 
## Convergence test\

#IMPORTANT TO COPY A MATRIX USE u_old = [row[:] for row in u_arr] this will unlink them...

##Physical parameter
Pe = 1.0
L=2*pi

##################################################################################################################
############################			Useful functions		      ############################
##################################################################################################################
#L2 error using Simpson rule 
def errorL2(f,dx,dy,nx,ny):
	error = f[0][0]+f[0][ny-1]+f[nx-1][0]+f[nx-1][ny-1]
	for j in range(1,int((ny-1)/2)+1):
		error = error+ 4.*(f[0][2*j-1]+f[nx-1][2*j-1])
	for j in range(1,int((ny-1)/2)):
		error = error+ 2.*(f[0][2*j]+f[nx-1][2*j])
	for i in range(1,int((nx-1)/2)+1):
		error = error+ 4.*(f[2*i-1][0]+f[2*i-1][ny-1])
	for i in range(1,int((nx-1)/2)):
		error = error+ 2.*(f[2*i][0]+f[2*i][ny-1])
	for i in range(1,int((nx-1)/2)+1):
		for j in range(1,int((ny-1)/2)+1):
			error = error + 16.*f[2*i-1][2*j-1]
        for i in range(1,int((nx-1)/2)+1):
		for j in range(1,int((ny-1)/2)):
			error = error + 8.*f[2*i-1][2*j]
        for i in range(1,int((nx-1)/2)):
		for j in range(1,int((ny-1)/2)+1):
			error = error + 8.*f[2*i][2*j-1]
        for i in range(1,int((nx-1)/2)):
		for j in range(1,int((ny-1)/2)):
			error = error + 4.*f[2*i][2*j]		
	return (error*1./9.*dx*dy)

## Initial conditions
def iCond(x,y):
    return .5*cos(x) + .5*cos(y)


#Exact solution
def u_exact(x,y,t,Pe):
	return exp(-t/Pe)*(.5*cos(x) + .5*cos(y))

## Create a one row vector from the matrix u_arr
def matrixToVec(u,nx,ny):
	u_vec = [0 for i in xrange(ny*nx)]
	for i in range(0,nx):
		for j in range(0,ny):
			u_vec[i*ny+j] = u[i][j]
		
    	return u_vec

## Reverse
def vecToMatrix(u,nx,ny):
	Umat=[[0 for i in xrange(ny)] for i in xrange(nx)]
	for i in range(0,nx):
		for j in range(0,ny):
			Umat[i][j]=u[i*ny+j]
    	return Umat
#Periodic BC
def period(k,n):
	if k==-1:
		return n-2
	else:
		if k==n:
			return 1
		else:
			return k 
##################################################################################################################
############################				End			      ############################
##################################################################################################################
#Simulation parameters
tend = 0.5


error=[]
## Series of N to try
N = [25,51,101]
for k in range(0,3):
        t=0
	## Discretization parameter 
	Nx=N[k]
	Ny=N[k]
	dx=L/float(Nx-1)
	dy=L/float(Ny-1)

	## time parameters 
	dt =0.1*dx*dx
	#to make sure tend is reached and stability ok
	#dt = tend/float(int(tend/dt)+1)

	## scheme parameters
	# Laplacian terms
	lx = 1.*dt/dx/dx
	ly=1.*dt/dy/dy


	## Figure numbering 
	numb = 0
	uval = iCond(0,0)
	u_arr = [[0 for i in xrange(Ny)] for i in xrange(Nx)]
	for i in range(0,Nx):
		for j in range(0,Ny):
    			uval = iCond(i*dx,j*dy)
    			u_arr[i][j] = uval
       
	## Construction of the update matrix
	#A = [[0 for i in xrange(Ny*Nx)] for i in xrange(Nx*Ny)]
	
	#for i in range(0,Nx):
	#	for j in range(0,Ny):      
	#		A[i*Ny+j][i*Ny+j]=1-2.*(lx+ly)
	#		A[i*Ny+j][i*Ny+(j+1)%Ny] = ly
	#		A[i*Ny+j][i*Ny+(j-1)%Ny] = ly  
	#		A[i*Ny+j][((i+1)%Nx)*Ny+j] = lx
	#		A[i*Ny+j][((i-1)%Nx)*Ny+j] = lx
	#convert u_arr to u_vec
	#u_vec = matrixToVec(u_arr,Nx,Ny)
	##Time loop
	while t<tend-1e-6 :
		t=t+dt
		#print repr(t)
		#update u_vec
		#u_vec = np.dot(A,u_vec)
		u_old = [row[:] for row in u_arr]
		for i in range(0,Nx):
			for j in range(0,Ny):
				u_arr[i][j]=u_old[i][j]*(1.-2.*(lx + ly)) + (u_old[period(i-1,Nx)][j]+u_old[period(i+1,Nx)][j])*lx+(u_old[i][period(j-1,Ny)]+u_old[i][period(j+1,Ny)])*ly

        print repr(u_arr[Nx-1][Ny-1])			    
	#u_arr=vecToMatrix(u_vec,Nx,Ny)
	##Calculate the error
	
	f = [[0 for i in xrange(Ny)] for i in xrange(Nx)]
	for i in range(0,Nx):
		for j in range(0,Ny):
			f[i][j] =u_arr[i][j]- u_exact(i*dx,j*dy,t,Pe)
			f[i][j]=f[i][j]**2
	
        err=sqrt(errorL2(f,dx,dy,Nx,Ny))
	error.append(err)
	print 'Calculate the error  for ' + repr(Nx) + ' : ' + repr(err)

plt.loglog(N, error)
plt.title('Error')
plt.savefig('Error.png')

print 'Simulation Completed without error'

		
        	

       
	
	


