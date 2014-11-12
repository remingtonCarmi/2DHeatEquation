Solve the 1D heat equation
u_t-1/Pe* (u_xx+ u_yy) = 0 for x in 0  2pi y in 0 2pi
periodic BC 
u(x=0,y,t) = u(x=2pi,y,t)
u(x,y=0,t) = u(x,y=2pi,t)
u(x,y,0) = cos(x)+sin(y)

Depenencies: 
 matplotlib.pyplot as plt (for plotting)

***********************************************************************

Solver use a center schemes for the space differenciation and explicit euler for the time
to be stable require :  dt/min(dx*dx) <0.5Pe

save the plot in result
to convert in gif use the command : 
 convert result/*png 2DHeatEquation.gif
