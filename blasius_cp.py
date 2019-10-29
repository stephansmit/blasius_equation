import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#--------discretize
n=200000
x = np.append(np.array([0]),np.logspace(-3,1,n-1))
dx = np.diff(x)
f = np.zeros(n)
g = np.zeros(n)
h = np.zeros(n)

#--------solve
error =1
h[0] = .1
dh = .1;
while (error > 1.0e-6):
    #apply bc
    f[0] = 0.0
    g[0] = 0.0
    
    #solve equation
    for i in range(0,n-1):
        ip = i+1
        f[ip] =            g[i]*dx[i] + f[i]
        g[ip] =            h[i]*dx[i] + g[i]
        h[ip] =   -(1/2)*f[i]*h[i]*dx[i] + h[i]
    
    #update bc
    if (g[-1] < 1.0):
        h[0] = h[0] + dh;
    else:
        h[0] = h[0] - dh;
        dh = dh/2;
        h[0] = h[0] + dh;
    #calc error
    error = abs(1-g[-1])

#--------plot

fig, ax = plt.subplots()
ax.plot(f,x, label='f')
ax.plot(g,x, label='g=u')
ax.plot(h,x, label='h')
ax.set_xlabel("x", fontsize=14)
ax.set_ylabel("y", fontsize=14)
ax.tick_params(labelsize=12)
ax.legend()
ax.set_xbound(0,1.5)
ax.set_ybound(0,10)
plt.savefig("blasius_cp.png")

#--------store
df = pd.DataFrame(np.array([x,f,h,g]).T,columns=['x','f','g','h'])
df.to_csv("blasius_cp.csv",index=None, sep="\t", header=True)
