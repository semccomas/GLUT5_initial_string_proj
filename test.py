import numpy as np
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

xlist= [1,2,3]
ylist = [4,5,6,7]

X,Y = np.meshgrid(xlist, ylist)

print(xlist)
print(ylist)
print(X)
print(Y)
#ax.scatter(X, Y, color="green")


Z = (Y)
cp = ax.contour(X, Y, Z)
plt.colorbar(cp)

plt.clf()
###

import numpy as np

fig = plt.figure(figsize=(6,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height]) 


start, stop, n_values = 0, 3, 5

x_vals = np.linspace(start, stop, n_values)
y_vals = np.linspace(start, stop, n_values)
X, Y = np.meshgrid(x_vals, y_vals)

Z = X*Y

cp = plt.contourf(X, Y, Z)
plt.colorbar(cp)

ax.set_title('Contour Plot')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
plt.show()


plt.clf()



name = 'GLUT5_in'
ec = np.loadtxt('../gate_dists/extracellular/%s.EC.xvg' %name)[:,1]
ic = np.loadtxt('../gate_dists/intracellular/%s.IC.xvg' %name)[:,1]
time = np.loadtxt('../gate_dists/extracellular/%s.EC.xvg' %name)[:,0]

x,y = np.meshgrid(ec,ic)
z = np.zeros(np.shape(x))
z[0][0] = 0
z[1][1] = 1000
z[2][2] = 2000


plt.contourf(x,y,z)
plt.show()
plt.clf()


################
ec = [1, 1, 2, 2, 3, 4]
ic = [5, 5, 5, 6, 6, 6]
time =[0, 1, 2, 3, 4, 5]

x,y = np.meshgrid(ec,ic)
z = np.zeros(np.shape(x))
for n,val in enumerate(time):
    z[n][n] = val











