from numpy import pi, cos, sin, arccos, arange
import mpl_toolkits.mplot3d
import matplotlib.pyplot as pp

num_pts = 3  
indices = arange(0, num_pts, dtype=float) + 0.5

phi = arccos(1 - 2*indices/num_pts)
theta = pi * (1 + 5**0.5) * indices

x, y, z = cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi);

fig = pp.figure()
ax = fig.add_subplot(111, projection='3d')

print x,y,z
ax.scatter(x, y, z);

num_pts = 32-3
indices = arange(0, num_pts, dtype=float) + 0.5

phi = arccos(1 - 2*indices/num_pts)
theta = pi * (1 + 5**0.5) * indices

x1, y1, z1 = cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi);

ax.scatter(x1, y1, z1);

pp.show()
