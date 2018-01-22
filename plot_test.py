import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as numpy
X = numpy.asarray([1,2,3])
Y = numpy.asarray([1,2,3])
Z = numpy.asarray([1,2,3])

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(X, Y, Z, cmap= cm.jet)