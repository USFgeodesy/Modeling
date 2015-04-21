import numpy as np 
import matplotlib.pylab as plt 
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.mlab as ml
from mpl_toolkits.basemap import Basemap

'''
USF Geodesy 
Nick Voss
Functions for plotting OkadaPY outputs
'''

def plotHorizontal(xG,yG,xD,yD):
  '''
  plots the x,y component of the displacment field
  '''
  yD = np.array(yD)
  xD = np.array(xD)
  X,Y = np.meshgrid(xG,yG)
  xD= xD.reshape(np.shape(X),order='F')
  yD= yD.reshape(np.shape(Y),order='F')
  plt.quiver(X,Y,xD,yD,scale_units = 'xy',angles = 'xy',scale = .1)
  plt.title('Horizontal')
  plt.show()
  
def plotVertical(x,y,z,x0,y0):
  '''
  plots the vertical components of the vecot field
  '''
  z = np.array(z)
  X,Y = np.meshgrid(x,y)
  Z = z.reshape(np.shape(X),order = 'F')
  plt.contourf(X,Y,Z)
  plt.title('vertical')
  x = plt.colorbar()
  x.set_label('mm')
  plt.scatter(x0,y0,facecolors='none',edgecolors = (0.5,0.5,0.5))
  plt.show()
  
def plotSlip(x,y,z,latCenter,lonCenter,lons,lats,gLon,gLat,xD,yD,e,n):
  '''
  plots the slip on the fault on a map of the location
  '''
  for i in range(len(e)):
    print e[i],xD[i]
  m = Basemap(llcrnrlon=lonCenter-1.0,llcrnrlat=latCenter-1.0,urcrnrlon=lonCenter+1.0,urcrnrlat=latCenter+1.0,
            resolution='l',projection='tmerc',lon_0=lonCenter,lat_0=latCenter)
  m.drawcoastlines()
  plt.title("Slip in Costa Rica")
  
  z = np.array(z)
  x,y = m(lons,lats)
  X,Y = np.meshgrid(x,y)
  xi = np.linspace(min(x), max(x),100)
  yi = np.linspace(min(y), max(y),100)
  zi = ml.griddata(x, y, z, xi, yi)
  CS2 = m.pcolormesh(xi, yi, zi, cmap = plt.get_cmap('rainbow'))

  m.colorbar(CS2,"bottom", size="5%", pad="2%", format = '%.1f')
  m.scatter(x, y, marker = 'o', c = 'b', s = 5, zorder = 10)
  x,y = m(gLon,gLat)
  e = [east*1000 for east in e]
  n = [north*1000 for north in n]
  #Q = m.quiver(gLon,gLat,xD,yD,scale= 10000,latlon = True)
  for i in range(len(e)):
    print e[i],xD[i]
  R = m.quiver(gLon,gLat,e,n,scale = 1000,latlon = True,color = 'r')
  Q = m.quiver(gLon,gLat,xD,yD,scale= 1000,latlon = True)
  plt.show()

def plotGPS(gLon,gLat,xD,yD,latCenter,lonCenter):
  '''
  plots the GPS vectors generated by the slip on the fault
  on a map
  '''
  
  m = Basemap(llcrnrlon=lonCenter-1.0,llcrnrlat=latCenter-1.0,urcrnrlon=lonCenter+1.0,urcrnrlat=latCenter+1.0,
            resolution='l',projection='tmerc',lon_0=lonCenter,lat_0=latCenter)
  m.drawcoastlines()
# now plot.
  x,y = m(gLon,gLat)
  Q = m.quiver(gLon,gLat,xD,yD,scale= 10,latlon = True)
  
  plt.title("Slip in Costa Rica")
  plt.show()

  