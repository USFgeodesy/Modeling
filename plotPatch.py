import numpy as np 
from matplotlib import colorbar
import matplotlib.pylab as plt 
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.mlab as ml
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import matplotlib.cm as cmx
import matplotlib.colors as colors
from matplotlib.collections import PatchCollection
def readFault(filename,nh,nv):
  #this is for if you have just a x,y,z file it is not working right now 
    x =np.loadtxt(filename)
    x = x.transpose()
    Lo1 = x[0]
    La1 = x[1]
    z = x[2]
    Lo1 = np.array(Lo1)
    La1 = np.array(La1)
    z = np.array(z)
    # this is due to the file format of Yan it is the Length and width of the grid
    #nh=19
    #nv=17
    Lats = []
    Lons = []
    zs = []
    for i in range(nh-1):
        for j in range(nv-1):
            ii=(j)*nh+i
            Lats.append([La1[ii],La1[ii+nh],La1[ii+1+nh],La1[ii+1]])
            Lons.append([Lo1[ii],Lo1[ii+nh],Lo1[ii+1+nh],Lo1[ii+1]])
            zs.append([z[ii],z[ii+nh],z[ii+1],z[ii+1+nh]])
            
    return Lats,Lons,zs
  
  
def draw_screen_poly( lats, lons, m,slip):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, facecolor= slip, alpha=0.4 )
    return poly      
def plotSlipPatch(latCenter,lonCenter,lons,lats,gLon,gLat,xD,yD,e,n,slip):
  '''
  plots the slip on the fault on a map of the location on polygons representing patch 
  '''
  fig = plt.figure(figsize = (10,10))
  m=Basemap(llcrnrlon=lonCenter-1.0,llcrnrlat=latCenter-1.0,urcrnrlon=lonCenter+1.0,urcrnrlat=latCenter+1.0,resolution='h',projection='tmerc',lon_0=lonCenter,lat_0=latCenter)
  m.drawcoastlines()
  m.drawparallels(np.arange(9,11,0.5),labels=[1,0,0,1], linewidth=1.0)
  m.drawmeridians(np.arange(-86.6,-64.5,0.5),labels=[1,0,0,1], linewidth=1.0)
  plt.title("SSE 14")
  cm = plt.get_cmap('rainbow')
  cNorm  = colors.Normalize(vmin = min(slip), vmax=max(slip))
  #cNorm  = colors.Normalize(vmin = 0, vmax=150)
  scalarMap = cmx.ScalarMappable(norm=cNorm, cmap= cm)
  array = scalarMap.set_array
  patches = []
  for i in range(len(lons)):
    colorVal = scalarMap.to_rgba(slip[i])
    patch = draw_screen_poly(lats[i],lons[i],m,colorVal)
    patches.append(patch)
  p = PatchCollection(patches,cmap = cm, alpha = 0.9)
  p.set_array(np.array(slip))
  ax = plt.gca()
  ax.add_collection(p)
  cb1 = plt.colorbar(p)
  #cb1.set_clim(0,100)
  x,y = m(gLon,gLat)
  e = [east for east in e]
  n = [north for north in n]
  #for i in range(len(e)):
   # print e[i],xD[i]
  R = m.quiver(gLon,gLat,e,n,scale = 50,latlon = True,color = 'r')
  Q = m.quiver(gLon,gLat,xD,yD,scale= 50,latlon = True,color = 'black')
  qk1 = plt.quiverkey(Q, 0.55, 0.95, 10, 'Modeled (10 mm)', labelpos='E',fontproperties={'weight': 'bold'})
  qk2 = plt.quiverkey(R, 0.55, 0.9, 10, 'Observed (10 mm)', labelpos='E',fontproperties={'weight': 'bold'})
  #l = [float(x) for x in lons]
  #la = [float(y) for y in lats]
  #x, y = m( l, la )
  #cset2 = m.contour(x,y,slip,latlon = True)
  cb1.set_label('mm')
  plt.savefig('2014sseLats.png')
  plt.show()
  
def plotPatch(latCenter,lonCenter,gLon,gLat,xD,yD,e,n,slip):
  Lats,Lons,zs = readFault('/home/nvoss/Geodesy/OkadaPY/geometry.txt',19,17)
  plotSlipPatch(latCenter,lonCenter,Lons,Lats,gLon,gLat,xD,yD,e,n,slip)
  #plt.colorbar()