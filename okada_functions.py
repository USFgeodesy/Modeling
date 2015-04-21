import numpy as np 
import dC3D as dC
import readData as rd 
from okada_wrapper1 import dc3d0wrapper, dc3dwrapper
'''
Functions needed for okadaPy
USF Geodesy Lab Nick Voss
October 20, 2014
'''

#compute the North,East,and Up around the reference point
def neu(lat,lon):
  earthRadius = 6371.0 #in km 
  latR = np.radians(lat) #convert lat and lon to radians 
  lonR = np.radians(lon)
  z0 = earthRadius*np.sin(latR)
  x0 = earthRadius*np.cos(latR)*np.cos(lonR)
  y0 = earthRadius*np.cos(latR)*np.sin(lonR)
  ux = x0/earthRadius
  uy = y0/earthRadius
  uz = 0.0 
  r1 = np.sqrt(x0**2+y0**2)
  ex = -y0/r1
  ey = x0/r1
  ez = 0.0
  nx=-uz*ey
  ny=uz*ex
  nz=ux*ey-uy*ex
  r2=np.sqrt(nx*nx+ny*ny+nz*nz)
  nx=nx/r2
  ny=ny/r2
  nz=nz/r2
  return x0,y0,z0,nx,ny,nz,ex,ey,ez,ux,uy,uz

# Read in text file with GPS locations,offsets and errors
def readGPS(gps_file):
  data = np.loadtxt(file_name,skiprows = 1)
  data = data.transpose()
  gLat = data[0]
  gLon = data[1]
  dispX = data[2]
  dispY = data[3]
  dispZ = data[4]
  errX = data[5]
  errY = data[6]
  errZ = data[7]
  gpsData = [gLat,gLon,dispX,dispY,dispZ,errX,errY,errZ]
  return gpsData

#create grid to compute displacements at
def createGrid(minX,minY,maxX,maxY,nX,nY):
  x = np.linspace(minX,maxX,nX)
  y = np.linspace(minY,maxY,nY)
  return x,y

#for rotating by strike relative to EAST
def rotate(angle,x,y):
  theta = np.radians(angle)
  xRot = x*np.cos(theta)-y*np.sin(theta)
  yRot = x*np.sin(theta)+y*np.cos(theta)
  return xRot,yRot

#compute the compenent of east and north with respect to reference
def lalo2xy(la,lo,LaCen,LonCen):
  x0,y0,z0,nx,ny,nz,ex,ey,ez,ux,uy,uz = neu(LaCen,LonCen)
  earthRadius = 6371.0 #in m 
  laR = np.radians(la)
  loR = np.radians(lo)
  z1=earthRadius*np.sin(laR)
  x1=earthRadius*np.cos(laR)*np.cos(loR)
  y1=earthRadius*np.cos(laR)*np.sin(loR)
  vx = x1-x0
  vy = y1-y0
  vz = z1-z0 
  e = ex*vx+ey*vy
  n = nx*vx+ny*vy+nz*vz
  u = ux*vx+uy*vy+uz*vz
  return e,n,u 

#go back to lat lon coordinates
def xy2lalo(x,y,z):
  r1=sqrt(x*x+y*y)
  ex=-y/r1
  ey=x/r1
  ez=0.
  vx=(x*ex+y*nx+z*ux)
  vy=(x*ey+y*ny+z*uy)
  vz=(x*ez+y*nz+z*nz)
  x1=vx+x
  y1=vy+y
  z1=vz+z
  r1=sqrt(x1*x1+y1*y1)
  la=np.degrees(np.arctan2(z1,r1))
  lo=np.degrees(np.arctan2(y1,x1))
  return la,lo

def faultcen2ref(lat,lon):
  #lat and lon are the reference lat and lon for the grid 
  #reutrn lat and lon of reference in local reference frame 
    y,x = neu(lat,lon)
    return x,y 

def faultPatch2ref(lat,lon):
  East,North,Vert = [],[],[]
  #lat and lon are arrays of fault patch centers
  #returns patches in local reference frame
  for i in range(len(lat)):
      e,n,u = lalo2xy(lat[i],lon[i])
      East.append(e)
      North.append(n)
      Vert.append(u)
  return East,North,Vert 

def strike2east(strike):
  #converts Strike to relative to east instead of North for 
  #easy rotation 
  angle = [90.0-s for s in strike]
  return angle

#python wrapper of Okada DC3D fortran code (Okada 1992)
def okada(alpha,x,y,z,depth,dip,al1,al2,aw1,aw2,ssDis,dsDis,tDis):
    #ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret = dC.dc3d(alpha,x,y,z,depth,dip,al1,al2,aw1,aw2,ssDis,dsDis,tDis)
    success,u,grad= dc3dwrapper(alpha,[x,y,z],depth,dip,[al1,al2],[aw1,aw2],[ssDis,dsDis,tDis])
    iret = success
    ux = u[0] 
    uy = u[1]
    uz = u[2]
    uxx = grad[0]
    uyx = grad[0]
    uzx = grad[0]
    uxy = grad[1]
    uyy = grad[1]
    uzy = grad[1]
    uxz = grad[2]
    uyz = grad[2]
    uzz = grad[2] 
    return ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret

#computes easting and northing from reference to obervation or fault 
def eastNorth(x0,y0,cX,cY):
  ''' computes easting and northing to center of fault'''

  eF = cX-x0 
  nF = cY-y0
  return eF,nF

def computeDisp(faultData,x,y,z,ss,ds,t):
 # print faultData
 #[alpha,depth,dip,strike,xCen,yCen,L,W]
  dx = float(0.0) 
  dy = float(0.0)
  dz = float(0.0)
  alpha = faultData[0]
  depth = faultData[1]
  dip = faultData[2]
  strike = faultData[3]
  x0 = faultData[4]
  y0 = faultData[5]
  L = faultData[6]
  W = faultData[7]
  angle = strike2east(strike)
  for i in range(len(depth)):
      #compute easting and northing to the fault patch 
      ef,nf = eastNorth(x0[i],y0[i],0.0,0.0)
      #shift grid point 
      x1 = x+ef 
      y1 = y+nf
      #rotate grid point
      x2,y2 = rotate(-angle[i],x1,y1)
      #compute output
     # (alpha,x,y,z,depth,dip,al1,al2,aw1,aw2,ssDis,dsDis,tDis)
      output = okada(alpha,x2,y2,z,depth[i],dip[i],-L[i]*0.5,L[i]*0.5,-W[i]*0.5,W[i]*0.5,ss[i],ds[i],t[i])

      #rotate displacement back to north east west south coordinated
      xOut,yOut = rotate(angle[i],output[0],output[1])
      #sum outputs for every patch 
      dx = dx + xOut
      dy = dy + yOut
      dz = dz + output[2]
  return dx,dy,dz 
   
def gridDisp(xRan,yRan,faults,ss,ds,t):
  #Compute Displacement from given fault on a uniform Grid
  x = np.linspace(-xRan/2.0-200.0,xRan/2.0+200.0,25.0)
  y = np.linspace(-yRan/2.0-200.0,yRan/2.0+200.0,25.0)
  X,Y,Z = [],[],[]
  for i in range(len(x)):
    for j in range(len(y)):
      dx,dy,dz = computeDisp(faults,x[i],y[j],0.0,ss,ds,t)
      X.append(dx)
      Y.append(dy)
      Z.append(dz)
  return X,Y,Z,x,y  

def gpsDisp(x,y,z,faults,ss,ds,t):
  #Compute Displacement at GPS coordinates 
  X,Y,Z = [],[],[]
  for i in range(len(x)):
    dx,dy,dz = computeDisp(faults,x[i],y[i],0.0,ss,ds,t)  
    X.append(dx)
    Y.append(dy)
    Z.append(dz)
  return X,Y,Z
