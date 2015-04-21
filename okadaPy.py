import numpy as np 
import sys 
import okada_functions as fun 
import readData as rD
import outputData as oD
import plotData as pD 
import plotPatch as pP
import matplotlib.pylab as plt

'''
okadaPy 
USF Geodesy Lab
Nick Voss October 20, 2014

okadaPy is the main program of of a collection of scripts:
okada_functions contains the neccessary routines
  -uses dC3D fortran code from Okada 1992
readData reads in fault and station 
outputData outputs the computed displacements to text file
plotData plots the outputData

to use: python okadaPy.py faultFile gps1File slipFile 
'''

#read in Fault
faultFile = sys.argv[1]

 #x y  depth      dip  strike    L       W 
xCen,yCen,depth,dip,strike,L,W,latCenter,lonCenter,lats,lons = rD.readDataFile(faultFile)


#read in GPS 
gpsFile = sys.argv[2]
gpsX,gpsY,gLon,gLat,gZ = rD.readGPS(gpsFile,latCenter,lonCenter)
#read in Slip 
slipFile = sys.argv[3]
index,ss,ds,t = rD.readSlip(slipFile)
#paramater = np.genfromtxt('/home/nvoss/Geodesy/PEST_OKADAPY/co_12_svd_wt.par',skiprows = 1,usecols = (1,2),dtype = None)
#value = paramater.transpose()[0]
#ds = []
#ss = []
#for i in range(len(value)):
 # if i%3 ==  0:
  #  ss.append(value[i])
  #elif value[i]!=0:
   # ds.append(value[i]-1.0)
#define elastic constants
alpha = 2.0/3.0  #lambda=mu

#set up fault data array
faults = [alpha,-depth,dip,strike,xCen,yCen,L,W]

#compute displacement 

#compute displacement on Grid 
xRange = max(xCen) - min(xCen)
yRange = max(yCen) - min(yCen)
#dxGrid,dyGrid,dzGrid,xGrid,yGrid = fun.gridDisp(xRange,yRange,faults,ss,ds,t)

#compute displacement at gps locations
dxGPS,dyGPS,dzGPS = fun.gpsDisp(gpsX,gpsY,gZ,faults,ss,ds,t)
stationData = np.genfromtxt('SSE14.txt',skiprows = 0, usecols = (0,4,5,6),dtype = None,delimiter = " " )
station,e,n,v = [],[],[],[]
for i in range(len(stationData)):
    station.append(stationData[i][0])
    e.append(stationData[i][1])
    n.append(stationData[i][2])
    v.append(stationData[i][3])
#print e
#print n 
#print v 
#print e , dxGPS
#output displacement
#at Grid
#oD.outputGrid(xGrid,yGrid,dxGrid,dyGrid,dzGrid)
#at GPS 
oD.outputGPS(gLon,gLat,dxGPS,dyGPS,dzGPS)
#plotDisplacement 
#pD.plotHorizontal(xGrid,yGrid,dxGrid,dyGrid,211)
#plt.figure()
#pD.plotVertical(xGrid,yGrid,dzGrid,212,xCen,yCen)
#pD.plotSlip(xCen,yCen,ds,latCenter,lonCenter,lons,lats,gLon,gLat,dxGPS,dyGPS,e,n)#
slip = []
for i in range(len(ds)):
	slip.append(np.sqrt(ss[i]**2+ds[i]**2)) 
def calcMo(l,w,slip,shearMod):
  mo = (np.array(l)*1000.0)*(np.array(w)*1000.0)*(np.array(slip)/1000.0)*shearMod
  return np.sum(mo)
def moMag(mo):
  mw = (2.0/3.0)*np.log10(mo)-6.0
  return mw 
shearMod = 32000000000 #Pa
moment = calcMo(L,W,slip,shearMod)
mwMag = moMag(moment)
plate = []
for i in range(len(slip)):
	plate.append(80*62)

print 'Magnitude : ',mwMag
accMo = calcMo(L,W,plate,shearMod)
mwAcc = moMag(accMo)
print mwAcc

pP.plotPatch(latCenter,lonCenter,gLon,gLat,dxGPS,dyGPS,e,n,slip)
#plt.figure()
#pD.plotGPS(gLon,gLat,dxGPS,dyGPS,latCenter,lonCenter)
plt.show()
