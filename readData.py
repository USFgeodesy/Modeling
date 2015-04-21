import sys #for getting comand line data 
import numpy as np
from math import sin, cos, sqrt, atan2, radians
from geopy.distance import great_circle
import okada_functions as fun
'''
USF Geodesy Lab
Nick Voss
Functions for reading Fault, and Slip files for OkadaPy
'''

def distance(lat1,lon1,lat2,lon2):
    R = 6373.0

    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)

    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (np.sin(dlat/2))**2 + np.cos(lat1) * np.cos(lat2) * (np.sin(dlon/2))**2
    c = np.degrees(2 * np.arctan2(np.sqrt(a), np.sqrt(1-a)))
    distance = c 
   # print distance 
    dLon = lon2 - lon1
    y = np.sin(dLon) * np.cos(lat2)
    x = np.cos(lat1) * np.sin(lat2) \
        - np.sin(lat1) * np.cos(lat2) * np.cos(dLon)
    azimuth = np.degrees(np.arctan2(y, x))
    return distance,azimuth
  
def readFault(filename,nh,nv):
  #this is for if you have just a x,y,z file it is not working right now 
    x =np.loadtxt(filename)
    x = x.transpose()
    Lo1 = x[0]
    La1 = x[1]
    z = x[2]
    Lo1 = np.array(Lo1)
    La1 = np.array(La1)
    #print len(La1)
    z = np.array(z)
    #print z 
    # this is due to the file format of Yan it is the Length and width of the grid
    #nh=19
    #nv=17
    # number of km per degree
    dist=(6371.0*np.pi)/180.0
    LatCF = []
    LonCF = []
    zCF = []
    for i in range(nh-1):
        for j in range(nv-1):
            ii=(j)*nh+i
            LatCF.append(0.25*(La1[ii]+La1[ii+1]+La1[ii+nh]+La1[ii+1+nh]))
            LonCF.append(0.25*(Lo1[ii]+Lo1[ii+1]+Lo1[ii+nh]+Lo1[ii+1+nh]))
            zCF.append(0.25*(z[ii]+z[ii+nh]+z[ii+1]+z[ii+1+nh]))
    #print LatCF
        
    L = []
    Strike  = []
    W = []
    nY = nh-1
    nX = nv-1
    xx = np.zeros((nX,nY))
    yy = np.zeros((nX,nY))
    zz = np.zeros((nX,nY))
    Strike1 = np.zeros((nX,nY))
    Dip2 = np.zeros((nY,nX))
    L3 = np.zeros((nX,nY))
    W4 = np.zeros((nX,nY))
    Dip = []
    k = 0
    Len,Width = [],[]
    for i in range(nh-1):
        for j in range(nv-1):
            ii=(j)*nh+i
            l1,az1=distance(La1[ii],Lo1[ii],La1[ii+1],Lo1[ii+1])
            #print l1
            l2,az2=distance(La1[ii+nh],Lo1[ii+nh],La1[ii+1+nh],Lo1[ii+1+nh])
            L.append(0.5*dist*(l1+l2))
            Len.append(0.5*dist*(l1+l2))
            Strike.append((0.5*(az1+az2))+360.0)
            l1,az1=distance(La1[ii],Lo1[ii],La1[ii+nh],Lo1[ii+nh])
            #print (l1*111.1)
            l2,az2=distance(La1[ii+1],Lo1[ii+1],La1[ii+1+nh],Lo1[ii+1+nh])
            Y=(0.5*dist*(l1+l2))
            #print Y
            Z=0.5*(z[ii]-z[ii+nh]+z[ii+1]-z[ii+1+nh])
            Width.append(Z)
            Dip.append(((atan2(Z,Y)))*180/np.pi)
           # print (((atan2(Z,Y)))*180/np.pi)
            W.append(sqrt(Y*Y+Z*Z))
    widthMat = np.zeros((nX,nY))
    lenMat = np.zeros((nX,nY))
    for i in range(nX):
        for j in range(nY):
            xx[i,j]=LonCF[k]
            yy[i,j]=LatCF[k]
            zz[i,j]=zCF[k]
            Strike1[i,j]=Strike[k]
            #Dip2[i,j]=Dip[k]
            L3[i,j]=L[k]
            W4[i,j]=W[k]
            widthMat[i,j] = Width[k]
            lenMat[i,j] = Len[k]
            k = k + 1
    k = 0
    for i in range(nY):
        for j in range(nX):
            Dip2[i,j]=Dip[k]
            k = k + 1
    #print zz
    startLat = np.zeros((nX,nY))
    startLon = np.zeros((nX,nY))
    finLon = np.zeros((nX,nY))
    finLat = np.zeros((nX,nY))
    k = 0
    for i in range(nX):
        for j in range(nY):
           startLat[i,j] = La1[k] 
           
           startLon[i,j] = Lo1[k]
           k = k+1
        k = k+1
    k = 0
    for i in range(nX):
        k = k+1
        for j in range(nY):
           finLat[i,j] = La1[k]
           #print i,k
           #print finLat[i,j], i 
           finLon[i,j] = Lo1[k]
           k = k+1
   # print finLat
    #print startLon
        
# find center of grid i.e km 0,0
    lonCenter = -85.3470000
    latCenter = 10.1199000
    startX = np.zeros((nX,nY))
    startY = np.zeros((nX,nY))
    finX = np.zeros((nX,nY))
    finY = np.zeros((nX,nY))
    k = 0
# computes the x,y coordinates with repect to the center of the grid in km
    for i in range(nX):
        for j in range(nY):
           startX[i,j] = (startLon[i,j]-lonCenter)*111.1
           startY[i,j] = (startLat[i,j]-latCenter)*111.1
           #print startY[i,j], i 
           k = k+1
        #k = k+1
    k = 0
    for i in range(nX):
        for j in range(nY):
           finX[i,j] = (finLon[i,j]-[lonCenter])*111.1
           finY[i,j] = (finLat[i,j]-[latCenter])*111.1
           #print finY[i,j]
           #print finY[i,j],i
           
           k = k+1
        #k = k+1   
    #print finY
#computes top and bottom of each slab
    k = 0
    faultTop = np.zeros((nX,nY))
    faultBottom = np.zeros((nX,nY)) 
    for i in range(nX):
        for j in range(nY):
           faultTop[i,j] = (z[k])
           k = k+1
        k = k+1
        
    k = nh
    for i in range(nX):
        for j in range(nY):
           faultBottom[i,j] = (z[k])
           k = k+1       
        k = k+1
    #print startLat
    #print np.shape(finY)
    Dip2 = Dip2.transpose()
    file = open('faultData.txt','w')
    file.write("x y  depth      dip  strike    L       W \n ")
    for i in range(len(Strike1)):
      for j in range(len(Strike1[0])):
	x = (startX[i,j]+finX[i,j])/2.0
	y = (startY[i,j]+finY[i,j])/2.0
	z = (faultTop[i,j]+faultBottom[i,j])/2.0
	file.write(str(x)+ '  ' + str(y) + '  '+ str(z)+'  '+str(Dip2[i,j]) + '  ' + str(Strike1[i,j]) + '  ' + str(L3[i,j])+ '  ' + str(W4[i,j]) + "\n")
    file.close()
    #return latCenter,lonCenter 
    return (Strike1,Dip2,L3,W4,startX,startY,finX,finY,latCenter,lonCenter,faultTop,faultBottom)
def readDataFile(filename):
   data = np.loadtxt(filename,skiprows = 1)
   data = data.transpose()
   #lon lat  depth      dip  strike    L       W 
   lon = data[0]
   lat = data[1]
   lonCenter = np.mean(lon)
   latCenter = np.mean(lat)
   center = (latCenter,lonCenter)
   x,y = [],[]
   #calculate unit vector around lat and lon center
   for i in range(len(lon)):
    e,n,u = fun.lalo2xy(lat[i],lon[i],latCenter,lonCenter)
    x.append(e)
    y.append(n)
   z = data[2]
   dip = data[3]
   strike = data[4]
   L = data[5]
   W = data[6]
   return x,y,z,dip,strike,L,W ,latCenter,lonCenter,lat,lon
  
def readGPS(filename,latCenter,lonCenter):
    '''
    reads in the GPS lat, lon, and z height from text file
    format of text file is longitude latitude altitude
    input = .txt
    output = gps Lat gps Lon gps Height 
    '''
    x = np.loadtxt(filename,usecols = (1,2,3))
    x = x.transpose()
    gLo1 = x[0]
    gLa1 = x[1] 
    xLoc,yLoc = [],[]
    #convert to local reference frame 
    for i in range(len(gLo1)):
      e,n,u = fun.lalo2xy(gLa1[i],gLo1[i],latCenter,lonCenter)
      xLoc.append(e)
      yLoc.append(n)
    gz = x[2]
    return xLoc,yLoc,gLo1,gLa1,gz
  
def readSlip(filename):
  
    '''
    reads in the Faultindex,  ss , ds, t   from text file
    format of text file is longitude latitude altitude stikeslip dipslip tensile
    input = .txt
    output = strike slip dipslip tensile 
    '''
    
    x = np.loadtxt(filename)
    x = x.transpose()
    index = x[0]
    ss = x[1]
    ds = x[2]
    t = x[3]
    return index, ss, ds, t 
    
