import numpy as np 

'''
Controls outputs from okadaPY
USF Geodesy Lab Nick Voss
Written Oct 3, 2014
'''

def outputGrid(gridX,gridY, dx,dy,dz):
  '''
  outputs the displacement on the grid points calculated by okadaPY
  to a text file names outputGrid.txt
  '''
  file = open("outputGrid.txt", "w")  #output file 
  file.write("Gridx  Gridy  dx  dy dz\n")
  k = 0 
  for i in range(len(gridX)):
    for j in range(len(gridY)):
      file.write(str(gridX[i]) + "  "  +  str(gridY[j]) + "  " + str(dx[k]) + "  "  +  str(dy[k]) + " " + str(dz[k]) + "\n")
      k  = k+1
  file.close()
  
def outputGPS(gpsX,gpsY,dx,dy,dz):
  '''
  outputs the displacement at the GPS point calculated by okadaPY
  to a text file names outputGPS.txt
  '''
  file = open("/home/nvoss/Geodesy/PEST_OKADAPY/outputGPS.txt", "w")  #output file 
  file.write("GPSx  GPSy  dx  dy dz\n")
  for i in range(len(gpsX)):
    file.write(str(gpsX[i]) + "  "  +  str(gpsY[i]) + " " + str(dx[i])+" "  +  str(dy[i]) + " " + str(dz[i]) + "\n")
  file.close()
	       