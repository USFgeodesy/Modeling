import numpy as np 
import dC3D as dC

def okada(alpha,x,y,z,depth,dip,al1,al2,aw1,aw2,ssDis,dsDis,tDis):
    ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret = dC.dc3d(alpha,x,y,z,depth,dip,al1,al2,aw1,aw2,ssDis,dsDis,tDis)
    return ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret
