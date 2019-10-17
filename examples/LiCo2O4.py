'''
Created on Aug 22, 2009

@author: javierbareno
'''
from numpy import *
import math
from Structure import Structure

if __name__ == '__main__':
    ifname="LiCo2O4-spinel.ACE"
    ofname="LiCo2O4-spinel"
    
    print "Reading input file"
    ifile = open(ifname)
    Li2Mn2O4=Structure()
    Li2Mn2O4.read_cell(ifile)
    ifile.close()
    ofile = open(ofname + ".xyz", 'w')
    Li2Mn2O4.write_xyz(ofile, "Converted from Carine cell file")
    ofile.close()
    G=array([[0.5,-1,0.5],[0.5,0,-0.5],[1,1,1]]).transpose()
    
    print "Changing to rhomb base"
    rhom = Li2Mn2O4.new_basis(G)
    a,b,c = rhom.basis.a, rhom.basis.b, rhom.basis.c
    alpha, beta, gamma = rhom.basis.alpha, rhom.basis.beta, rhom.basis.gamma
    rhom2 = Structure(a,b,c,alpha,beta,gamma)
    rhom.basis = rhom2.basis
    
    print "Writing rhomb output"
    ofile= open(ofname +"-rhom.dat",'w')
    rhom.write(ofile)#, "Converted to hex base from standard cubic spinel def")
    ofile.close()
    ofile= open(ofname +"-rhom.xyz",'w')
    rhom.write_xyz(ofile,"Converted to hex base from standard cubic spinel def")
    ofile.close()
    ofile = open(ofname+"-rhom.emp",'w')
    rhom.write_emaps(ofile)
    ofile.close()
    
    
    rhom2 = rhom.extend(1.01, 1.01, 1.01)
    ofile= open(ofname +"-rhom.pov",'w')
    rhom2.write_povray(ofile, "Converted to hex base from standard cubic spinel def")
    ofile.close()
    
    print "splitting into plane stack"
    plane_stack = rhom2.cluster(array([0,0,1]), 2./49)
    i=0
    
    print "writing stack output"
    for plane in plane_stack:
        ofile = open(ofname + "-rhom-plane-%02d.dat" %(i),'w')
        plane.write(ofile)
        ofile.close()
        ofile = open(ofname + "-rhom-plane-%02d.xyz" %(i),'w')
        plane.write_xyz(ofile)
        ofile.close()        
        ofile = open(ofname + "-rhom-plane-%02d.pov" %(i),'w')
        plane.write_povray(ofile)
        ofile.close()
        i+=1
        
       
    print "That's all folks"