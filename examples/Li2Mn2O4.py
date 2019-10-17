'''
Created on Aug 22, 2009

@author: javierbareno
'''
from numpy import *
import math
from Structure import Structure

if __name__ == '__main__':
    ifname="Li2Mn2O4-spinel.ACE"
    ofname="Li2Mn2O4-spinel"
    
    print "Reading input file"
    ifile = open(ifname)
    Li2Mn2O4=Structure()
    Li2Mn2O4.read_cell(ifile)
    ifile.close()
    ofile = open(ofname + ".xyz", 'w')
    Li2Mn2O4.write_xyz(ofile, "Converted from Carine cell file")
    ofile.close()
    #G=array([[-1,0,1],[0,1,0],[2.65168315,0,1]]).transpose()
    
    pov=Li2Mn2O4.extend(1.01, 1.01, 1.01)
    ofile = open(ofname + ".pov",'w')
    pov.write_povray(ofile)
    ofile.close()
    
    #print "Changing to rhomb base"
    big = Li2Mn2O4.extend(3,3,3)
    
    
    print "Writing rhomb output"
    ofile= open(ofname +"-big.dat",'w')
    big.write(ofile)#, "Converted to hex base from standard cubic spinel def")
    ofile.close()
    ofile= open(ofname +"-big.xyz",'w')
    big.write_xyz(ofile,"Converted to hex base from standard cubic spinel def")
    ofile.close()
    """ofile = open(ofname+"-rhom.emp",'w')
    rhom.write_emaps(ofile)
    ofile.close()"""
    
    """
    rhom2 = rhom.extend(1.01, 1.01, 1.01)
    ofile= open(ofname +"-rhom.pov",'w')
    rhom2.write_povray(ofile, "Converted to hex base from standard cubic spinel def")
    ofile.close()"""
    
    print "splitting into plane stack"
    plane_stack = big.cluster(array([1,0,1]), 1./30)
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