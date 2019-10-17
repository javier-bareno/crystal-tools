'''
Created on Aug 12, 2009

@author: bareno
'''
from numpy import *
#import scipy
import math

def v_angle(a,b):
    """
    returns angle between two vectors coded as numpy.array
    """
    ab=dot(a,b)
    ab/=linalg.norm(a)*linalg.norm(b)
    return math.degrees(math.acos(ab))


class Basis:
    '''
    Class to handle crystallographic basis and basis conversions
    '''


    def __init__(self, a=1., b=1., c=1., alpha=90., beta=90., gamma=90.):
        '''
        Basis constructor from standard description.
        Generates and stores reciprocal lattice basis.
        Matrix M to convert to cartesian coords. Mr in reciprocal space.
        '''
        
        bx =math.cos(math.radians(gamma))
        by = math.sin(math.radians(gamma))
        cx = math.cos(math.radians(beta))
        cy = (math.cos(math.radians(alpha)) - \
            math.cos(math.radians(beta))*math.cos(math.radians(gamma))) \
            / math.sin(math.radians(gamma))        
        cz = math.sqrt(1-cx**2-cy**2)
        self.av=array([a,0.,0.])
        self.bv=array([b*bx,b*by,0.])
        self.cv=array([c*cx,c*cy,c*cz])
        self.a=a
        self.b=b
        self.c=c
        self.alpha=alpha
        self.beta=beta
        self.gamma=gamma
        self.vol=dot(self.av,cross(self.bv,self.cv))
        self.M = array([self.av, self.bv, self.cv]).transpose()
        self.iM = linalg.inv(self.M)
        self.avr = cross(self.bv,self.cv)/self.vol
        self.bvr = cross(self.cv,self.av)/self.vol
        self.cvr = cross(self.av,self.bv)/self.vol
        self.ar = linalg.norm(self.avr)
        self.br = linalg.norm(self.bvr)
        self.cr = linalg.norm(self.cvr)
        self.alphar = v_angle(self.bvr,self.cvr)
        self.betar = v_angle(self.cvr,self.avr)
        self.gammar = v_angle(self.bvr,self.avr)
        self.volr = 1./self.vol
        self.Mr = array([self.avr, self.bvr, self.cvr]).transpose()
        self.iMr = linalg.inv(self.Mr)
        
    def __cmp__(self, other):
        if not (self.a == other.a):
            return False
        elif not (self.b == other.b):
            return False
        elif not (self.c == other.c):
            return False
        elif not (self.alpha == other.alpha):
            return False
        elif not (self.beta == other.beta):
            return False
        elif not (self.gamma == other.gamma):
            return False        
        else:
            return(True)
    
    def __repr__(self):
        return (self.M)
    
    def __str__(self):
        ret_str = "Direct basis:\n"
        ret_str += "a= %.4f\tb= %.4f\tc= %.4f\n" %(self.a, self.b, self.c)
        ret_str += "alpha= %.2f\tbeta= %.2f\tgamma= %.2f\n" %(self.alpha, self.beta, self.gamma)
        ret_str += "vol= %.4f\n\n" %(self.vol)
        ret_str += "Reciprocal basis:\n"
        ret_str += "a= %.4f\tb= %.4f\tc= %.4f\n" %(self.ar, self.br, self.cr)
        ret_str += "alpha= %.2f\tbeta= %.2f\tgamma= %.2f\n" %(self.alphar, self.betar, self.gammar)
        ret_str += "vol= %.4f\n\n" %(self.volr)
        return ret_str
    
    def from_basis_vectors(self,av,bv,cv):
        '''
        Creates a Basis from three basis vectors.
        Specified as numpy.array's
        '''
        a = linalg.norm(av)
        b = linalg.norm(bv)
        c = linalg.norm(cv)
        alpha = v_angle(bv, cv)
        beta = v_angle(av, cv)
        gamma = v_angle(av, bv)
        ret_basis = Basis(a,b,c, alpha, beta, gamma)
        ret_basis.av=av
        ret_basis.bv=bv
        ret_basis.cv=cv        
        ret_basis.vol=dot(ret_basis.av,cross(ret_basis.bv,ret_basis.cv))
        ret_basis.M = array([ret_basis.av, ret_basis.bv, ret_basis.cv]).transpose()
        ret_basis.iM = linalg.inv(ret_basis.M)
        ret_basis.avr = cross(ret_basis.bv,ret_basis.cv)/ret_basis.vol
        ret_basis.bvr = cross(ret_basis.cv,ret_basis.av)/ret_basis.vol
        ret_basis.cvr = cross(ret_basis.av,ret_basis.bv)/ret_basis.vol
        ret_basis.ar = linalg.norm(ret_basis.avr)
        ret_basis.br = linalg.norm(ret_basis.bvr)
        ret_basis.cr = linalg.norm(ret_basis.cvr)
        ret_basis.alphar = v_angle(ret_basis.bvr,ret_basis.cvr)
        ret_basis.betar = v_angle(ret_basis.cvr,ret_basis.avr)
        ret_basis.gammar = v_angle(ret_basis.bvr,ret_basis.avr)
        ret_basis.volr = 1./ret_basis.vol
        ret_basis.Mr = array([ret_basis.avr, ret_basis.bvr, ret_basis.cvr]).transpose()
        ret_basis.iMr = linalg.inv(ret_basis.Mr)
        return ret_basis
        
    
    def extend(self, xa=1, xb=1, xc=1):
        '''
        extends unit cell along basis vectors times xa, xb, xc
        '''
        retBasis = Basis(xa*self.a, xb*self.b, xc*self.c, self.alpha, self.beta, self.gamma)
        retBasis.av = xa * self.av
        retBasis.bv = xb * self.bv
        retBasis.cv = xc * self.cv        
        retBasis.vol=dot(retBasis.av,cross(retBasis.bv,retBasis.cv))
        retBasis.M = array([retBasis.av, retBasis.bv, retBasis.cv]).transpose()
        retBasis.iM = linalg.inv(retBasis.M)
        retBasis.avr = cross(retBasis.bv,retBasis.cv)/retBasis.vol
        retBasis.bvr = cross(retBasis.cv,retBasis.av)/retBasis.vol
        retBasis.cvr = cross(retBasis.av,retBasis.bv)/retBasis.vol
        retBasis.ar = linalg.norm(retBasis.avr)
        retBasis.br = linalg.norm(retBasis.bvr)
        retBasis.cr = linalg.norm(retBasis.cvr)
        retBasis.alphar = v_angle(retBasis.bvr,retBasis.cvr)
        retBasis.betar = v_angle(retBasis.cvr,retBasis.avr)
        retBasis.gammar = v_angle(retBasis.bvr,retBasis.avr)
        retBasis.volr = 1./retBasis.vol
        retBasis.Mr = array([retBasis.avr, retBasis.bvr, retBasis.cvr]).transpose()
        retBasis.iMr = linalg.inv(retBasis.Mr)
        return retBasis       
            
    def dir_to_cart(self, v):
        return dot(self.M,v)
    
    def dir_from_cart(self, v):
        return dot(self.iM,v)
    
    def plane_to_cart(self,v):
        return dot(self.Mr,v)
    
    def plane_from_cart(self,v):
        return dot(self.iMr, v)
    
    def dir_from_plane(self, v):
        return self.dir_from_cart(self.plane_to_cart(v))
    
    def plane_from_dir(self, v):
        return self.plane_from_cart(self.dir_to_cart(v))
    
    def dir_norm(self, v):
        return linalg.norm(self.dir_to_cart(v))
    
    def plane_dist(self, p):
        return 1./linalg.norm(self.plane_to_cart(p))
    
    def dir_dot(self, d1, d2):
        return dot(self.dir_to_cart(d1), self.dir_to_cart(d2) )
    
    def plane_dot(self, p1, p2):
        return dot(self.plane_to_cart(p1), self.plane_to_cart(p2))
    
    def mix_dot(self, p, v):
        return dot(p,v)
    
    def plane_angle(self, p1, p2):
        dot = self.plane_dot(p1, p2)
        cosang = dot * self.plane_dist(p1) * self.plane_dist(p2)
        if abs(cosang)>1:
            return 0.
        else:
            return math.degrees(math.acos(cosang))
        
    
    def dir_angle(self, v1, v2):
        dot = self.dir_dot(v1, v2)
        cosang = dot / (self.dir_norm(v1)* self.dir_norm(v2))
        if abs(cosang)>1:
            return 0.
        else:
            return math.degrees(math.acos(cosang))
        
    
    def common_plane(self, d1, d2):
        cart_plane = cross(self.dir_to_cart(d1),self.dir_to_cart(d2))
        return self.plane_from_cart(cart_plane)
    
    def common_dir(self, p1, p2):
        cart_dir = cross(self.plane_to_cart(p1),self.plane_to_cart(p2))
        return self.dir_from_cart(cart_dir)
    
    def new_basis(self, G):
        """ 
        Creates a new base from matrix G.
        Columns of G taken as coords of new basis vectors in self basis.
        """
        nB = Basis(1,1,1)
        nB.av = self.dir_to_cart(G.transpose()[0])
        nB.bv = self.dir_to_cart(G.transpose()[1])
        nB.cv = self.dir_to_cart(G.transpose()[2])
        nB.M = array([nB.av, nB.bv, nB.cv]).transpose()
        nB.iM = linalg.inv(nB.M)
        nB.a = linalg.norm(nB.av)
        nB.b = linalg.norm(nB.bv)
        nB.c = linalg.norm(nB.cv)
        nB.alpha = v_angle(nB.bv, nB.cv)
        nB.beta = v_angle(nB.cv, nB.av)
        nB.gamma = v_angle(nB.av, nB.bv)
        nB.vol = dot(nB.av, cross(nB.bv, nB.cv))
        nB.avr = cross(nB.bv, nB.cv)/nB.vol
        nB.bvr = cross(nB.cv, nB.av)/nB.vol
        nB.cvr = cross(nB.av, nB.bv)/nB.vol
        nB.ar = linalg.norm(nB.avr)
        nB.br = linalg.norm(nB.bvr)
        nB.cr = linalg.norm(nB.cvr)
        nB.alphar = v_angle(nB.bvr, nB.cvr)
        nB.betar = v_angle(nB.avr, nB.cvr)
        nB.gammar = v_angle(nB.avr, nB.bvr)
        nB.Mr = array([nB.avr, nB.bvr, nB.cvr]).transpose()
        nB.iMr = linalg.inv(nB.Mr)
        nB.volr = dot(nB.avr, cross(nB.bvr, nB.cvr))
        return nB
          
    def seek_angle(self, target, tol, rg=range(-3,4)):
        s=[]
        for h1 in rg:
            for k1 in rg:
                for l1 in rg:
                    for h2 in rg:
                        for k2 in rg:
                            for l2 in rg:
                                x = self.plane_angle(array([h1,k1,l1]), array([h2,k2,l2]))
                                if abs(x-target)< tol:
                                    st1 ="(%d, %d, %d)\t" %(h1,k1,l1)
                                    st1 += "(%d, %d, %d)\t%.2f" %(h2,k2,l2,x)
                                    s.append(st1)
        return s
    
    def update_basis(self, **update_dir):
        '''
        Auxiliary for XRD 2theta optimization
        Gnerates new cell with modified lattice params
        Needs to be called with named arguments
        If param not in arg_list, use same value as self'''
        cell_params = {}
        for key in ['a', 'b', 'c', 'alpha', 'beta', 'gamma']:
            cell_params[key]= getattr(self,key)
            if key in update_dir:
                #this is a guess param
                cell_params[key]=update_dir[key]
                
        
        return Basis(**cell_params)
                                
    
if __name__ == '__main__':
    A=Basis(9.8028,16.9226,4.9885,90,109.09,90)
    #m=array([[0.5,0,0],[-0.25,0.25,0],[0,0,1]]).transpose()
    #P=A.new_basis(m)
    rgh = range(-4,3)
    rgk = range(-4,3)
    rgl=range(-4,3)
    for h1 in rgh:
        for k1 in rgk:
            for l1 in rgk:
                for h2 in rgh:
                    for k2 in rgk:
                        for l2 in rgk:
                            x1=array([h1,k1,l1])
                            x2 = array([h2,k2,l2])
                            ang =A.plane_angle(x1, x2)
                            if abs(ang - 63) < 1:
                                rat =(A.plane_dist(x1) / A.plane_dist(x2))
                                if abs(rat - 1) < 0.1:
                                    #z = cross(x1,x2)
                                    #x1 = A.plane_from_cart(P.plane_to_cart(x1))
                                    #x2 = A.plane_from_cart(P.plane_to_cart(x2))
                                    #z= A.dir_from_cart(P.dir_to_cart(z))
                                    z = cross(x1,x2)
                                    s= "%.2f\t%.3f\t" %(ang, rat)
                                    s+= "(%.3f, %.3f, %.3f)\t" %(x1[0], x1[1], x1[2])
                                    s+= "(%.3f, %.3f, %.3f)\t" %(x2[0], x2[1], x2[2])
                                    s+= "[%.3f, %.3f, %.3f]\t" %(z[0], z[1], z[2])
                                    print s
    print "Done :)"
    
                
        