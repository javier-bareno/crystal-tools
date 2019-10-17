'''
Created on Jan 25, 2010

@author: bareno
'''
from numpy import *
import math
from Structure import Structure
from Base import Basis
from scipy import optimize
import profile

class XRD_pattern:
    '''
    class to handle cell refinement from known (indexed) XRD peaks'''
    
    def __init__(self, lattice =Basis(), pattern=[], wavelength = 1.54056):
        '''
        Constructor
        wavelength defaults to Cu K alpha
        pattern is list of peaks, each a tuple (h, k, l, 2theta)'''
        self.lattice = lattice
        self.pattern = pattern
        self.wavelength = wavelength
    
    def error(self, lattice = 0):
        '''
        returns residual of peaks using lattice
        if no lattice especified, uses self.lattice '''
        if type(lattice) != type(Basis()):
            lattice =self.lattice
            
        residual = 0
        for peak in self.pattern:
            indices = array(peak[0:3])
            error = lattice.plane_dist(indices)
            error = 2*degrees(arcsin( 0.5*self.wavelength / error ))
            error -= peak[3]
            residual += error**2
        return residual
    
    def powell(self, *params):
        ''' 
        returns an XRD_pattern object with a new cell with optimized lattice constants params
        Non optimized and guesses taken from self.lattice
        params is list of lattice constant names to optimize'''
        
        def guess_Basis(*x0):
            lattice_params ={}
            for key in ['a','b','c','alpha','beta','gamma']:
#                print key
                if key in params:
                    which_param = params.index(key)
                    lattice_params[key]=x0[which_param]
                else:
                    lattice_params[key] = getattr(self.lattice,key)
            
#            print lattice_params
            return Basis(**lattice_params)
            
        def func(x0):           
            kk=guess_Basis(*x0)
#            print kk
#            print x0, self.error(kk)
            return self.error(kk)
        
        guess=[]
        lower =[]
        upper =[]
        for param in params:
            guess.append(getattr(self.lattice,param))
            if param in ['a','b','c']:
                lower.append(1.0)
                upper.append(30.0)
            else:
                lower.append(10.0)
                upper.append(170.0)
               
        xf = optimize.fmin_powell(func, guess)
        
        return XRD_pattern(guess_Basis(*xf), self.pattern, self.wavelength)
    
    def fmin(self, *params):
        ''' 
        returns an XRD_pattern object with a new cell with optimized lattice constants params
        Non optimized and guesses taken from self.lattice
        params is list of lattice constant names to optimize'''
        
        def guess_Basis(*x0):
            lattice_params ={}
            for key in ['a','b','c','alpha','beta','gamma']:
#                print key
                if key in params:
                    which_param = params.index(key)
                    lattice_params[key]=x0[which_param]
                else:
                    lattice_params[key] = getattr(self.lattice,key)
            
#            print lattice_params
            return Basis(**lattice_params)
            
        def func(x0):           
            kk=guess_Basis(*x0)
#            print kk
#            print x0, self.error(kk)
            return self.error(kk)
        
        guess=[]
        for param in params:
            guess.append(getattr(self.lattice,param))
               
        xf = optimize.fmin(func, guess)
        
        return XRD_pattern(guess_Basis(*xf), self.pattern, self.wavelength)

    def twotheta(self,h,k,l,wavelength = 1.54056):
        ref = array([h,k,l])
        d=self.lattice.plane_dist(ref)
        sinth = 0.5*wavelength/d
        return 2*degrees(arcsin(sinth)) 
                    

if __name__ == '__main__':
    def test():
        mono_cell = Basis(5.0, 9.0, 4.5, beta=90.0)
        mono_pattern=[(0,0,1,18.92), (0,2,0,21.036), (1,1,0,21.964), (-1,1,1,24.456)]
        model = XRD_pattern(mono_cell, mono_pattern)
        A=model.fmin('a','b','c','beta')
        B=model.powell('a','b','c','beta')
        print 'Starting guess'
        print 'a=%.5f b=%.5f c=%.5f beta=%.5f' %(model.lattice.a, model.lattice.b, model.lattice.c, model.lattice.beta)
        print 'Residual: %.5f' %(model.error()**.5)
        print 'fmin model:'
        print 'a=%.5f b=%.5f c=%.5f beta=%.5f' %(A.lattice.a, A.lattice.b, A.lattice.c, A.lattice.beta)
        print 'Residual: %.5f' %(A.error()**.5)
        print 'Powell model'
        print 'a=%.5f b=%.5f c=%.5f beta=%.5f' %(B.lattice.a, B.lattice.b, B.lattice.c, B.lattice.beta)
        print 'Residual: %.5f' %(B.error()**.5)
#        print (1,3,-2), A.twotheta(1,3,-2)
#        print (1,1,2), A.twotheta(1,1,2)
#        print (2,2,-2), A.twotheta(2,2,-2)
#        print (-2,2,2), A.twotheta(-2,2,2)
#        print (2,2,1), A.twotheta(2,2,1)
    
    profile.run('test()')
    
    
    