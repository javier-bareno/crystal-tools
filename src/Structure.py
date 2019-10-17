'''
Created on Aug 21, 2009

@author: bareno
'''

from numpy import *
#import scipy
import math
from Base import Basis

class Structure():
    """Store, read, and write unit cell info.
    All units in Angstrom """
    def __init__(self, a=1., b=1., c=1., alpha=90., beta=90., gamma=90.):
        self.basis= Basis( a, b, c, alpha, beta, gamma)        
        self.num_atoms=0
        self.atoms=[]
        self.atom_properties =[] # acts as storage of other props user may want to keep track of
        self.center=array([0.5,0.5,0.5])
    

    def __len__(self): return(self.num_atoms)
    
    def __cmp__(self, other):
        """ cmp method compares Bases for compatibility """
        return self.basis == other.basis
    
    def __getitem__(self, num):
        return self.atoms[num]
    
    def __repr__(self):
        return [self.basis.__repr__(), self.num_atoms, self.atoms, self.atom_properties, self.center]
    
    def __str__(self): 
        retstr = "a=%.3f, b=%.3f, c=%.3f\n" %(self.basis.a, self.basis.b, self.basis.c)
        retstr+= "alpha=%.3f, beta=%.3f, gamma=%.3f\n" %(self.basis.alpha, self.basis.beta, self.basis.gamma)
        retstr+= "%d atoms in the structure" %(self.num_atoms)       
        return retstr

    def from_basis(self, basis):
        '''
        Creates empty structure from basis
        '''
        retstruct = UnitCell(basis.a, basis.b, basis.c, basis.alpha, basis.beta, basis.gamma)
        retstruct.basis = basis
        return retstruct
        
        
    def adatom(self, symbol, x, y, z):
        self.num_atoms +=1
        self.atoms.append((symbol, array([x, y, z])))
        
    def write(self, ofile):
        #print "write called with ofile:"
        try:
            ofile.write('Unit cell\n')
            ofile.write('Unit vectors (Angstroms)\n')
            ofile.write('%f\t%f\t%f\n' % tuple(self.basis.av))
            ofile.write('%f\t%f\t%f\n' % tuple(self.basis.bv))
            ofile.write('%f\t%f\t%f\n' % tuple(self.basis.cv))
            ofile.write(str(self.num_atoms) + ' atoms\n')
            for atom in self.atoms:
                ofile.write("%s\t" %(atom[0]))
                ofile.write("%f\t%f\t%f\n" % tuple(atom[1]))
            if len(self.atom_properties)!=0:
                for prop in self.atom_properties:
                    ofile.write(str(prop))
        except IOError:
            raise

    def read(self, ifile):
        #To do cleanup
        def str_to_float(st):
            st2 = st.split("/")
            if len(st2)==2:
                return float(st2[0])/float(st2[1])
            else:
                return float(st)
            
        def str_to_array(st):
            """ Aux funtion to transform "%f %f %f" into array of floats.
            Understands fractions, e.g. if %f ="1/4"""            
            v = []
            for st2 in st.split():
                v.append(str_to_float(st2))
            return array(v)
            
        try:
            ifile.readline()
            ifile.readline()
            av= str_to_array(ifile.readline())
            bv= str_to_array(ifile.readline())
            cv= str_to_array(ifile.readline())
            self.basis = Basis().from_basis_vectors(av,bv,cv)
            self.num_atoms = int(ifile.readline().split()[0])
            for atn in range(self.num_atoms):
                atn2 = ifile.readline().split()
                atred = str_to_array(atn2[1:])
                self.atoms.append((atn2[0], atred))
            #To do: add code to read/write atom_properties
            prop = ifile.readline()
            while prop != '':
                self.atom_properties.append(prop)
                prop = ifile.readline()
            ifile.close()
            

        except IOError:
            raise
    
    def write_povray(self, ofile, comment="Blank comment"):
        """Writes atoms in PovRay #include file"""
        
        try:
            ofile.write("//" + comment +"\n")
            for atom in self.atoms:
                ofile.write("PutAtom_" + atom[0] )
                atwc= dot(self.basis.M,(atom[1] - self.center))
                atwc*=array([-1,1,1])
                ofile.write("(%f, %f, %f)\n" %tuple(atwc))
        except IOError:
            print "Inner exception in write_povray"
            raise
        
    def write_xyz(self, ofile, comment = "\n"):
        try:
            ofile.write("%d\n" %(self.num_atoms))
            ofile.write(comment)
            if comment[-1] != "\n":
                ofile.write("\n")
            for atom in self.atoms:
                atwc= dot(self.basis.M,(atom[1] - self.center))
                ofile.write("%s\t" %(atom[0]) +"%.3f\t%.3f\t%.3f\n" %tuple(atwc))
        except IOError:
            raise

    def extend(self,xa=1,xb=1,xc=1):
        """Extends unit cell along unit vectors xa, xb, xc"""
        retUnitCell = Structure()
        retUnitCell.basis = self.basis.extend(xa,xb,xc)
        retUnitCell.atom_properties = self.atom_properties
        for xxa in range(int(math.ceil(xa))):
            for xxb in range(int(math.ceil(xb))):
                for xxc in range(int(math.ceil(xc))):
                    for atom in self.atoms:
                        new_at_pos = array([1.,1.,1.]) *atom[1] + array([xxa,xxb,xxc])
                        new_at_pos /= array([xa,xb,xc])
                        if new_at_pos[0] <=1:
                            if new_at_pos[1] <=1:
                                if new_at_pos[2] <=1:
                                    retUnitCell.adatom(atom[0], new_at_pos[0], new_at_pos[1], new_at_pos[2])
        return retUnitCell

    def translate(self,xa=0, xb=0, xc=0):
        """Displaces unit cell (xa, xb, xc) reduced  coords"""
        retStruct=Structure()
        retStruct.basis =self.basis
        retStruct.atom_properties = self.atom_properties
        retStruct.num_atoms = self.num_atoms
        for atom in self.atoms:
            retStruct.atoms.append((atom[0], atom[1]+array([xa, xb, xc])))
        return retStruct

    def add(self, uc):
        """ Incorporates contents of UnitCell uc into self.
        Does nothing if uc is different shape"""
        if self == uc:
            for at in uc.atoms:
                self.num_atoms +=1
                self.atoms.append(at)
            return self
        else:
            print "Non-compatible strcutures"
                
    def write_emaps(self, ofile, comment="Blank comment"):
        """ Writes atom list to upload (manually) to webEmaps.
        Assumes Occupancy 1 and D-W = 0 for all atoms"""
        try:
            i=0
            ofile.write(comment+'\n')
            for at in self.atoms:
                i+=1
                str ="%d:  EL:%s,  " %(i, at[0])
                str+= "x:%.6f,  y:%.6f,  z:%.6f,  " %tuple(at[1])
                str +="D-W:%.6f,  occ:%.6f\n" %(0.,1.)
                ofile.write(str)
        except IOError:
            raise    
    
    def carve_plane(self, n, p):
        """ Carves out atoms sitting to the positive side of a plane normal to n passing through p.
        n and p are of type array(). Reciprocal and direct space, respectively"""
        new_at = []
        new_cent=array([0.,0.,0.])
        for atom in self.atoms:
            x = atom[1]    
            x-=p
            if dot(x,n) <= 0:
                new_at.append(atom)
                new_cent +=atom[1]
            self.atoms = new_at
            self.num_atoms=len(self.atoms)
            self.center = new_cent/self.num_atoms
        return self
    
    def read_cell(self, ifile):
        """reads cell definition ascii file from Carine into Structure"""
        
        def str_to_float(st):
            """aux to hndle rational number notation"""
            st2 = st.split("/")
            if len(st2)==2:
                return float(st2[0])/float(st2[1])
            else:
                return float(st)
            
        try:
            rl = ifile.readline()
            while rl[:2] != "a=":
                rl = ifile.readline()
            abc=rl.split()
            a=float(abc[0].split("=")[1])
            b=float(abc[1].split("=")[1])
            c=float(abc[2].split("=")[1])
            rl = ifile.readline()
            abc=rl.split()
            alpha = float(abc[0].split("=")[1])
            beta  = float(abc[1].split("=")[1])
            gamma = float(abc[2].split("=")[1])
            self.basis=Basis(a,b,c,alpha,beta,gamma)
            self.num_atoms =0
            self.atoms =[]
            self.atom_properties=[] 
            self.center = array([0.5, 0.5, 0.5])
            rl = ifile.readline()
            while rl[:28] != "Number of positions in Cell=":
                rl = ifile.readline()
            num_atoms = int(rl[28:])                
            while rl[:4] != "Atom":
                rl = ifile.readline()
            for i in range(num_atoms):               
                rl = ifile.readline().split()
                self.adatom(rl[0], str_to_float(rl[2]), str_to_float(rl[3]) \
                    , str_to_float(rl[4]))
                self.atom_properties.append((rl[1],rl[4],rl[5]))
            return self
        except IOError:
            raise 
    
    def inc_atom(self, symbol, pos, tol):
        """ Folds atom into 0<= reduced coords < 1 and
        adds atom to structure if it is farther than tol angstroms from existing atoms
        Takes into account unit cell periodic boundary conditions"""
        #1st fold into uc
        #Note -0%1 can be = 1. Need to check
        pos = pos % array([1,1,1])
        #print pos
        for i in range(3):
            if abs(pos[i]- 1.)<1e-3:
                pos[i]=0
        for atom in self.atoms:
            dist = 2*tol
            for h in range(2):
                for k in range(2):
                    for l in range(2):
                        dist = min(linalg.norm(self.basis.dir_to_cart(atom[1]+array([h,k,l])-pos)), dist)
            if dist < tol:                
                return self
        self.adatom(symbol, pos[0], pos[1], pos[2])
        return self
        
    def new_basis(self,G):
        """ Creates a new basis taken G columns as coords of new basis vectors in current basis.
        Then folds atoms into new basis avoiding duplicates. Cells share origin """
        nst = Structure()
        nst.basis = self.basis.new_basis(G)
        #1st make sure that I have enough atoms to fill new structure
        v_rat = abs(nst.basis.vol / self.basis.vol)
        n_at_target = int("%.0f" %(v_rat * self.num_atoms))
        v_sc = math.ceil(v_rat**0.333)
        temp = self.extend(v_sc+1, v_sc+1, v_sc+1)
        for atom in temp.atoms:
            symbol = atom[0]
            pos = nst.basis.dir_from_cart(temp.basis.dir_to_cart(atom[1]))
            
                           
            nst.inc_atom(symbol, pos, 0.3)
            #if nst.num_atoms == n_at_target:
            #    return nst
        return nst
    
    def sort(self, plane):
        """ Sorts cell along normal to plane, taken as array([]) in reciprocal coords"""
        i=0
        sort_list=[]
        for atom in self.atoms:
            dist = dot(atom[1],plane)
            sort_list.append([i,dist])
            i+=1
        sort_list.sort(lambda x,y: cmp(x[1],y[1]))
        retStr = Structure()
        retStr.basis = self.basis
        for el in sort_list:
            retStr.addatom(self.atoms[el[0]])
        return ret_str
    
    def cluster(self, plane, width):
        """ Returns a list of sub-structures from self.
        Atoms within each are less than width apart from each other along dir"""
        i=0
        sort_list=[]
        for atom in self.atoms:
            dist = dot(atom[1],plane)
            sort_list.append([i,dist])
            i+=1
        sort_list.sort(lambda x,y: cmp(x[1],y[1]))
        structs =[]
        curr_str = Structure()
        curr_str.basis = self.basis
        curr_dist = sort_list[0][1]
        for el in sort_list:
            if abs(el[1]-curr_dist) < width:
                curr_str.num_atoms +=1
                curr_str.atoms.append(self.atoms[el[0]])
            else:
                structs.append(curr_str)
                curr_dist = el[1]
                curr_str = Structure()
                curr_str.basis = self.basis
                curr_str.num_atoms +=1
                curr_str.atoms.append(self.atoms[el[0]])             
        return structs
        
        

if __name__ == "__main__":
    a=Structure()
    b=Structure()
    b.adatom("Li",0,0,0)
    c = Structure.add(b,a)
    print c.__str__()
    print a.__str__()