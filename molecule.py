from rdkit import Chem
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse

def transpose(list_):
    return list(zip(*list_))

def isInTol(value,target,tol):
    return target - tol <= value <= target + tol

def distance(p1,p2):
    return np.sqrt(np.power(p1.x - p2.x,2) + np.power(p1.y - p2.y,2))

class Point:
    '''
        Creates a point in 3D
    '''
    def __init__(self,point_id, x:float, y:float, z:float, **kwargs):
        self.x = x
        self.y = y
        self.z = z
        self.point_id = point_id
        self.kwargs = kwargs

    def __repr__(self):
        return f"({self.x}, {self.y}, {self.z})"

    def __add__(self, other):
        return Point(self.point_id, self.x + other.x, self.y + other.y, self.z + other.z)

    def __radd__(self,other):
        if other == 0:
            return self
        else:
            return self.__add__(other)
    
    def __sub__(self, other):
        return Point(self.point_id, self.x - other.x, self.y - other.y, self.z - other.z)
    
    def to_list(self):
        '''
            Returns a 1D numpy array: [x,y,z]
        '''
        return np.array([self.x,self.y,self.z])

class Molecule:
    '''
        Creates a molecule from a mol file and find rings.
        Takes file_path (path of the mol file) as parameter.
    '''
    def __init__(self,file_path:str,**kwargs):
        self.file_path = file_path
        self.kwargs = kwargs
        self.expected_distance = self.kwargs["expected_distance"]
        self.tolerance = self.kwargs["tolerance"]

        self.molecule = Chem.MolFromMolFile(self.file_path)

        self.rdkit_atoms = list(self.molecule.GetAtoms())
        self.atoms_mass = [a.GetMass() for a in self.rdkit_atoms]
        self.atoms_symbols = [a.GetSymbol() for a in self.rdkit_atoms]

        self.total_atoms = self.molecule.GetNumAtoms()

        self.molecule_conformer = self.molecule.GetConformer()
        self.atoms_positions = self.molecule_conformer.GetPositions()
        self.rdkit_rings = list(Chem.GetSymmSSSR(self.molecule))
        
        self.create()
    
    def create(self):

        self.atoms = [
            Atom(
                atom_id,
                Point(atom_id,*self.atoms_positions[atom_id],**self.kwargs),
                self.atoms_mass[atom_id],
                self.atoms_symbols[atom_id],
                **self.kwargs
            )
            for atom_id in range(self.total_atoms)
        ]

        self.atoms_of_rings = [
            [
                Atom(
                    atom_id,
                    Point(atom_id,*self.atoms_positions[atom_id],**self.kwargs),
                    self.atoms_mass[atom_id],
                    self.atoms_symbols[atom_id],
                    **self.kwargs
                )
                for atom_id in ring
            ] for ring in self.rdkit_rings
        ]

        self.rings = [
            Ring(ring_id,ring_atoms,**self.kwargs) for ring_id,ring_atoms in enumerate(self.atoms_of_rings)
        ]

        self.center_of_mass = DummyAtom(
            0,
            Point(0,*sum((np.array(self.atoms_positions).T * self.atoms_mass).T)/sum(self.atoms_mass),**self.kwargs),
            **self.kwargs
        )

        self.prepare()
    def prepare(self):
        # we do not know how the molecule is positioned initially
        # to calculate the best we need to center the molecule to origin
        # TODO: handle rotation

        for atom_index, atom in enumerate(self.atoms):
            atom.coordinate -= self.center_of_mass.coordinate

        for ring in self.rings:
            
            for ring_atom in ring.atoms:
                ring_atom.coordinate -= self.center_of_mass.coordinate
            ring.calculate()

    def energy(self):
        pass

    def translation(self,axis,delta):
        ''' apply translation '''

        if axis == "x":
            translation_point = Point(0,delta,0,0)
        
        elif axis == "y":
            translation_point = Point(0,0,delta,0)
        elif axis == "z":
            translation_point = Point(0,0,0,delta)

        for atom_index, atom in enumerate(self.atoms):
            atom.coordinate += translation_point

        for ring in self.rings:
            
            for ring_atom in ring.atoms:
                ring_atom.coordinate += translation_point
            ring.calculate()
        return self

    def rotate(self,rotation_point: Point, rotation_axis, rotation_angle):
        ''' rotates the molecule in given axis by a rotation angle around a point.'''
        if rotation_axis == "x":
            rotation_matrix = np.array([
                [1, 0, 0],
                [0, np.cos(rotation_angle), -np.sin(rotation_angle)],
                [0, np.sin(rotation_angle), np.cos(rotation_angle)]
            ])
        elif rotation_axis == "y":
            rotation_matrix = np.array([
                [np.cos(rotation_angle), 0, np.sin(rotation_angle)],
                [0, 1, 0],
                [-np.sin(rotation_angle), 0, np.cos(rotation_angle)]
            ])
        elif rotation_axis == "z":
            rotation_matrix = np.array([
                [np.cos(rotation_angle), -np.sin(rotation_angle),0],
                [np.sin(rotation_angle), np.cos(rotation_angle),0],
                [0, 0, 1]
            ])
            
        for atom_index, atom in self.atoms:
            atom_dist = atom.coordinate - rotation_point # atom to rotation point distance

            atom.coordinate = Point(
                    atom.coordinate.point_id,
                    *np.matmul(rotation_matrix, atom_dist.to_list().T)
                    )

        for ring in self.rings:
            
            for ring_atom in ring.atoms:
                ring_atom_dist = ring_atom.coordinate - rotation_point # atom to rotation point distance
                ring_atom.coordinate = Point(
                    ring_atom.coordinate.point_id,
                    *np.matmul(rotation_matrix, ring_atom_dist.to_list().T)
                    )

            ring.calculate()
        return self
    
    def bestTranslation(self,graphene,trange,dr):
        ''' apply translations to the molecule in a given range and the step in both axis (x, y)
            then returns the distances with the number of translations: [(y,x), distance]
        '''
    
        n = int(trange//dr)
        cmdists = []
        for ydelta in range(n):
            self.translation("y",dr)
            for xdelta in range(n):
                self.translation("x",dr)
                cmd = [(ydelta,xdelta),self.centerAtomsDist(graphene)]
                # print(cmd)
                cmdists.append(cmd)

        return min(cmdists,key=lambda x: x[1]) # TODO: maybe change the min function.

    def getBestTranslation(self,graphene,trange,dr):
        self.create() # reset all translations and rotations
        best = self.bestTranslation(graphene,trange,dr) # find best translation
        self.create() # reset again
        # get the best number of translation in y axis (best[0][0]) and multiply with the step (dr) to translate to best position
        self.translation("y",best[0][0] * dr)
        # get the best number of translation in x axis (best[0][1]) and multiply with the step (dr) to translate to best position
        self.translation("x",best[0][1] * dr)
        return best
        
    def recursiveBestTranslation(self,graphene,trange,dr,recursion_size):
        # TODO: optimize and re think it

        self.create()
        recrange = trange
        delta = dr
        bests = []
        for r in range(recursion_size):
            best = self.bestTranslation(graphene,recrange,delta)
            # self.create()
            self.translation("y",(best[0][0] - 1) * delta)
            self.translation("x",(best[0][1] - 1) * delta)
            bests.append([best[0],recrange,delta])
            recrange /= 5
            delta /= 10
        return bests

    # TODO: implement for rotation
    # def rotateInRange(molecule,rotation_point,axis,rrange,dt):
        
    #     n = drange//dt
    #     molecule_states = []

    #     for delta in range(n):

    #         molecule = molecule.rotate(rotation_point,axis,dt)

    #         molecule_states.append(molecule)
    #     return molecule_states
    def centerAtomsDist(self,graphene):
        ''' returns the miniumum distances of center of mass of rings to graphene atoms.'''

        cm_min_dists = [] # distances of center of mass to atoms
        for ring in self.rings:
            cm = ring.center_of_mass
            cm_dists = []
            for atom in graphene.atoms:
                cm_dists.append(distance(cm.coordinate,atom.coordinate))
            cm_min_dists.append(min(cm_dists))
        return cm_min_dists

    def draw(self,ax):
        '''
            Draw molecule rings.
        '''
        ax.set_title(self.file_path)
        _ = [r.draw(ax) for r in self.rings]
class Atom:
    '''
        Creates an atom.
        Takes id, coordinate, mass and symbol as parameters.
    '''
    def __init__(self,atom_id:int,coordinate:Point,mass:float,symbol:str,**kwargs):
        self.atom_id = atom_id
        self.coordinate = coordinate
        self.mass = mass
        self.symbol = symbol
        self.kwargs = kwargs
    def __repr__(self):
        return f"Atom: {self.atom_id}, {self.symbol}"
    def draw(self,ax,atom_color = "red"):
        ax.scatter(*self.coordinate.to_list(),c=atom_color)

class DummyAtom(Atom):
    def __init__(self,dummy_atom_id:int,coordinate:Point,**kwargs):
        super().__init__(dummy_atom_id,coordinate,0,"Dummy Atom",**kwargs)

    def __repr__(self):
        return f"Dummy Atom: {self.atom_id}, {self.symbol}"

class Ring:
    def __init__(self,ring_id:int,atoms:list,**kwargs):
        self.ring_id = ring_id
        self.atoms = atoms
        self.kwargs = kwargs
        self.calculate()

    def calculate(self):

        self.atoms_coordinates = np.array([atom.coordinate.to_list() for atom in self.atoms])
        self.atoms_mass = np.array([atom.mass for atom in self.atoms])
        self.total_mass = sum(self.atoms_mass)
        self.center_of_mass = DummyAtom(
            self.ring_id,
            Point(self.ring_id,*sum((self.atoms_coordinates.T * self.atoms_mass).T)/self.total_mass,**self.kwargs),
            **self.kwargs
        )

        self.geometric_center = DummyAtom(
            self.ring_id,
            Point(self.ring_id,*sum(self.atoms_coordinates)/len(self.atoms_coordinates),**self.kwargs),
            **self.kwargs
        )

    def draw(self,ax):
        if self.kwargs["show_ring_atoms"]:

            [a.draw(ax) for a in self.atoms]

        if self.kwargs["show_ring_center"]:

            self.geometric_center.draw(ax,atom_color="green")

        if self.kwargs["show_ring_center_of_mass"]:
            self.center_of_mass.draw(ax,atom_color="blue")

        if self.kwargs["show_ring"]:

            rings_xyz = transpose(np.concatenate((self.atoms_coordinates,[self.atoms_coordinates[0]])))
            ax.plot(rings_xyz[0],rings_xyz[1],rings_xyz[2],"k")