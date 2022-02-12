from rdkit import Chem
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from itertools import combinations
import os
import argparse

def transpose(list_):
    return list(zip(*list_))

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
    def create(self):
        pass
    def energy(self):
        pass

    def translation(self,axis,delta):
        if axis == "x":
            translation_point = Point(0,delta,0,0)
        
        elif axis == "y":
            translation_point = Point(0,0,delta,0)
        elif axis == "z":
            translation_point = Point(0,0,0,delta)

        for atom_index, atom in enumerate(self.atoms):
            atom.coordinate += translation_point

        # for ring_index, ring in enumerate(self.rings):
        #     for ring_atom in enumerate(ring.atoms):
        #         ring_atom.coordinate += translation_point
            
    def draw(self,ax):
        '''
            Draw molecule rings.
        '''
        ax.set_title(self.file_path)
        _ = [r.draw(ax) for r in self.rings]
        # for t in range(len(self.ccombinations)):
        #     c1 = self.ccombinations[t][0]
        #     c2 = self.ccombinations[t][1]
        #     dist = distance(c1.coordinate.xyz,c2.coordinate.xyz)

        #     if in_tolerance(dist,self.explen,self.tolerance) and self.kwargs["show_center_distance"]:
        #         ax.plot([c1.coordinate.x,c2.coordinate.x],[c1.coordinate.y,c2.coordinate.y],color = "g")
        #         if self.kwargs["show_center_distance_text"]:
        #             mftext = ax.text((c1.coordinate.x+c2.coordinate.x)/2,(c1.coordinate.y+c2.coordinate.y)/2,f"C{c1.id} - C{c2.id} -> {dist}")
        #             mftext.set_bbox(dict(facecolor='green', alpha = textbox_alpha, edgecolor='green'))
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

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--file", default = "./mol_files/polypyrrole.mol")
    parser.add_argument("--center_number",default=False)
    parser.add_argument("--center_distance",default=False)
    parser.add_argument("--center_distance_text",default=False)
    parser.add_argument("--atom_number",default=False)
    parser.add_argument("--ring",default=True)
    parser.add_argument("--ring_atoms",default=False)
    parser.add_argument("--ring_center",default=False)
    parser.add_argument("--ring_center_of_mass",default=True)
    parser.add_argument("--dump_center_info",default=False)
    parser.add_argument("--draw_file",default=True)
    parser.add_argument("--expected_distance",default=1.42)
    parser.add_argument("--tolerance",default=1e-2)


    args = parser.parse_args()

    show_settings = {
        "show_center_number":args.center_number,
        "show_center_distance":args.center_distance,
        "show_center_distance_text":args.center_distance_text,
        "show_atom_number":args.atom_number,
        "show_ring":args.ring,
        "show_ring_atoms":args.ring_atoms,
        "show_ring_center":args.ring_center,
        "show_ring_center_of_mass":args.ring_center_of_mass,
        "expected_distance": args.expected_distance,
        "tolerance": args.tolerance
    }

    mf = Molecule(args.file,**show_settings)

    if args.dump_center_info:
        mf.dump()

    if args.draw_file:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # fig, ax = plt.subplots()
        ax.set_xlabel("X - AXIS")
        ax.set_ylabel("Y - AXIS")
        ax.set_zlabel("Z - AXIS")
        mf.draw(ax)
        mf.translation("x",1)
        mf.draw(ax)
        plt.show()

