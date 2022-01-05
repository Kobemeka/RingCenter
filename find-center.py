from rdkit import Chem
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import os
import argparse
# cc_length = 1.4289*np.sqrt(3)

textbox_alpha = 0.75
class Point:
    def __init__(self, x:float, y:float, z:float, id_ = None, **kwargs):
        self.x = x
        self.y = y
        self.z = z
        self.id = id_
        self.kwargs = kwargs
        self.xyz = np.array([x,y,z])

    def __repr__(self):
        return f"({self.x}, {self.y}, {self.z})"

class MolFile:
    def __init__(self,file_, **kwargs):
        self.file = file_
        self.file_name, self.file_extension = os.path.splitext(self.file)
        self.kwargs = kwargs
        self.mol = Chem.MolFromMolFile(self.file)
        self.atom_num = self.mol.GetNumAtoms()
        self.explen = self.kwargs["expected_length"]
        self.tolerance = self.kwargs["tolerance"]
        self.atoms = [
            Atom(
                    i,
                    Point(*list(self.mol.GetConformer().GetAtomPosition(i)),**self.kwargs),
                    **self.kwargs
                ) for i in range(self.atom_num)
            ]
        self.rings = [
            Ring(
                    [
                        Atom(j,Point(*list(self.mol.GetConformer().GetAtomPosition(j)),**self.kwargs),**self.kwargs) for j in i
                    ],
                    rid,
                    **self.kwargs
                ) for rid,i in enumerate(Chem.GetSymmSSSR(self.mol))
            ]
        self.ring_centers = [r.center for r in self.rings]
        self.ccombinations = list(combinations(self.ring_centers,2))

    def compute_center_distances(self):
        distances = []
        for t in range(len(self.ccombinations)):
            c1 = self.ccombinations[t][0]
            c2 = self.ccombinations[t][1]
            dist = distance(c1.coordinate.xyz,c2.coordinate.xyz)
            distances.append(f"{c1.id} {c2.id} {dist} {in_tolerance(dist,self.explen,self.tolerance)} {dist/self.explen}")
        return distances

    def draw(self,figure):
        figure.set_title(self.file)
        _ = [r.draw(figure) for r in self.rings]
        for t in range(len(self.ccombinations)):
            c1 = self.ccombinations[t][0]
            c2 = self.ccombinations[t][1]
            dist = distance(c1.coordinate.xyz,c2.coordinate.xyz)
            if in_tolerance(dist,self.explen,self.tolerance) and self.kwargs["show_center_distance"]:
                figure.plot([c1.coordinate.x,c2.coordinate.x],[c1.coordinate.y,c2.coordinate.y],color = "g")
                if self.kwargs["show_center_distance_text"]:
                    mftext = figure.text((c1.coordinate.x+c2.coordinate.x)/2,(c1.coordinate.y+c2.coordinate.y)/2,f"C{c1.id} - C{c2.id} -> {dist}")
                    mftext.set_bbox(dict(facecolor='green', alpha = textbox_alpha, edgecolor='green'))
    def dump(self):
        with open(f"{self.file_name}-center-info.log","w") as dump_file:
            dump_file.write(f"Expected length: {self.explen}\nTolerance: {self.tolerance}\n")
            dump_file.write(f"First Center id - Second Center id - Distance - As Expected - Ratio\n")
            dump_file.write("\n".join(self.compute_center_distances()))
        print(f"{self.file_name} center distance info is dumped")
class Atom:
    def __init__(self, id_:int,coordinates:list,**kwargs):
        self.id = id_
        self.coordinates = coordinates
        self.kwargs = kwargs

    def __repr__(self):
        return f"({self.id}, {self.coordinates})"

    def draw(self, figure):
        figure.scatter(self.coordinates.x,self.coordinates.y,c = "k")
        if self.kwargs["show_atom_number"]:
            atext = figure.text(self.coordinates.x,self.coordinates.y,self.id)
            atext.set_bbox(dict(facecolor='aqua', alpha = textbox_alpha, edgecolor='aqua'))

class DummyAtom:
    def __init__(self, coordinate, id_ = None, **kwargs):
        self.coordinate = coordinate
        self.id = id_
        self.kwargs = kwargs

    def __repr__(self):
        return f"({self.coordinate})"

    def draw(self,figure):

        figure.scatter(self.coordinate.x,self.coordinate.y, c = "r")

        if self.kwargs["show_center_number"]:

            ptext = figure.text(self.coordinate.x,self.coordinate.y,f"C-{self.id}")
            ptext.set_bbox(dict(facecolor='red', alpha = textbox_alpha, edgecolor='red'))

class Ring:

    def __init__(self,atoms:list,id_,**kwargs):
        self.atoms = atoms
        self.id = id_
        self.kwargs = kwargs
        self.atoms_coordinates = [a.coordinates.xyz for a in self.atoms]
        self.center = DummyAtom(Point(*sum(self.atoms_coordinates)/len(self.atoms_coordinates)),self.id,**self.kwargs)
    
    def __repr__(self):
        return f"({self.atoms})"


    def draw(self,figure):
        if self.kwargs["show_ring_atoms"]:

            _ = [a.draw(figure) for a in self.atoms]

        if self.kwargs["show_ring_center"]:

            self.center.draw(figure)

        if self.kwargs["show_ring"]:
            tac = transpose(self.atoms_coordinates + [self.atoms_coordinates[0]])
            figure.plot(tac[0],tac[1],"k")


def distance(p1,p2):
    # calculates the distance between two points: p1 and p2
    # p1 and p2 are tuples or list

    return np.sqrt(sum(list(map(lambda x: abs(x[0]-x[1])**2,list(zip(*[p1,p2]))))))

def transpose(list_):
    return list(zip(*list_))

def in_tolerance(value, actual_value, tolerance):
    return actual_value + tolerance >= value >= actual_value - tolerance


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--file", default = "graphene.mol")
    parser.add_argument("--center_number",default=False)
    parser.add_argument("--center_distance",default=True)
    parser.add_argument("--center_distance_text",default=False)
    parser.add_argument("--atom_number",default=False)
    parser.add_argument("--ring",default=True)
    parser.add_argument("--ring_atoms",default=False)
    parser.add_argument("--ring_center",default=True)
    parser.add_argument("--dump_center_info",default=False)
    parser.add_argument("--draw_file",default=True)
    parser.add_argument("--expected_length",default=1.42)
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
        "expected_length": args.expected_length,
        "tolerance": args.tolerance
    }

    mf = MolFile(args.file,**show_settings)
    if args.dump_center_info:
        mf.dump()
    if args.draw_file:
        fig, ax = plt.subplots()
        mf.draw(ax)
        plt.show()

