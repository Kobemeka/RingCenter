from rdkit import Chem
from rdkit.Chem import rdMolTransforms
import numpy as np
from matplotlib import cm
from copy import deepcopy
from fractions import Fraction
from datetime import datetime
import transform
'''
TODO!!: ntz ve bse dosyalarini duzelt
TODO: Point id kaldirilabilir
TODO: graphene, molekul ve best molecule icin farkli renklendirme -> FIXME: **kwargs
TODO: ayni anda hem rotation hem translation yapan fonksiyon ve bunun recursive olani
TODO: variable isimlerini daha anlasilir olanlarla degistir
TODO: molekule gore graphene olustur
'''


def transpose(list_):
    return list(zip(*list_))


def isInTol(value, target, tol):
    return target - tol <= value <= target + tol

def distance(p1, p2):
    return np.sqrt(np.power(p1.x - p2.x, 2) + np.power(p1.y - p2.y, 2))

def createFinal(poly,graphene,dx,dy,drot):

        rdMolTransforms.TransformConformer(poly.molecule.GetConformer(0), transform.transform(dx - poly.center_of_mass.coordinate.x,dy- poly.center_of_mass.coordinate.y,0,drot))
        
        rdMolTransforms.TransformConformer(graphene.molecule.GetConformer(0), transform.transform( - graphene.center_of_mass.coordinate.x, - graphene.center_of_mass.coordinate.y,0,0))
       
        poly_name = poly.file_path.split(".")[-2].split("/")[-1]
        newFileName = f"./final_mol_files/{poly_name}_final_{datetime.now().strftime('%d_%m_%Y-%H_%M_%S')}.sdf"
        writer = Chem.SDWriter(newFileName)
        mol = Chem.rdmolops.CombineMols(graphene.molecule,poly.molecule)
        writer.write(mol)
        writer.close()
        return newFileName

class Point:
    '''
        Creates a point in 3D
    '''

    def __init__(self, point_id, x: float, y: float, z: float, **kwargs):
        self.x = x
        self.y = y
        self.z = z
        self.point_id = point_id
        self.kwargs = kwargs

    def __repr__(self):
        return f"({self.x}, {self.y}, {self.z})"

    def __add__(self, other):
        return Point(self.point_id, self.x + other.x, self.y + other.y, self.z + other.z)

    def __radd__(self, other):
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
        return np.array([self.x, self.y, self.z])


class MolFile:
    '''
        Creates a molecule from a mol file and find rings.
        Takes file_path (path of the mol file) as parameter.
    '''

    def __init__(self, file_path: str, **kwargs):
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
                Point(atom_id, *self.atoms_positions[atom_id], **self.kwargs),
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
                    Point(
                        atom_id, *self.atoms_positions[atom_id], **self.kwargs),
                    self.atoms_mass[atom_id],
                    self.atoms_symbols[atom_id],
                    **self.kwargs
                )
                for atom_id in ring
            ] for ring in self.rdkit_rings
        ]

        self.rings = [
            Ring(ring_id, ring_atoms, **self.kwargs) for ring_id, ring_atoms in enumerate(self.atoms_of_rings)
        ]

        self.center_of_mass = DummyAtom(
            0,
            Point(0, *sum((np.array(self.atoms_positions).T *
                  self.atoms_mass).T)/sum(self.atoms_mass), **self.kwargs),
            **self.kwargs
        )

        self.prepare()

    def prepare(self):
        # we do not know how the molecule is positioned initially
        # to calculate the best we need to center the molecule to origin
        # TODO: consider implementing this function in Molecule class not MolFile
        # TODO: handle rotation

        for atom_index, atom in enumerate(self.atoms):
            atom.coordinate -= self.center_of_mass.coordinate

        for ring in self.rings:

            for ring_atom in ring.atoms:
                ring_atom.coordinate -= self.center_of_mass.coordinate
            ring.calculate()


class Molecule:
    def __init__(self, atoms, rings, center_of_mass) -> None:
        self.atoms = atoms
        self.rings = rings
        self.center_of_mass = center_of_mass

    def translation(self, axis, delta):
        ''' apply translation '''

        if axis == "x":
            translation_point = Point(0, delta, 0, 0)
        elif axis == "y":
            translation_point = Point(0, 0, delta, 0)
        elif axis == "z":
            translation_point = Point(0, 0, 0, delta)

        new_atoms = deepcopy(self.atoms)
        new_rings = deepcopy(self.rings)
        new_center_of_mass = deepcopy(self.center_of_mass)
        for atom in new_atoms:
            atom.coordinate += translation_point

        for ring in new_rings:

            for ring_atom in ring.atoms:
                ring_atom.coordinate += translation_point
            ring.calculate()
        
        new_center_of_mass.coordinate += translation_point
        return Molecule(new_atoms, new_rings, new_center_of_mass)

    def rotation(self, rotation_point: Point, rotation_axis, rotation_angle):
        ''' rotates the molecule in given axis by a rotation angle around a point.'''
        
        rotation_matrix = transform.xyzRotation(rotation_angle)[rotation_axis]

        new_atoms = deepcopy(self.atoms)
        new_rings = deepcopy(self.rings)
        new_center_of_mass = deepcopy(self.center_of_mass)

        for atom in new_atoms:
            atom_dist = atom.coordinate - rotation_point  # atom to rotation point distance

            atom.coordinate = Point(
                atom.coordinate.point_id,
                *np.matmul(rotation_matrix, atom_dist.to_list().T)
            )

        for ring in new_rings:

            for ring_atom in ring.atoms:
                ring_atom_dist = ring_atom.coordinate - \
                    rotation_point  # atom to rotation point distance
                ring_atom.coordinate = Point(
                    ring_atom.coordinate.point_id,
                    *np.matmul(rotation_matrix, ring_atom_dist.to_list().T)
                )

            ring.calculate()
        #####
        new_center_of_mass.coordinate = Point(
                self.center_of_mass.coordinate.point_id,
                *np.matmul(rotation_matrix, (self.center_of_mass.coordinate - rotation_point).to_list().T)
            )
        return Molecule(new_atoms, new_rings,new_center_of_mass)

    def translation2d(self, translation_range, delta):

        n = Fraction(str(float(translation_range))
                     ) // Fraction(str(float(delta))) + 1
        translated_molecules = []

        for ydelta in range(n):
            ydist = delta * ydelta
            ytranslation = self.translation("y", ydist)
            for xdelta in range(n):
                xdist = delta * xdelta
                xtranslation = ytranslation.translation("x", xdist)
                translated_molecules.append([(xdist, ydist), xtranslation])

        return translated_molecules

    def translation_rotation(self, translation_range, tdelta, rotation_range, rdelta):
        translation_count = Fraction(str(float(translation_range))) // Fraction(str(float(tdelta))) + 1

        rotation_count = Fraction(str(float(rotation_range))) // Fraction(str(float(rdelta))) + 1

        tr_molecules = []

        for ydelta in range(translation_count):
            ydist = tdelta * ydelta
            ytranslation = self.translation("y", ydist)

            for xdelta in range(translation_count):
                xdist = tdelta * xdelta
                xtranslation = ytranslation.translation("x", xdist)
                
                for r in range(rotation_count):
                    rotation_angle = rdelta * r
                    zrotation = xtranslation.rotation(self.center_of_mass.coordinate, "z", rotation_angle)
                    tr_molecules.append([(xdist, ydist, rotation_angle), zrotation])

        return tr_molecules

    def rotationZ(self, center, rrange, delta):
        # ?: center should be the center of mass of this molecule

        n = Fraction(str(float(rrange))) // Fraction(str(float(delta))) + 1
        rotated_molecules = []
        for r in range(n):
            rotation_angle = delta * r
            zrotation = self.rotation(
                center, "z", rotation_angle)  # rotate around origin
            rotated_molecules.append([rotation_angle, zrotation])
        return rotated_molecules

    def averageMinDistTranslationGraph(self, graphene, translation_range, delta, ax):

        n = Fraction(str(float(translation_range))
                     ) // Fraction(str(float(delta))) + 1

        avgdst = []

        vals = np.arange(0, translation_range, delta)
        xvals, yvals = np.meshgrid(vals, vals)

        for ydelta in range(n):
            ydist = delta * ydelta
            ytranslation = self.translation("y", ydist)
            yds = []
            for xdelta in range(n):
                xdist = delta * xdelta
                xtranslation = ytranslation.translation("x", xdist)
                yds.append(xtranslation.getAverageMinDist(graphene))
            avgdst.append(yds)
        ax.plot_surface(xvals, yvals, np.array(avgdst), cmap=cm.coolwarm)

    def getBestRotation(self, graphene, center, rrange, dr):
        return min(self.rotationZ(center, rrange, dr), key=lambda e: e[1].getAverageMinDist(graphene))

    def getBestTranslation(self, graphene, trange, dr):
        return min(self.translation2d(trange, dr), key=lambda e: e[1].getAverageMinDist(graphene))

    def getBestTranslationRotation(self,graphene,trange,rrange,dt,dr):
        return min(self.translation_rotation(trange, dt, rrange, dr), key=lambda e: e[1].getAverageMinDist(graphene))

    def getAverageMinDist(self, graphene):
        return np.array(self.centerAtomsDist(graphene)).mean()

    def recursiveBestTranslation(self, graphene, trange, dr, recursion_size):
        # TODO: reimplement it with new method or create a new function with new method

        temp_molecule = self.translation("x", 0)
        final_best = [0, 0]
        for r in range(recursion_size):
            best = temp_molecule.getBestTranslation(graphene, trange, dr)
            temp_molecule = self.translation(
                "x", best[0][0] - dr/2).translation("y", best[0][1] - dr/2)
            final_best[0] += best[0][0] - dr/2
            final_best[1] += best[0][1] - dr/2
            trange /= 10
            dr /= 10

        final_best[0] += best[0][0]
        final_best[1] += best[0][1]
        return final_best, self.translation("x", final_best[0]).translation("y", final_best[1])

    def recursiveBestRotation(self, graphene, center, rrange, delta, recursion_size):
        temp_molecule = self.rotation(center, "z", 0)
        final_best = 0
        for r in range(recursion_size):
            best = temp_molecule.getBestRotation(
                graphene, center, rrange, delta)
            temp_molecule = self.rotation(center, "z", best[0] - delta/2)
            final_best += best[0] - delta/2
            rrange /= 10
            delta /= 10

        final_best += best[0]
        return final_best, self.rotation(center, "z", final_best)

    def centerAtomsDist(self, graphene):
        ''' returns the miniumum distances of center of mass of rings to graphene atoms.'''

        cm_min_dists = []  # distances of center of mass to atoms
        for ring in self.rings:
            cm = ring.center_of_mass
            cm_dists = []
            for atom in graphene.atoms:
                cm_dists.append(distance(cm.coordinate, atom.coordinate))
            cm_min_dists.append(min(cm_dists))
        return cm_min_dists

    def draw(self, ax):
        '''
            Draw molecule rings.
        '''
        # ax.set_title(self.file_path)
        _ = [r.draw(ax) for r in self.rings]


class Atom:
    '''
        Creates an atom.
        Takes id, coordinate, mass and symbol as parameters.
    '''

    def __init__(self, atom_id: int, coordinate: Point, mass: float, symbol: str, **kwargs):
        self.atom_id = atom_id
        self.coordinate = coordinate
        self.mass = mass
        self.symbol = symbol
        self.kwargs = kwargs

    def __repr__(self):
        return f"Atom: {self.atom_id}, {self.symbol}"

    def draw(self, ax, atom_color="red"):
        ax.scatter(*self.coordinate.to_list(), c=atom_color)


class DummyAtom(Atom):
    def __init__(self, dummy_atom_id: int, coordinate: Point, **kwargs):
        super().__init__(dummy_atom_id, coordinate, 0, "Dummy Atom", **kwargs)

    def __repr__(self):
        return f"Dummy Atom: {self.atom_id}, {self.symbol}"


class Ring:
    def __init__(self, ring_id: int, atoms: list, **kwargs):
        self.ring_id = ring_id
        self.atoms = atoms
        self.kwargs = kwargs
        self.calculate()

    def calculate(self):

        self.atoms_coordinates = np.array(
            [atom.coordinate.to_list() for atom in self.atoms])
        self.atoms_mass = np.array([atom.mass for atom in self.atoms])
        self.total_mass = sum(self.atoms_mass)
        self.center_of_mass = DummyAtom(
            self.ring_id,
            Point(self.ring_id, *sum((self.atoms_coordinates.T *
                  self.atoms_mass).T)/self.total_mass, **self.kwargs),
            **self.kwargs
        )

        self.geometric_center = DummyAtom(
            self.ring_id,
            Point(self.ring_id, *sum(self.atoms_coordinates) /
                  len(self.atoms_coordinates), **self.kwargs),
            **self.kwargs
        )

    def draw(self, ax):
        if self.kwargs["show_ring_atoms"]:

            [a.draw(ax) for a in self.atoms]

        if self.kwargs["show_ring_center"]:

            self.geometric_center.draw(ax, atom_color="green")

        if self.kwargs["show_ring_center_of_mass"]:
            self.center_of_mass.draw(ax, atom_color="blue")

        if self.kwargs["show_ring"]:
            if "ring_color" not in self.kwargs:
                self.kwargs["ring_color"] = "k"
            rings_xyz = transpose(np.concatenate(
                (self.atoms_coordinates, [self.atoms_coordinates[0]])))
            ax.plot(rings_xyz[0], rings_xyz[1],
                    rings_xyz[2], self.kwargs["ring_color"])
