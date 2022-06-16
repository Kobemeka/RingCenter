from molecule import *
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from writeMolFile import MolFileWriter

show_settings = {
    "show_center_number":False,
    "show_center_distance":False,
    "show_center_distance_text":False,
    "show_atom_number":False,
    "show_ring":True,
    "show_ring_atoms":False,
    "show_ring_center":False,
    "show_ring_center_of_mass":True,
    "expected_distance": 1.42,
    "tolerance": 1e-2
}
graphene_settings = {
    "show_center_number":False,
    "show_center_distance":False,
    "show_center_distance_text":False,
    "show_atom_number":False,
    "show_ring":True,
    "show_ring_atoms":True,
    "show_ring_center":False,
    "show_ring_center_of_mass":False,
    "expected_distance": 1.42,
    "tolerance": 1e-2
}


gf = MolFile("mol_files/graphene-piece.mol",**graphene_settings)
graphene = Molecule(gf.atoms,gf.rings,gf.center_of_mass,gf.bonds_ids,True)
allMolFiles = [f for f in glob.glob("./correct/*.mol")]

test_range = 2
test_rotation = np.pi*2
test_ratio = 15
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

molfile = "./correct/pth24.mol"
molfileid = 0
print("="*24)
print(f"{molfileid}/{len(allMolFiles)}")
print(molfile)

mf = MolFile(molfile,**show_settings)

molecule = Molecule(mf.atoms,mf.rings,mf.center_of_mass,mf.bonds_ids,True)

gbtr = molecule.getBestTranslationRotation(graphene,test_range,test_rotation,test_range/test_ratio,test_rotation/test_ratio)

writer = MolFileWriter([gbtr[1],graphene])
writer.write(f"./final/{molfile.split('.')[-2].split('/')[-1]}-final-{test_ratio}.mol")
print(" "*10+"DONE"+" "*10+"\n")
print("="*24)
# for molfileid, molfile in enumerate(allMolFiles):
#     print("="*24)
#     print(f"{molfileid}/{len(allMolFiles)}")
#     print(molfile)

#     mf = MolFile(molfile,**show_settings)

#     molecule = Molecule(mf.atoms,mf.rings,mf.center_of_mass,mf.bonds_ids,True)

#     gbtr = molecule.getBestTranslationRotation(graphene,test_range,test_rotation,test_range/test_ratio,test_rotation/test_ratio)

#     writer = MolFileWriter([gbtr[1],graphene])
#     writer.write(f"./final/{molfile.split('.')[-2].split('/')[-1]}-final-{test_ratio}.mol")
#     print(" "*10+"DONE"+" "*10+"\n")
#     print("="*24)