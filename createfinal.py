from molecule import *
import glob


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


gf = MolFile("./mol_files/graphene-big-square.mol",**graphene_settings)
graphene = Molecule(gf.atoms,gf.rings,gf.center_of_mass)
print(f"{graphene.center_of_mass.coordinate=}")
allMolFiles = [f for f in glob.glob("./test-files/*.mol")]
for molfile in [allMolFiles[0]]:
    print(molfile)
    # print("="*24)
    try:

        mf = MolFile(molfile,**show_settings)

        molecule = Molecule(mf.atoms,mf.rings,mf.center_of_mass)
        print(f"{molecule.center_of_mass.coordinate=}")

        test_range = 1.42
        test_rotation = np.pi/2
        test_ratio = 10

        gbtr = molecule.getBestTranslationRotation(graphene,test_range,test_rotation,test_range/test_ratio,test_rotation/test_ratio)
        print(f"{gbtr=}")
        finalmol = createFinal(mf,gf,*gbtr[0],"test-final")
        # finalmol = createFinal(mf,gf,5,5,np.pi/2,"final-test")
        print("="*10+"DONE"+"="*10)

    except Exception as e:
        print(f"{molfile} HATA!")