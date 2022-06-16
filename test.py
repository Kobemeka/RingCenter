from molecule import *
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--file", default = "./correct/fdtbt-A.mol")
    # parser.add_argument("--file", default = "./mol_fi-ales/ntz.mol")
    # parser.add_argument("--file", default = "mol_files/bo.mol")
    # parser.add_argument("--file", default = "mol_files/a33.mol")
    # parser.add_argument("--file", default = "mol_files/bse.mol")
    
    parser.add_argument("--center_number",default=False)
    parser.add_argument("--center_distance",default=False)
    parser.add_argument("--center_distance_text",default=False)
    parser.add_argument("--atom_number",default=False)
    parser.add_argument("--ring",default=True)
    parser.add_argument("--ring_atoms",default=False)
    parser.add_argument("--ring_center",default=False)
    parser.add_argument("--ring_center_of_mass",default=True)
    parser.add_argument("--dump_center_info",default=False)
    parser.add_argument("--draw_file",default=False)
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
        "tolerance": args.tolerance,
        "isGraphene": False

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
        "expected_distance": args.expected_distance,
        "tolerance": args.tolerance,
        "isGraphene": True
    }

    mf = MolFile(args.file,**show_settings)
    gf = MolFile("./mol_files/graphene-big-square.mol",**graphene_settings)
    molecule = Molecule(mf.atoms,mf.rings,mf.center_of_mass,mf.bonds_ids)
    graphene = Molecule(gf.atoms,gf.rings,gf.center_of_mass,gf.bonds_ids)


    if args.dump_center_info:
        mf.dump()

    if args.draw_file:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        ax.set_xlabel("X - AXIS")
        ax.set_ylabel("Y - AXIS")
        ax.set_zlabel("Z - AXIS")

        test_range = 1.42
        test_rotation = np.pi/2
        test_ratio = 5

        
        # gbtr = molecule.getBestTranslationRotation(graphene,test_range,test_rotation,test_range/test_ratio,test_rotation/test_ratio)
        gbtr = molecule.recursiveBestTransform(graphene,molecule.center_of_mass.coordinate,test_range,test_rotation,test_range/test_ratio,test_rotation/test_ratio,4)
        # gbtr = molecule.recursiveBestTranslation(graphene,test_rotation,test_range/test_ratio,4)
        print(gbtr)
        # print(gbtr[0])
        
        # graphene.draw(ax)
        gbtr[1].draw(ax)

        '''
            STSRT FINAL SDF FILE
        '''
        # finalmol = createFinal(mf,gf,*gbtr[0],ax,"final-test-15")
        # fmol = MolFile(finalmol,**show_settings)
        # fm = Molecule(fmol.atoms,fmol.rings,fmol.center_of_mass)
        # fm.draw(ax)

        '''
            END FINAL SDF FILE
        '''        

        plt.show()