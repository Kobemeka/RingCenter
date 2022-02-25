from molecule import *
import argparse
import time
if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--file", default = "./mol_files/polypyrrole.mol")
    # parser.add_argument("--file", default = "./mol_files/pedot24_new.mol")
    # parser.add_argument("--file", default = "./mol_files/ppp24_new.mol")
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

    mf = MolFile(args.file,**show_settings)
    gf = MolFile("./mol_files/graphene_long.mol",**show_settings)
    
    molecule = Molecule(mf.atoms,mf.rings)
    graphene = Molecule(gf.atoms,gf.rings)

    # translated_molecule = molecule.translation("x",10)

    tr2d = molecule.translation2d(2,1)


    if args.dump_center_info:
        mf.dump()

    if args.draw_file:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        ax.set_xlabel("X - AXIS")
        ax.set_ylabel("Y - AXIS")
        ax.set_zlabel("Z - AXIS")

        test_range = 2
        test_d = 0.5

        molecule.draw(ax)
        # graphene.draw(ax)
        # translated_molecule.draw(ax)
        # [i[1].draw(ax) for i in tr2d]

        # print(molecule.centerAtomsDist(graphene))

        # besttr = molecule.getBestTranslation(graphene,test_range,test_d)
        # ax.scatter(besttr[0][0],besttr[0][1])

        # # print(besttr[0])
        # # # besttr[1].draw(ax)
        
        # rbst = molecule.recursiveBestTranslation(graphene,test_range,test_d,4)
        # print(rbst)
        # ax.scatter([rbst[1][0] + rbst[2][0]],[rbst[1][1] + rbst[2][1]])

        # molecule.averageMinDistGraph(graphene,test_range,test_d,ax)

        # rotm = molecule.rotationZ(Point(0,0,0,0),np.pi/2,np.pi/10)
        # [i[1].translation("z",2).draw(ax) for i in rotm]
        gbr = molecule.getBestRotation(graphene,Point(0,0,0,0),np.pi/2,np.pi/20)
        print(gbr[0])
        gbr[1].translation("z",2).draw(ax)
        plt.show()