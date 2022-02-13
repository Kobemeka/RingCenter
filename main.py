from molecule import *
if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    # parser.add_argument("--file", default = "./mol_files/polypyrrole.mol")
    parser.add_argument("--file", default = "./mol_files/ppp24_new.mol")
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

    graphene = Molecule("./mol_files/graphene_long.mol",**show_settings)

    if args.dump_center_info:
        mf.dump()

    if args.draw_file:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        ax.set_xlabel("X - AXIS")
        ax.set_ylabel("Y - AXIS")
        ax.set_zlabel("Z - AXIS")

        # graphene.create()
        graphene.draw(ax)

        mf.draw(ax)
        # mf.create()
        # mf.draw(ax)

        # cmds = mf.translateCheck(graphene,1,0.1)
        # print(cmds)
        # gbt = mf.getBestTranslation(graphene,2,0.1)
        # print("gbt: ", gbt)
        rbt = mf.recursiveBestTranslation(graphene,2,0.1,2)
        print(rbt)
        mf.draw(ax)
        plt.show()