import sys
from LHEclass import LHEparse

def main(input_lhe,output_root,outtree_name):
    
    P = LHEparse(input_lhe)
    P.build()
    P.write_kinematics_to_ROOT(output_root,outtree_name)    
    
if __name__ == "__main__":
    input_lhe    = sys.argv[1]
    output_root  = sys.argv[2]
    outtree_name = sys.argv[3] if sys.argv[3] else "tree"
    main(input_lhe,output_root,outtree_name)