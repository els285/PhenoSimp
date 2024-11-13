import sys
import uproot
import awkward as ak
import vector

from truth_tools import parse_tops_and_Ws, parse_decays, parse_meta, parse_reco

def main(inname,outname):
    # Load
    print("Loading file")
    tree = uproot.open(f"{inname}:Delphes")

    print("Parsing truth")
    # Initialise truth dictionary
    d = {}

    print("Parsing event metadata")
    d = parse_meta(tree,d)
    
    print("Parsing top and W information")
    d = parse_tops_and_Ws(tree,d)
    assert ak.count_nonzero(ak.count(d["top_id"],axis=1) !=3 ) == 0, "There are events with other than four top quarks"
    assert ak.count_nonzero(ak.count(d["W_id"],axis=1)   !=3 ) == 0, "There are events with other than four on-shell W bosons"

    print("Parsing decay information")
    d = parse_decays(tree,d)
    assert ak.count_nonzero(ak.count(d["W_decay_id"],axis=1) !=6 ) == 0, "There are events with other than four bottom quarks"

    print("Truth particle parsing complete")

    ### Reco-level
    print("Writing reco-level trees")
    r = {}

    r = parse_reco(tree,r)

    print("Writing file")
    with  uproot.recreate(f"{outname}") as file:
        file["Truth"] = d
        file["Reco"]  = r
    print("Complete")
    
    
if __name__ == "__main__":
    inname =  sys.argv[1]
    outname = sys.argv[2]
    main(inname,outname)