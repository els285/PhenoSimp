import sys
import uproot
import awkward as ak
import vector

from truth_tools import parse_tops_and_Ws, parse_decays, parse_meta, parse_reco, do_matching, examine_matched_indices

def main(inname,outname):
    # Load
    print("Loading file")
    tree = uproot.open(f"{inname}:Delphes")

    print("Parsing truth")
    # Initialise truth dictionary
    truth_dict = {}

    print("Parsing event metadata")
    truth_dict = parse_meta(tree,truth_dict)
    
    print("Parsing top and W information")
    truth_dict = parse_tops_and_Ws(tree,truth_dict)
    assert ak.count_nonzero(ak.count(truth_dict["top_id"],axis=1) !=4 ) == 0, "There are events with other than four top quarks"
    assert ak.count_nonzero(ak.count(truth_dict["W_id"],axis=1)   !=4 ) == 0, "There are events with other than four on-shell W bosons"

    print("Parsing decay information")
    truth_dict = parse_decays(tree,truth_dict)
    assert ak.count_nonzero(ak.count(truth_dict["W_decay_id"],axis=1) !=8 ) == 0, "There are events with other than four bottom quarks"

    print("Truth particle parsing complete")

    ### Reco-level
    print("Writing reco-level trees")
    reco_dict = {}

    # Copy basic reco branches
    reco_dict = parse_reco(tree,reco_dict)
    
    # Find the matched indices 
    print("Performing jet and lepton matching")
    reco_dict = do_matching(tree,truth_dict,reco_dict)
    
    list_of_matched_indices = [reco_dict["jet_matched_indices"],
                             reco_dict["electron_matched_indices"],
                             reco_dict["muon_matched_indices"]]
    
    # Examine matched indices
    reco_dict = examine_matched_indices(list_of_matched_indices,reco_dict,12)
    
    print("Writing file")
    with  uproot.recreate(f"{outname}") as file:
        file["Truth"] = truth_dict
        file["Reco"]  = reco_dict
    print("Complete")
    
if __name__ == "__main__":
    inname =  sys.argv[1]
    outname = sys.argv[2]
    main(inname,outname)