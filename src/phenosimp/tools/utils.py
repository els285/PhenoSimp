import uproot
import numpy as np
import awkward as ak
import vector
import h5py 
import sys

def awkward_to_dict(arr:ak.Array) -> dict:
    return {k:arr[k] for k in arr.fields}

def save_awkward_to_ROOT(arr: ak.Array , output_file: str , treename: str = "tree"):
    
    """
    Takes an ragged awkward array and saves it to a ROOT file
    """
    print(f"Writing cut array to {output_file}")
    d = {k:arr[k] for k in arr.fields}
    with  uproot.recreate(output_file) as file:
        file[treename] = d


def combine_awkward_arrays(list_of_arrays):
    
    """
    Takes a list of arrays with whatever fields.
    Returns one combined array with all thsoe fields
    """
    main_array = list_of_arrays[0]
    for array in list_of_arrays[1:]:
        for field in array.fields: 
            main_array[field] = array[field] 
            
    return main_array
    
    
def build_object_arrays(tree):
    
    """
    Parses the reco-tree and returns individual object trees
    """
    
    electron_array = tree.arrays(["el_pt","el_eta","el_phi","el_charge","electron_matched_indices"])
    electron_mass = 5.1110e-3
    electron_array["el_mass"] = electron_mass*ak.ones_like(tree["el_pt"].array())
    
    muon_array = tree.arrays(["mu_pt","mu_eta","mu_phi","mu_charge","muon_matched_indices"])
    muon_mass = 0.10566
    muon_array["mu_mass"] = muon_mass*ak.ones_like(tree["mu_pt"].array())
    
    jet_array = tree.arrays(["jet_pt","jet_eta","jet_phi","jet_mass","jet_btag","jet_tautag","jet_matched_indices"])
    
    met_array = tree.arrays(["met_met","met_phi","met_eta"])
    
    return {"electrons" : electron_array,
            "muons"     : muon_array,
            "jets"      :jet_array,
            "met"       :met_array}

required_branches = ['EventNumber', 
                     'el_pt', 'el_eta', 'el_phi', 'el_charge', 
                     'mu_pt', 'mu_eta', 'mu_phi', 'mu_charge', 
                     'jet_pt', 'jet_eta', 'jet_phi', 'jet_mass', 'jet_btag'
                     'met_met', 'met_eta', 'met_phi', 
                     'jet_matched_indices','electron_matched_indices', 'muon_matched_indices', 
                     'contains_duplicates', 'fully_matched']


def pad_variable(variable, max_len, pad_to = np.nan):
    padded_variable = ak.pad_none(variable, max_len, axis=1, clip=True)
    return ak.fill_none(padded_variable, pad_to)



        
