import uproot
import awkward as ak
import vector
import numpy as np 

def parse_tops_and_Ws(tree,d):
    _22mask = tree["Particle.Status"].array()==22
    top_mask = abs(tree["Particle.PID"].array())==6
    W_mask   = abs(tree["Particle.PID"].array())==24

    d["top_pt"]     = tree["Particle.PT"].array()[_22mask & top_mask]
    d["top_eta"]    = tree["Particle.Eta"].array()[_22mask & top_mask]
    d["top_phi"]    = tree["Particle.Phi"].array()[_22mask & top_mask]
    d["top_e"]      = tree["Particle.E"].array()[_22mask & top_mask]
    d["top_mass"]   = tree["Particle.Mass"].array()[_22mask & top_mask]
    d["top_id"]     = tree["Particle.PID"].array()[_22mask & top_mask] 

    d["W_pt"]       = tree["Particle.PT"].array()[_22mask & W_mask]
    d["W_eta"]      = tree["Particle.Eta"].array()[_22mask & W_mask]
    d["W_phi"]      = tree["Particle.Phi"].array()[_22mask & W_mask]
    d["W_e"]        = tree["Particle.E"].array()[_22mask & W_mask]
    d["W_mass"]     = tree["Particle.Mass"].array()[_22mask & W_mask]
    d["W_id"]       = tree["Particle.PID"].array()[_22mask & W_mask] 
    
    return d


def parse_decays(tree,d):
        
    # Final-state particles
    status23    = tree["Particle.Status"].array()==23
    # Particles with a W mother
    has_a_W_mother = abs(tree["Particle.PID"].array()[tree["Particle.M1"].array()])==24
    # Particles which are bs
    bottom_mask = abs(tree["Particle.PID"].array())==5
    # Is a quark or lepton
    is_fermion = abs(tree["Particle.PID"].array())<17

    # B quarks
    print("Extracting b-quark information")

    d["b_id"]       = tree["Particle.PID"].array()[status23     & ~has_a_W_mother & bottom_mask]
    d["b_pt"]       = tree["Particle.PT"].array()[status23      & ~has_a_W_mother & bottom_mask]
    d["b_eta"]      = tree["Particle.Eta"].array()[status23     & ~has_a_W_mother & bottom_mask]
    d["b_phi"]      = tree["Particle.Phi"].array()[status23     & ~has_a_W_mother & bottom_mask]
    d["b_e"]        = tree["Particle.E"].array()[status23       & ~has_a_W_mother & bottom_mask]
    d["b_mass"]     = tree["Particle.Mass"].array()[status23    & ~has_a_W_mother & bottom_mask]

    # Wdecays
    print("Extracting W decay information")

    d["W_decay_id"]     = tree["Particle.PID"].array()[has_a_W_mother   & is_fermion]
    d["W_decay_pt"]     = tree["Particle.PT"].array()[has_a_W_mother    & is_fermion]
    d["W_decay_eta"]    = tree["Particle.Eta"].array()[has_a_W_mother   & is_fermion]
    d["W_decay_phi"]    = tree["Particle.Phi"].array()[has_a_W_mother   & is_fermion]
    d["W_decay_e"]      = tree["Particle.E"].array()[has_a_W_mother     & is_fermion]
    d["W_decay_mass"]   = tree["Particle.Mass"].array()[has_a_W_mother  & is_fermion]
    
    d["truth_decay_id"]   = ak.concatenate([d["b_id"][:,-4:]   , d["W_decay_id"]],   axis=1)
    d["truth_decay_pt"]   = ak.concatenate([d["b_pt"][:,-4:]   , d["W_decay_pt"]],   axis=1)
    d["truth_decay_eta"]  = ak.concatenate([d["b_eta"][:,-4:]  , d["W_decay_eta"]],  axis=1)
    d["truth_decay_phi"]  = ak.concatenate([d["b_phi"][:,-4:]  , d["W_decay_phi"]],  axis=1)
    d["truth_decay_e"]    = ak.concatenate([d["b_e"][:,-4:]    , d["W_decay_e"]],    axis=1)
    d["truth_decay_mass"] = ak.concatenate([d["b_mass"][:,-4:] , d["W_decay_mass"]], axis=1)

    return d 


def parse_meta(tree,d):
    
    d["EventNumber"]      = tree["Event.Number"].array()
    
    return d

def parse_reco(tree,r):

    r["EventNumber"]      = tree["Event.Number"].array()

    r["jet_pt"]     =  tree["Jet.PT"].array()
    r["jet_eta"]    =  tree["Jet.Eta"].array()
    r["jet_phi"]    =  tree["Jet.Phi"].array()
    r["jet_mass"]   =  tree["Jet.Mass"].array()
    r["jet_btag"]   = tree["Jet.BTag"].array()

    r["el_pt"]      =  tree["Electron.PT"].array()
    r["el_eta"]     =  tree["Electron.Eta"].array()
    r["el_phi"]     =  tree["Electron.Phi"].array()
    r["el_charge"]  =  tree["Electron.Charge"].array()

    r["mu_pt"]      =  tree["Muon.PT"].array()
    r["mu_eta"]      =  tree["Muon.Eta"].array()
    r["mu_phi"]      =  tree["Muon.Phi"].array()
    r["mu_charge"]  =  tree["Muon.Charge"].array()

    r["met_met"]    =  tree["MissingET.MET"].array()
    r["met_eta"]    =  tree["MissingET.Eta"].array()
    r["met_phi"]    =  tree["MissingET.Phi"].array()
    
    return r


def DR_matching(truth_vectors,reco_vectors,DR):
    
    # Perform matching
    delta_r_values = truth_vectors[:, None].deltaR(reco_vectors)

    # Find the minimum value in each row
    min_values = ak.min(delta_r_values, axis=2)

    # Find the indices of the minimum values in each row
    min_indices = ak.argmin(delta_r_values, axis=2)

    # Apply the condition: if no value is less than 0.4, set the result to -9
    matched_indices = ak.where(min_values < DR, min_indices, -9)
    
    return matched_indices

def do_matching(tree,truth_dict,reco_dict):
    
    """
    Performs truth-matching based on DR<some-min separately for jets, electrons
    and muons
    """

    # Building hadronic truth objects
    truth_objects = vector.zip({"pt":truth_dict["truth_decay_pt"],
                            "eta":truth_dict["truth_decay_eta"],
                            "phi":truth_dict["truth_decay_phi"],
                            "e":truth_dict["truth_decay_e"]})
    
    # truth_Bs = vector.zip({"pt":truth_dict["b_pt"],
    #                         "eta":truth_dict["b_eta"],
    #                         "phi":truth_dict["b_phi"],
    #                         "e":truth_dict["b_e"]})[:,-4:]

    # truth_W_decays = vector.zip({"pt":truth_dict["W_decay_pt"],
    #                             "eta":truth_dict["W_decay_eta"],
    #                             "phi":truth_dict["W_decay_phi"],
    #                             "e":truth_dict["W_decay_e"]})

    # truth_objects = ak.concatenate([truth_Bs,truth_W_decays],axis=1)

    # Building jets
    all_jets = vector.zip({"pt":tree["Jet.PT"].array(),
            "eta":tree["Jet.Eta"].array(),
            "phi":tree["Jet.Phi"].array(),
            "m":tree["Jet.Mass"].array()})
    
    # Build charged leptons
    electrons = vector.zip({"pt":tree["Electron.PT"].array(),
           "eta":tree["Electron.Eta"].array(),
           "phi":tree["Electron.Phi"].array(),
           "m":0.511e-3})
    
    muons = vector.zip({"pt":tree["Muon.PT"].array(),
           "eta":tree["Muon.Eta"].array(),
           "phi":tree["Muon.Phi"].array(),
           "m":0.105}) 
    
    # reco_charged_leptons = ak.concatenate([electrons,muons],axis=1)    
    
    # The fill_none is purely for the awkward data type to be int64 not ?int64
    # which is not writeable to ROOT 
    reco_dict["jet_matched_indices"]      = ak.fill_none(DR_matching(truth_objects,all_jets,0.4),-99)
    reco_dict["electron_matched_indices"] = ak.fill_none(DR_matching(truth_objects,electrons,0.1),-99)
    reco_dict["muon_matched_indices"]     = ak.fill_none(DR_matching(truth_objects,muons,0.1),-99)

    return reco_dict


def compute_if_fully_matched(matched_indices,max_index):
    
    """
    Computes whether an event is fully-matched based on each index only
    appearing once (but with the possibility for additional jets)
    Args:
    - matched_indices = the combined matched_indices
    - max_index = how many truth final states there should be e.g. for 4tops, 12
    (4b + 8 decay)
    Returns:
        An ak.array with entry 1 if event fully matched else 0
    """
    
    max_event_size = ak.max(ak.count(matched_indices,axis=1))
    padded_indices = ak.fill_none(ak.pad_none(matched_indices,target=max_event_size),np.nan).to_numpy()
    
    index_appears_mask = [] 
    for i in range(max_index):
        index_appears_mask.append(np.count_nonzero(np.isin(padded_indices, [i]), axis=1) == 1)
        
    return  ak.Array(1*np.logical_and.reduce(index_appears_mask, axis=0))

def compute_if_duplicate(matched_indices,max_index):
    
    """
    Computes whether an event contains duplicates based on index appearing more
    than once
    Args:
    - matched_indices = the combined matched_indices
    - max_index = how many truth final states there should be e.g. for 4tops, 12
    (4b + 8 decay)
    Returns:
        An ak.array with entry 1 if event contains duplicates else 0
    """
    
    max_event_size = ak.max(ak.count(matched_indices,axis=1))
    padded_indices = ak.fill_none(ak.pad_none(matched_indices,target=max_event_size),np.nan).to_numpy()
    
    index_appears_mask = [] 
    for i in range(max_index):
        index_appears_mask.append(np.count_nonzero(np.isin(padded_indices, [i]), axis=1) <= 1)
        
    return  ak.Array(1*~np.logical_and.reduce(index_appears_mask, axis=0))


def examine_matched_indices(list_of_matched_indices_arrays,reco_dict,max_index):
    
    all_decay_matched_indices = ak.concatenate(list_of_matched_indices_arrays,axis=1)
    
    reco_dict["contains_duplicates"] = compute_if_duplicate(all_decay_matched_indices,max_index)
    reco_dict["fully_matched"]       = compute_if_fully_matched(all_decay_matched_indices,max_index)
    
    return reco_dict

    
    