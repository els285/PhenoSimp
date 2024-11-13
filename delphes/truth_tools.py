import uproot
import awkward as ak

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

    d["b_id"]     = tree["Particle.PID"].array()[status23 & ~has_a_W_mother & bottom_mask]
    # assert ak.count_nonzero(ak.count(d["b_id"],axis=1) !=4 ) == 0, "There are events with other than four bottom quarks"

    d["b_pt"]       = tree["Particle.PT"].array()[status23 & ~has_a_W_mother & bottom_mask]
    d["b_eta"]      = tree["Particle.Eta"].array()[status23 & ~has_a_W_mother & bottom_mask]
    d["b_phi"]      = tree["Particle.Phi"].array()[status23 & ~has_a_W_mother & bottom_mask]
    d["b_e"]        = tree["Particle.E"].array()[status23 & ~has_a_W_mother & bottom_mask]
    d["b_mass"]     = tree["Particle.Mass"].array()[status23 & ~has_a_W_mother & bottom_mask]

    # Wdecays
    print("Extracting W decay information")

    d["W_decay_id"]     = tree["Particle.PID"].array()[has_a_W_mother & is_fermion]
    d["W_decay_pt"]     = tree["Particle.PT"].array()[has_a_W_mother & is_fermion]
    d["W_decay_eta"]    = tree["Particle.Eta"].array()[has_a_W_mother & is_fermion]
    d["W_decay_phi"]    = tree["Particle.Phi"].array()[has_a_W_mother & is_fermion]
    d["W_decay_e"]      = tree["Particle.E"].array()[has_a_W_mother & is_fermion]
    d["W_decay_mass"]   = tree["Particle.Mass"].array()[has_a_W_mother & is_fermion]
    
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
    r["mu_pt"]      =  tree["Muon.Eta"].array()
    r["mu_pt"]      =  tree["Muon.Phi"].array()
    r["mu_charge"]  = tree["Muon.Charge"].array()

    r["met_met"]    =  tree["MissingET.MET"].array()
    r["met_eta"]    =  tree["MissingET.Eta"].array()
    r["met_phi"]    =  tree["MissingET.Phi"].array()
    
    return r