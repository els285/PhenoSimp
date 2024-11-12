import uproot
import awkward as ak
import vector

import sys

# Load
print("Loading file")
fname = sys.argv[1]
tree = uproot.open(f"{fname}:Delphes")

# Initialise dictionary
d = {}

print("Writing top and W data")
# Partons
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

status23    = tree["Particle.Status"].array()==23

# B quarks
print("Extracting b-quark information")
bottom_mask = abs(tree["Particle.PID"].array())==5
d["b_pt"]       = tree["Particle.PT"].array()[status23 & bottom_mask]
d["b_eta"]      = tree["Particle.Eta"].array()[status23 & bottom_mask]
d["b_phi"]      = tree["Particle.Phi"].array()[status23 & bottom_mask]
d["b_e"]        = tree["Particle.E"].array()[status23 & bottom_mask]
d["b_mass"]     = tree["Particle.Mass"].array()[status23 & bottom_mask]
d["b_id"]       = tree["Particle.PID"].array()[status23 & bottom_mask]


# Hadronic W decays
"""
Hadronic W decays are taken from status==23 with PDGID filter

"""
print("Extracting hadronic W decay information")
only_quarks_mask = abs(tree["Particle.PID"].array()[status23])<=6

W_hadronic_decays_id    = tree["Particle.PID"].array()[status23][only_quarks_mask][:,4:]
W_hadronic_decays_pt    = tree["Particle.PT"].array()[status23][only_quarks_mask][:,4:]
W_hadronic_decays_eta   = tree["Particle.Eta"].array()[status23][only_quarks_mask][:,4:]
W_hadronic_decays_phi   = tree["Particle.Phi"].array()[status23][only_quarks_mask][:,4:]
W_hadronic_decays_e     = tree["Particle.E"].array()[status23][only_quarks_mask][:,4:]
W_hadronic_decays_m     = tree["Particle.Mass"].array()[status23][only_quarks_mask][:,4:]


# Leptons
"""
Leptons are found based on PID and having a mother with PID==24
"""
print("Extracting leptonic W decay information")

# Define a lepton mask based on PID
m11 = abs(tree["Particle.PID"].array())==11
m12 = abs(tree["Particle.PID"].array())==12
m13 = abs(tree["Particle.PID"].array())==13
m14 = abs(tree["Particle.PID"].array())==14
m15 = abs(tree["Particle.PID"].array())==15
m16 = abs(tree["Particle.PID"].array())==16
all_lepton_mask = m11 | m12 | m13 | m14 | m15 | m16

# Filter for only leptons which come from Ws
lepton_mothers = tree["Particle.M1"].array()[all_lepton_mask]
leptons_from_W = abs(tree["Particle.PID"].array()[lepton_mothers])==24

W_leptonic_decays_id    = tree["Particle.PID"].array()[all_lepton_mask][leptons_from_W]
W_leptonic_decays_pt    = tree["Particle.PT"].array()[all_lepton_mask][leptons_from_W]
W_leptonic_decays_eta   = tree["Particle.Eta"].array()[all_lepton_mask][leptons_from_W]
W_leptonic_decays_phi   = tree["Particle.Phi"].array()[all_lepton_mask][leptons_from_W]
W_leptonic_decays_e     = tree["Particle.E"].array()[all_lepton_mask][leptons_from_W]
W_leptonic_decays_m     = tree["Particle.Mass"].array()[all_lepton_mask][leptons_from_W]


# Combine hadronic and leptonic decay products
Wdecay_combined_id  = ak.concatenate([W_hadronic_decays_id,W_leptonic_decays_id],axis=1)
Wdecay_combined_pt  = ak.concatenate([W_hadronic_decays_pt,W_leptonic_decays_pt],axis=1)
Wdecay_combined_eta = ak.concatenate([W_hadronic_decays_eta,W_leptonic_decays_eta],axis=1)
Wdecay_combined_phi = ak.concatenate([W_hadronic_decays_phi,W_leptonic_decays_phi],axis=1)
Wdecay_combined_e   = ak.concatenate([W_hadronic_decays_e,W_leptonic_decays_e],axis=1)
Wdecay_combined_m   = ak.concatenate([W_hadronic_decays_m,W_leptonic_decays_m],axis=1)

print("Performing delta-R matching")
# Build vectors of W bosons
"""
The best performing matching is through delta-R matching of candidate Ws to the
true Ws found through status==22
{W1,W2,W3} are candidate W bosons built from the decay products
"""
W1 = vector.zip({"pt":Wdecay_combined_pt[:,0],
                 "eta":Wdecay_combined_eta[:,0],
                 "phi":Wdecay_combined_phi[:,0],
                 "e":Wdecay_combined_e[:,0]}) + vector.zip({"pt":Wdecay_combined_pt[:,1],
                 "eta":Wdecay_combined_eta[:,1],
                 "phi":Wdecay_combined_phi[:,1],
                 "e":Wdecay_combined_e[:,1]})
                 
W2 = vector.zip({"pt":Wdecay_combined_pt[:,2],
                 "eta":Wdecay_combined_eta[:,2],
                 "phi":Wdecay_combined_phi[:,2],
                 "e":Wdecay_combined_e[:,2]}) + vector.zip({"pt":Wdecay_combined_pt[:,3],
                 "eta":Wdecay_combined_eta[:,3],
                 "phi":Wdecay_combined_phi[:,3],
                 "e":Wdecay_combined_e[:,3]})
                 
W3 = vector.zip({"pt":Wdecay_combined_pt[:,4],
                 "eta":Wdecay_combined_eta[:,4],
                 "phi":Wdecay_combined_phi[:,4],
                 "e":Wdecay_combined_e[:,4]}) + vector.zip({"pt":Wdecay_combined_pt[:,5],
                 "eta":Wdecay_combined_eta[:,5],
                 "phi":Wdecay_combined_phi[:,5],
                 "e":Wdecay_combined_e[:,5]})
                 
dW1 = vector.zip({"pt":d["W_pt"][:,0],                  
                 "eta":d["W_eta"][:,0],
                 "phi":d["W_phi"][:,0],
                 "e":d["W_e"][:,0]}) 

dW2 = vector.zip({"pt":d["W_pt"][:,1],                  
                 "eta":d["W_eta"][:,1],
                 "phi":d["W_phi"][:,1],
                 "e":d["W_e"][:,1]})

dW3 = vector.zip({"pt":d["W_pt"][:,2],                  
                 "eta":d["W_eta"][:,2],
                 "phi":d["W_phi"][:,2],
                 "e":d["W_e"][:,2]})  

DR_dW1 = ak.concatenate([ak.unflatten(dW1.deltaR(W1),1,axis=0),
                ak.unflatten(dW1.deltaR(W2),1,axis=0),
                ak.unflatten(dW1.deltaR(W3),1,axis=0)],axis=1)

DR_dW2 = ak.concatenate([ak.unflatten(dW2.deltaR(W1),1,axis=0),
                ak.unflatten(dW2.deltaR(W2),1,axis=0),
                ak.unflatten(dW2.deltaR(W3),1,axis=0)],axis=1)

DR_dW3 = ak.concatenate([ak.unflatten(dW3.deltaR(W1),1,axis=0),
                ak.unflatten(dW3.deltaR(W2),1,axis=0),
                ak.unflatten(dW3.deltaR(W3),1,axis=0)],axis=1)


print("Computing indices")
min_DR_dW1 = ak.argmin(abs(DR_dW1),axis=1)
min_DR_dW2 = ak.argmin(abs(DR_dW2),axis=1)
min_DR_dW3 = ak.argmin(abs(DR_dW3),axis=1)

indices = ak.concatenate([ak.unflatten(min_DR_dW1,1,axis=0),
                          ak.unflatten(min_DR_dW2,1,axis=0),
                          ak.unflatten(min_DR_dW3,1,axis=0)],axis=1)

"""
The delta-R matching is not always successful: matched branch indiciates whether
there is a unique set of indices
"""

"""
Sorting the indices to reflect the desired order
"""

final_indices = ak.concatenate([2*indices,2*indices+1],axis=1)[:,[0,3,1,4,2,5]]

d["W_decay_pt"]       = ak.fill_none(Wdecay_combined_pt[final_indices],float("nan"))
d["W_decay_eta"]      = ak.fill_none(Wdecay_combined_eta[final_indices],float("nan"))
d["W_decay_phi"]      = ak.fill_none(Wdecay_combined_phi[final_indices],float("nan"))
d["W_decay_e"]        = ak.fill_none(Wdecay_combined_e[final_indices],float("nan"))
d["W_decay_mass"]     = ak.fill_none(Wdecay_combined_m[final_indices],float("nan"))
d["W_decay_id"]       = ak.fill_none(Wdecay_combined_id[final_indices],float("nan"))

d["EventNumber"]      = tree["Event.Number"].array()

d["duplicate_matched"]=ak.fill_none(((ak.sum(indices,axis=1))!=3)*1,float("nan"))


all_minima = ak.concatenate([DR_dW1[ak.unflatten(indices[:,0],1,axis=0)],
                             DR_dW2[ak.unflatten(indices[:,1],1,axis=0)],
                             DR_dW3[ak.unflatten(indices[:,2],1,axis=0)]],axis=1)


d["greater_than_0p4"] = ak.count_nonzero(all_minima>0.4,axis=1)


print("Truth particle parsing complete")


### Reco-level
print("Writing reco-level trees")
r = {}

r["EventNumber"]      = tree["Event.Number"].array()

r["jet_pt"]   =  tree["Jet.PT"].array()
r["jet_eta"]  =  tree["Jet.Eta"].array()
r["jet_phi"]  =  tree["Jet.Phi"].array()
r["jet_mass"] =  tree["Jet.Mass"].array()
r["jet_btag"] = tree["Jet.BTag"].array()

r["el_pt"] =  tree["Electron.PT"].array()
r["el_eta"] =  tree["Electron.Eta"].array()
r["el_phi"] =  tree["Electron.Phi"].array()
r["el_charge"] =  tree["Electron.Charge"].array()

r["mu_pt"] =  tree["Muon.PT"].array()
r["mu_pt"] =  tree["Muon.Eta"].array()
r["mu_pt"] =  tree["Muon.Phi"].array()
r["mu_charge"] = tree["Muon.Charge"].array()

r["met_met"]  =  tree["MissingET.MET"].array()
r["met_eta"]  =  tree["MissingET.Eta"].array()
r["met_phi"]  =  tree["MissingET.Phi"].array()

print("Writing file")
fname = sys.argv[2]
with  uproot.recreate(f"{fname}") as file:
    file["Truth"] = d
    file["Reco"]  = r
print("Complete")