import awkward as ak 
import uproot 
import numpy as np
import vector
import h5py
import sys


def pad_variable(variable, max_len, pad_to = 0):
    padded_variable = ak.pad_none(variable, int(max_len), axis=1, clip=True)
    return ak.fill_none(padded_variable, pad_to)

def parse_jets(reco,selection):
    
    """
    Build jet vectors and apply the selection to them
    Build jet_data object
    """
    
    jet_vectors = vector.zip({"pt":reco["jet_pt"].array(),
                          "eta":reco["jet_eta"].array(),
                          "phi":reco["jet_phi"].array(),
                          "m":reco["jet_mass"].array()})[selection]

    jet_dt  = np.dtype([('pt', np.float32), 
                        ('eta', np.float32), 
                        ('phi', np.float32), 
                        ('energy', np.float32),
                        ('is_tagged', np.float32)])
    
    max_jets = ak.max(ak.count(reco["jet_pt"].array(),axis=1))
    jet_data = np.zeros((ak.count_nonzero(selection), max_jets), dtype=jet_dt)

    # Jets
    jet_data["pt"]          = pad_variable(jet_vectors.pt           , max_jets)
    jet_data["eta"]         = pad_variable(jet_vectors.eta          , max_jets)
    jet_data["phi"]         = pad_variable(jet_vectors.phi          , max_jets)
    jet_data["energy"]      = pad_variable(jet_vectors.e            , max_jets)
    jet_data["is_tagged"]   = pad_variable(reco["jet_btag"].array()[selection] , max_jets)
    
    jet_data = jet_data.astype([('pt', '<f4'), ('eta', '<f4'), ('phi', '<f4'), ('energy', '<f4'), ('is_tagged', '<?')])

    return jet_data


def parse_leptons(reco,selection):
    
    """
    Parses leptons - still needs the lepton type
    """
    
    electron_vectors = vector.zip({"pt":reco["el_pt"].array(),
                          "eta":reco["el_eta"].array(),
                          "phi":reco["el_phi"].array(),
                          "m":0.5110e-3})[selection]


    muon_vectors = vector.zip({"pt":reco["mu_pt"].array(),
                          "eta":reco["mu_eta"].array(),
                          "phi":reco["mu_phi"].array(),
                          "m":105.66e-3})[selection]
    
    lepton_vectors = ak.concatenate([electron_vectors,muon_vectors],axis=1)
    lepton_charge  = ak.concatenate([reco["el_charge"].array(),reco["mu_charge"].array()],axis=1)
    lepton_type    = ak.concatenate([abs(reco["el_charge"].array())*0,abs(reco["mu_charge"].array())*1],axis=1)

    lepton_pTs = ak.concatenate([reco["el_pt"].array(),reco["mu_pt"].array()],axis=1)[selection]
    sorted_indices =  ak.argsort(lepton_pTs,ascending=False)
    
    lepton_vectors = lepton_vectors[sorted_indices]
    lepton_charge  = lepton_charge[selection][sorted_indices]
    lepton_type    = lepton_type[selection][sorted_indices]
    
    leptons_dt = np.dtype( [("pt"       , np.float64),
                            ("eta"      , np.float64),
                            ("phi"      , np.float64),
                            ("energy"   , np.float64),
                            ("charge"   , np.float64),
                            ("type"     , np.float64)])
                          
    Nleptons = ak.max(ak.count(lepton_pTs,axis=1))

    # Leptons
    lepton_data = np.zeros((ak.count_nonzero(selection), Nleptons), dtype=leptons_dt)
    lepton_data["pt"]     = pad_variable(lepton_vectors.pt , Nleptons)
    lepton_data["eta"]    = pad_variable(lepton_vectors.eta , Nleptons)
    lepton_data["phi"]    = pad_variable(lepton_vectors.phi , Nleptons)
    lepton_data["energy"] = pad_variable(lepton_vectors.e   , Nleptons)
    lepton_data["charge"] = pad_variable(lepton_charge      , Nleptons)
    lepton_data["type"]   = pad_variable(lepton_type      , Nleptons)


    return lepton_data


def parse_neutrinos(truth,selection):
    
    truth_leptonic_mask = abs(truth["W_decay_id"].array())>10
    up_type_mask        = truth["W_decay_id"].array()%2==0
    
    neutrino_vectors = vector.zip({"pt":truth["W_decay_pt"].array(),
                          "eta":truth["W_decay_eta"].array(),
                          "phi":truth["W_decay_phi"].array(),
                          "m":0})[truth_leptonic_mask*up_type_mask]
    
    neutrino_vectors = neutrino_vectors[selection]
    
    neutrino_ids = truth["W_decay_id"].array()[truth_leptonic_mask*up_type_mask][selection]
    
    neutrino_ordering = ak.argsort(neutrino_ids,ascending=False) #This might break, it's to ensure the positive IDs come first (i.e.e the neutrinos produced alongside a negative lepton)
    
    sorted_neutrino_ids     = neutrino_ids[neutrino_ordering]
    sorted_neutrino_vectors = neutrino_vectors[neutrino_ordering]
    
    N_neutrinos = ak.max(ak.count(neutrino_ids,axis=1))

    neutrino_dt = np.dtype( { 'names':["PDGID","pt","eta","phi","mass"],
                             'formats':(np.float64, np.float64, np.float64, np.float64, np.float64, np.float64) } )

    # This is to be edited
    neutrino_data = np.zeros((ak.count_nonzero(selection), N_neutrinos), dtype=neutrino_dt)
    neutrino_data["pt"]     = pad_variable(sorted_neutrino_vectors.pt,  N_neutrinos)
    neutrino_data["eta"]    = pad_variable(sorted_neutrino_vectors.eta, N_neutrinos)
    neutrino_data["phi"]    = pad_variable(sorted_neutrino_vectors.phi, N_neutrinos)
    neutrino_data["mass"]   = pad_variable(sorted_neutrino_vectors.m,   N_neutrinos)
    
    neutrino_data["PDGID"] = pad_variable(sorted_neutrino_ids, N_neutrinos)
    
    return neutrino_data

def parse_MET(reco,selection):
    
    met_dt = np.dtype( { 'names':["MET","phi"],
                             'formats':(np.float64, np.float64) } )

    met_data = np.zeros((ak.count_nonzero(selection)), dtype=met_dt)
    met_data["MET"] = reco["met_met"].array().to_numpy()[selection][:,0]
    met_data["phi"] = reco["met_phi"].array().to_numpy()[selection][:,0]
    
    return met_data
    

def parse_event_level(reco,selection):
    
    """
    This assumes that the truth and reco trees are already matched
    """
    # Event-level data
    njets  = ak.count(reco["jet_btag"].array()[selection],axis=1)
    nbjets = ak.count_nonzero(reco["jet_btag"].array()[selection],axis=1)
    eventNumber = reco["EventNumber"].array()[selection][:,0]
    
    return njets,nbjets,eventNumber


def main(infile_name,outfile_name,split_ratio):
    
    logging.info("Reading input ROOT file and splitting into various lepton multiplicity channels")
    reco  = uproot.open(f"{infile_name}:Reco")
    truth = uproot.open(f"{infile_name}:Truth")
        
    Nelectrons = ak.count(reco["el_pt"].array(),axis=1)
    Nmuons     = ak.count(reco["mu_pt"].array(),axis=1)

    reco_1L_mask = (Nelectrons+Nmuons)==1
    reco_2L_mask = (Nelectrons+Nmuons)==2
    reco_3L_mask = (Nelectrons+Nmuons)==3
    reco_4L_mask = (Nelectrons+Nmuons)==4

    truth_leptonic_mask = abs(truth["W_decay_id"].array())>10

    truth_1L_mask = ak.count(truth["W_decay_id"].array()[truth_leptonic_mask],axis=1)==2
    truth_2L_mask = ak.count(truth["W_decay_id"].array()[truth_leptonic_mask],axis=1)==4
    truth_3L_mask = ak.count(truth["W_decay_id"].array()[truth_leptonic_mask],axis=1)==6
    truth_4L_mask = ak.count(truth["W_decay_id"].array()[truth_leptonic_mask],axis=1)==8

    tau_event_mask = ak.count_nonzero(abs(truth["W_decay_id"].array())==15,axis=1)!=0

    full_1L_event = reco_1L_mask*truth_1L_mask*~tau_event_mask
    full_2L_event = reco_2L_mask*truth_2L_mask*~tau_event_mask
    full_3L_event = reco_3L_mask*truth_3L_mask*~tau_event_mask
    full_4L_event = reco_4L_mask*truth_4L_mask*~tau_event_mask

    print("Parsing event-level data")
    njets , nbjets, eventNumber = parse_event_level(reco,full_1L_event)
    
    print("Parsing jet data")
    jet_data      = parse_jets(reco,full_1L_event)
    
    print("Parsing lepton data")
    lepton_data   = parse_leptons(reco,full_1L_event)
    
    print("Parsing neutrino data")
    neutrino_data = parse_neutrinos(truth,full_1L_event)
    
    print("Parsing MET data")
    met_data      = parse_MET(reco,full_1L_event)
    
    if any([len(x)!=jet_data.shape[0] for x in [lepton_data,neutrino_data,met_data]]):
        raise ValueError("The arrays are of different length")
    
    
    hf_train = h5py.File(f"{outfile_name}_train.h5", 'w')
    hf_test  = h5py.File(f"{outfile_name}_test.h5" , 'w')
    
    test_mask = np.arange(0,len(jet_data))%split_ratio==0
    train_mask = ~test_mask
    
    hf_train.create_group("data")
    hf_test.create_group("data")
    
    hf_train["data"].create_dataset('jets',data=jet_data[train_mask])
    hf_test["data"].create_dataset('jets',data=jet_data[test_mask])
    
    hf_train["data"].create_dataset('leptons',data=lepton_data[train_mask])
    hf_test["data"].create_dataset('leptons',data=lepton_data[test_mask])
    
    hf_train["data"].create_dataset('neutrinos',data=neutrino_data[train_mask])
    hf_test["data"].create_dataset('neutrinos',data=neutrino_data[test_mask])
    
    hf_train["data"].create_dataset('MET',data=met_data[train_mask])
    hf_test["data"].create_dataset('MET',data=met_data[test_mask])
    
    hf_train["data"].create_dataset('njets',data=njets[train_mask])
    hf_test["data"].create_dataset('njets',data=njets[test_mask])
    
    hf_train["data"].create_dataset('nbjets',data=nbjets[train_mask])
    hf_test["data"].create_dataset('nbjets',data=nbjets[test_mask])
    
    hf_train["data"].create_dataset('eventNumber',data=eventNumber[train_mask])
    hf_test["data"].create_dataset('eventNumber',data=eventNumber[test_mask])
    

if __name__ == "__main__":
    main(infile_name=sys.argv[1],
         outfile_name=sys.argv[2],
         split_ratio=sys.argv[3])