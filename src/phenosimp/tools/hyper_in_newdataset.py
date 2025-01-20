import uproot
import numpy as np
import awkward as ak
import vector
import h5py 
import sys

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


def pad_variable(variable, max_len, pad_to = 0):
    padded_variable = ak.pad_none(variable, max_len, axis=1, clip=True)
    return ak.fill_none(padded_variable, pad_to)


def parse(data_array                : ak.Array, 
          outfile                   : str, 
          do_individual_leptons     : str = False, 
          do_combined_leptons       : str = True,
          do_neutrinos              : str = False,
          train_test_split          : int = 10) -> None:
    
    """
    Takes a high-level awkward array with a selection pre-applied
    Writes the output to separate test and train h5 files
    Args:
        data_array - a high-level Awkward Array where the selection is pre-applied 
                        i.e. no further filtering will happen
        outfile    - the template name of the output h5 file
        do_individual_leptons - if electron and muon required separately
        do_combined_leptons - if charged leptons should be combined into one object
        do_neutrinos - if including neutrinos, presumably from a separate source
        train_test_split - integer ratio of number of (training events - 1):(testing events)
    Optional:

    """

    pad_to_jet = int(ak.max(ak.count(data_array["jet_pt"],axis=1)))
    Nevents = len(data_array["jet_pt"])
    
    print("Preparing data")
    
    # Globals 
    njets  = ak.count(data_array["jet_pt"],axis=1).to_numpy()
    nbjets = ak.count_nonzero(data_array["jet_btag"]==1,axis=1).to_numpy()
    
    global_dt   = np.dtype([('njet', np.float32), ('nbTagged', np.float32)])
    global_data = np.zeros((Nevents, 1), dtype=global_dt)
    
    global_data['njet']     = njets.reshape(-1,1)
    global_data['nbTagged'] = nbjets.reshape(-1,1)
    
    # Jets
    print("Parsing jets")
    jet_dt  = np.dtype([('e', np.float32), 
                        ('eta', np.float32), 
                        ('phi', np.float32), 
                        ('pt', np.float32), 
                        ('btag', np.int32), 
                        ('charge', np.float32),
                        ('id', np.float32)])
    jet_data = np.zeros((Nevents, pad_to_jet), dtype=jet_dt)

    jet_vectors = vector.zip({"pt"  : data_array["jet_pt"],
                            "eta" : data_array["jet_eta"],
                            "phi" : data_array["jet_phi"],
                            "m"   : data_array["jet_mass"]})
    
    jet_data['e']      = pad_variable(jet_vectors.e   , pad_to_jet)
    jet_data['eta']    = pad_variable(jet_vectors.eta , pad_to_jet)
    jet_data['phi']    = pad_variable(jet_vectors.phi , pad_to_jet)
    jet_data['pt']     = pad_variable(jet_vectors.pt  , pad_to_jet)
    jet_data['btag']   = pad_variable(data_array["jet_btag"] , pad_to_jet)
    jet_data['charge'] = np.zeros(Nevents).reshape(-1,1)
    jet_data['id']     =  1*(np.arange(pad_to_jet) < njets[:, None])
    
    jet_indices = pad_variable(data_array["jet_matched_indices"],       pad_to_jet,pad_to=-99).to_numpy()

    ################### LEPTONS #####################
    if do_individual_leptons or do_combined_leptons:
        
        print("Parsing electron and muon separately")
        
        ############ ELECTRONS ####################
        Nelectrons = ak.num(data_array["el_pt"])
        pad_to_el = int(ak.max(Nelectrons))
            
        el_dt  = np.dtype([('e', np.float32), 
                            ('eta', np.float32), 
                            ('phi', np.float32), 
                            ('pt', np.float32), 
                            ('btag', np.int32), 
                            ('charge', np.float32)])
        el_data = np.zeros((Nevents, pad_to_el), dtype=el_dt)

        el_vectors = vector.zip({"pt"  : data_array["el_pt"],
                                "eta" : data_array["el_eta"],
                                "phi" : data_array["el_phi"],
                                "m"   : data_array["el_mass"]})
        
        el_data['e']      = pad_variable(el_vectors.e   , pad_to_el)
        el_data['eta']    = pad_variable(el_vectors.eta , pad_to_el)
        el_data['phi']    = pad_variable(el_vectors.phi , pad_to_el)
        el_data['pt']     = pad_variable(el_vectors.pt  , pad_to_el)
        el_data['btag']   = np.zeros(Nevents).reshape(-1,1)
        el_data['charge'] = pad_variable(data_array["el_charge"] , pad_to_el)
        
        el_indices  = pad_variable(data_array["electron_matched_indices"],   pad_to_el,pad_to=-99).to_numpy()

        ################ MUONS ###################
        Nmuons = ak.num(data_array["mu_pt"])
        pad_to_mu = int(ak.max(Nmuons))
            
        # muectrons
        mu_dt  = np.dtype([('e', np.float32), 
                            ('eta', np.float32), 
                            ('phi', np.float32), 
                            ('pt', np.float32), 
                            ('btag', np.int32), 
                            ('charge', np.float32)])
        mu_data = np.zeros((Nevents, pad_to_mu), dtype=mu_dt)

        mu_vectors = vector.zip({"pt"  : data_array["mu_pt"],
                                "eta" : data_array["mu_eta"],
                                "phi" : data_array["mu_phi"],
                                "m"   : data_array["mu_mass"]})
        
        mu_data['e']      = pad_variable(mu_vectors.e   , pad_to_mu)
        mu_data['eta']    = pad_variable(mu_vectors.eta , pad_to_mu)
        mu_data['phi']    = pad_variable(mu_vectors.phi , pad_to_mu)
        mu_data['pt']     = pad_variable(mu_vectors.pt  , pad_to_mu)
        mu_data['btag']   = np.zeros(Nevents).reshape(-1,1)
        mu_data['charge'] = pad_variable(data_array["mu_charge"] , pad_to_mu)
        
        mu_indices  = pad_variable(data_array["muon_matched_indices"],  pad_to_mu,pad_to=-99).to_numpy()
        
    if do_combined_leptons:
                
        print("Parsing charged leptons as one object")
        Nleptons = ak.num(data_array["el_pt"]) + ak.num(data_array["mu_pt"])
        pad_to_lep = int(ak.max(Nleptons))

        lep_dt  = np.dtype([('e', np.float32), 
                            ('eta', np.float32), 
                            ('phi', np.float32), 
                            ('pt', np.float32), 
                            ('btag', np.int32), 
                            ('charge', np.float32)])
        lep_data = np.zeros((Nevents, pad_to_lep), dtype=lep_dt)
        
        # Combine the electron and muon 4-vectors horizontally and order by pt
        og_lep_vectors = ak.concatenate([el_vectors,mu_vectors],axis=1)
        lep_order = ak.argsort(og_lep_vectors.pt,ascending=False)
        lep_vectors = og_lep_vectors[lep_order]
        lep_charge  = ak.concatenate([data_array["el_charge"],data_array["mu_charge"]],axis=1)[lep_order]
        lep_indices = ak.concatenate([data_array["electorn_matched_indices"],data_array["muon_matched_indices"]],axis=1)[lep_order]
        
        lep_data['e']      = pad_variable(lep_vectors.e   , pad_to_lep)
        lep_data['eta']    = pad_variable(lep_vectors.eta , pad_to_lep)
        lep_data['phi']    = pad_variable(lep_vectors.phi , pad_to_lep)
        lep_data['pt']     = pad_variable(lep_vectors.pt  , pad_to_lep)
        lep_data['btag']   = np.zeros(Nevents).reshape(-1,1)
        lep_data['charge'] = pad_variable(lep_charge , pad_to_lep)
        
        lep_indices  = pad_variable(lep_indices,   pad_to_lep,pad_to=-99).to_numpy()
        
    FullyMatched = data_array["fully_matched"].to_numpy()
    
    """
    Writing to file
    """
    print("Writing to files")

    # Split into test and train    
    test_mask = np.arange(0,len(jet_data))%train_test_split==0
    train_mask = ~test_mask
    
    # Save to files
    train_file = outfile.replace(".h5","_train.h5")
    test_file  = outfile.replace(".h5","_test.h5")

    with h5py.File(train_file, 'w') as h5_file:
        inputs_group = h5_file.create_group('INPUTS')
        labels_group = h5_file.create_group('LABELS')
        metadata_group = h5_file.create_group("METADATA")
        
        if do_individual_leptons:
            inputs_group.create_dataset("ELECTRON", data=el_data[train_mask])
            labels_group.carate_dataset("ELECTRON", data=el_indices[train_mask])
            inputs_group.create_dataset("MUON",     data=mu_data[train_mask])
            labels_group.carate_dataset("MUON",     data=mu_indices[train_mask])

        if do_combined_leptons:
            inputs_group.create_dataset("LEPTON",     data=lep_data[train_mask])
            labels_group.carate_dataset("LEPTON",     data=lep_indices[train_mask])

        inputs_group.create_dataset("JET",          data=jet_data[train_mask])
        labels_group.create_dataset("JET",          data=jet_indices[train_mask])
        inputs_group.create_dataset("GLOBAL",       data=global_data[train_mask])
        
        metadata_group.create_dataset("FullyMatched",  data = np.array(FullyMatched[train_mask], dtype= np.int32))  
              
    print(f"File {train_file} written")
        
    with h5py.File(test_file, 'w') as h5_file:
        inputs_group = h5_file.create_group('INPUTS')
        labels_group = h5_file.create_group('LABELS')
        metadata_group = h5_file.create_group("METADATA")
        
        if do_individual_leptons:
            inputs_group.create_dataset("ELECTRON", data=el_data[test_mask])
            labels_group.carate_dataset("ELECTRON", data=el_indices[test_mask])
            inputs_group.create_dataset("MUON",     data=mu_data[test_mask])
            labels_group.carate_dataset("MUON",     data=mu_indices[test_mask])

        if do_combined_leptons:
            inputs_group.create_dataset("LEPTON",     data=lep_data[test_mask])
            labels_group.carate_dataset("LEPTON",     data=lep_indices[test_mask])

        inputs_group.create_dataset("JET",          data=jet_data[test_mask])
        labels_group.create_dataset("JET",          data=jet_indices[test_mask])
        inputs_group.create_dataset("GLOBAL",       data=global_data[test_mask])
        
        metadata_group.create_dataset("FullyMatched",  data = np.array(FullyMatched[test_mask], dtype= np.int32))  
        
    print(f"File {test_file} written")
        
def main(infile,outfile):
    
    """
    For parsing directly from a ROOT file 
    """
    tree = uproot.open(f"{infile}:Reco")
    all_branches = [br for br in tree.keys() if br[0]!="n"]
    event_array = tree.arrays(all_branches)
    parse(event_array)
    
if __name__ == "__main__":
    infile = sys.argv[1]
    outfile= sys.argv[2]
    main(infile,outfile)