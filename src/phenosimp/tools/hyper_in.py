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

def parse(data_array       : ak.Array, 
          outfile          : str, 
          do_leptons       : str = False, 
          do_neutrinos     : str = False,
          train_test_split : int = 10) -> None:
    
    """
    Takes a high-level awkward array with a selection pre-applied
    Writes the output to separate test and train h5 files
    Args:
        data_array - a high-level Awkward Array where the selection is pre-applied 
                        i.e. no further filtering will happen
        outfile    - the template name of the output h5 file
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

    
    if do_leptons:
        
        print("Parsing leptons")
        
        # Electrons
        Nelectrons = ak.num(data_array["el_pt"])
        pad_to_el = int(ak.max(Nelectrons))
            
        el_dt  = np.dtype([('e', np.float32), 
                            ('eta', np.float32), 
                            ('phi', np.float32), 
                            ('pt', np.float32), 
                            ('btag', np.int32), 
                            ('charge', np.float32),
                            ('id', np.float32)])
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
        el_data['id']     =  1*(np.arange(pad_to_el) < Nelectrons[:, None])
        
        el_indices  = pad_variable(data_array["electron_matched_indices"],   pad_to_el,pad_to=-99).to_numpy()

        
        # Muons
        Nmuons = ak.num(data_array["mu_pt"])
        pad_to_mu = int(ak.max(Nmuons))
            
        # muectrons
        mu_dt  = np.dtype([('e', np.float32), 
                            ('eta', np.float32), 
                            ('phi', np.float32), 
                            ('pt', np.float32), 
                            ('btag', np.int32), 
                            ('charge', np.float32),
                            ('id', np.float32)])
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
        mu_data['id']     =  1*(np.arange(pad_to_mu) < Nmuons[:, None])
        
        mu_indices  = pad_variable(data_array["muon_matched_indices"],  pad_to_mu,pad_to=-99).to_numpy()
        
    # Labels
    if do_leptons:
        VertexID    = ak.concatenate([jet_indices,el_indices,mu_indices],axis=1)
    else:
        VertexID    = jet_indices
        
    IndexSelect = data_array["fully_matched"].to_numpy()
    
    print("Writing to files")

    # Split into test and train    
    test_mask = np.arange(0,len(jet_data))%10==0
    train_mask = ~test_mask
    
    # Save to files
    train_file = outfile.replace(".h5","_train.h5")
    test_file  = outfile.replace(".h5","_test.h5")

    with h5py.File(train_file, 'w') as h5_file:
        inputs_group = h5_file.create_group('INPUTS')
        labels_group = h5_file.create_group('LABELS')
        
        if do_leptons:
            inputs_group.create_dataset("electron", data=el_data[train_mask])
            inputs_group.create_dataset("muon",     data=mu_data[train_mask])

        inputs_group.create_dataset("jet",          data=jet_data[train_mask])
        inputs_group.create_dataset("global",       data=global_data[train_mask])
        
        labels_group.create_dataset("VertexID",     data=np.array(VertexID[train_mask]))#, dtype=np.float32))
        labels_group.create_dataset("IndexSelect",  data = np.array(IndexSelect[train_mask], dtype= np.int32))  
              
    print(f"File {train_file} written")
        
    with h5py.File(test_file, 'w') as h5_file:
        inputs_group = h5_file.create_group('INPUTS')
        labels_group = h5_file.create_group('LABELS')
        
        if do_leptons:
            inputs_group.create_dataset("electron", data=el_data[test_mask])
            inputs_group.create_dataset("muon",     data=mu_data[test_mask])
        
        inputs_group.create_dataset("jet",          data=jet_data[test_mask])
        inputs_group.create_dataset("global",       data=global_data[test_mask])
        labels_group.create_dataset("VertexID",     data=np.array(VertexID[test_mask], dtype=np.int64))
        labels_group.create_dataset("IndexSelect",  data = np.array(IndexSelect[test_mask], dtype= np.int32))
        
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