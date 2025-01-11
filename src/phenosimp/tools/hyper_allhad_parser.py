import uproot
import numpy as np
import awkward as ak
import vector
import h5py 
import sys


def pad_variable(variable, max_len, pad_to = 0):
    padded_variable = ak.pad_none(variable, max_len, axis=1, clip=True)
    return ak.fill_none(padded_variable, pad_to)

def main():
    tree = uproot.open("4tops_ethan_allhad_forHyPER_0312.root:reco")
    pad_to_jet = int(ak.max(ak.count(tree["jet_pt"].array(),axis=1)))
    Nevents = len(tree["jet_pt"].array())
    
    # Jets
    jet_dt  = np.dtype([('e', np.float32), 
                        ('eta', np.float32), 
                        ('phi', np.float32), 
                        ('pt', np.float32), 
                        ('btag', np.int32), 
                        ('charge', np.float32)])
    jet_data = np.zeros((Nevents, pad_to_jet), dtype=jet_dt)

    jet_vectors = vector.zip({"pt"  : tree["jet_pt"].array(),
                            "eta" : tree["jet_eta"].array(),
                            "phi" : tree["jet_phi"].array(),
                            "m"   : tree["jet_mass"].array()})
    
    jet_data['pt']     = pad_variable(jet_vectors.pt  , pad_to_jet)
    jet_data['eta']    = pad_variable(jet_vectors.eta , pad_to_jet)
    jet_data['phi']    = pad_variable(jet_vectors.phi , pad_to_jet)
    jet_data['e']      = pad_variable(jet_vectors.e   , pad_to_jet)
    jet_data['btag']   = pad_variable(tree["jet_btag"].array() , pad_to_jet)
    jet_data['charge'] = np.zeros(Nevents).reshape(-1,1)
    
    # Globals 
    njets = ak.count(jet_vectors.pt,axis=1)
    nbjets = ak.count_nonzero(tree["jet_btag"].array()==1,axis=1)
    
    global_dt   = np.dtype([('njet', np.float32), ('nbTagged', np.float32)])
    global_data = np.zeros((Nevents, 1), dtype=global_dt)
    
    global_data['njet']     = njets.to_numpy().reshape(-1,1)
    global_data['nbTagged'] = nbjets.to_numpy().reshape(-1,1)
    
    # Matched Indices
    jet_indices = pad_variable(tree["jet_matched_indices"].array(),22,pad_to=-99).to_numpy()
    jet_indices[jet_indices==-9]=-8
    jet_indices[jet_indices==-99]=-9
    jet_truthmatch = jet_indices
    IndexSelect = tree["fully_matched"].array()

    h5_file = h5py.File(f"{outfile}", 'w')
    inputs_group = h5_file.create_group('INPUTS')
    labels_group = h5_file.create_group('LABELS')
    inputs_group.create_dataset("jet", data=jet_data)
    inputs_group.create_dataset("global", data=global_data)
    labels_group.create_dataset("VertexID", data=np.array(jet_truthmatch, dtype=np.int64))
    labels_group.create_dataset("IndexSelect", data = np.array(IndexSelect, dtype= np.int32))
    h5_file.close()
    
if __name__ == "__main__":
    infile = sys.argv[1]
    outfile= sys.argv[2]
    main(infile,outfile)

