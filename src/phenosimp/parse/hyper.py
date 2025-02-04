import uproot
import numpy as np
import awkward as ak
import vector
import h5py 
import argparse
import os

from phenosimp.tools.utils import pad_variable

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
                        # ('btag', np.int32), 
                        ('charge', np.float32)])
    jet_data = np.zeros((Nevents, pad_to_jet), dtype=jet_dt)

    jet_vectors = vector.zip({"pt"  : data_array["jet_pt"],
                            "eta" : data_array["jet_eta"],
                            "phi" : data_array["jet_phi"],
                            "m"   : data_array["jet_mass"]})
    
    jet_data['e']      = pad_variable(jet_vectors.e   , pad_to_jet)
    jet_data['eta']    = pad_variable(jet_vectors.eta , pad_to_jet)
    jet_data['phi']    = pad_variable(jet_vectors.phi , pad_to_jet)
    jet_data['pt']     = pad_variable(jet_vectors.pt  , pad_to_jet)
    # jet_data['btag']   = pad_variable(data_array["jet_btag"] , pad_to_jet, pad_to=np.nan)
    jet_data['charge'] = pad_variable(ak.zeros_like(data_array["jet_btag"]), pad_to_jet)
    
    jet_indices = pad_variable(data_array["jet_matched_indices"],       pad_to_jet,pad_to=-99).to_numpy()

    ################### LEPTONS #####################
    if do_individual_leptons or do_combined_leptons:
        
        # print("Parsing electron and muon separately")
        
        ############ ELECTRONS ####################
        Nelectrons = ak.num(data_array["el_pt"])
        pad_to_el = int(ak.max(Nelectrons))
            
        el_dt  = np.dtype([('e', np.float32), 
                            ('eta', np.float32), 
                            ('phi', np.float32), 
                            ('pt', np.float32), 
                            # ('btag', np.int32), 
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
        # el_data['btag']   = pad_variable(ak.zeros_like(el_vectors.pt), pad_to_el)
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
                            # ('btag', np.int32), 
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
        # mu_data['btag']   = pad_variable(ak.zeros_like(mu_vectors.pt), pad_to_mu)
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
                            # ('btag', np.int32), 
                            ('charge', np.float32)])
        lep_data = np.zeros((Nevents, pad_to_lep), dtype=lep_dt)
        
        # Combine the electron and muon 4-vectors horizontally and order by pt
        og_lep_vectors = ak.concatenate([el_vectors,mu_vectors],axis=1)
        lep_order = ak.argsort(og_lep_vectors.pt,ascending=False)
        lep_vectors = og_lep_vectors[lep_order]
        lep_charge  = ak.concatenate([data_array["el_charge"],data_array["mu_charge"]],axis=1)[lep_order]
        lep_indices = ak.concatenate([data_array["electron_matched_indices"],data_array["muon_matched_indices"]],axis=1)[lep_order]
        
        lep_data['e']      = pad_variable(lep_vectors.e   , pad_to_lep)
        lep_data['eta']    = pad_variable(lep_vectors.eta , pad_to_lep)
        lep_data['phi']    = pad_variable(lep_vectors.phi , pad_to_lep)
        lep_data['pt']     = pad_variable(lep_vectors.pt  , pad_to_lep)
        # lep_data['btag']   = pad_variable(ak.zeros_like(lep_charge), pad_to_lep)
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
    train_file = f"{outfile}_train.h5"
    test_file  = f"{outfile}_test.h5"

    with h5py.File(train_file, 'w') as h5_file:
        inputs_group = h5_file.create_group('INPUTS')
        labels_group = h5_file.create_group('LABELS')
        metadata_group = h5_file.create_group("METADATA")
        
        if do_individual_leptons:
            inputs_group.create_dataset("ELECTRON", data=el_data[train_mask])
            labels_group.create_dataset("ELECTRON", data=el_indices[train_mask])
            inputs_group.create_dataset("MUON",     data=mu_data[train_mask])
            labels_group.create_dataset("MUON",     data=mu_indices[train_mask])

        if do_combined_leptons:
            inputs_group.create_dataset("LEPTON",     data=lep_data[train_mask])
            labels_group.create_dataset("LEPTON",     data=lep_indices[train_mask])

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
            labels_group.create_dataset("ELECTRON", data=el_indices[test_mask])
            inputs_group.create_dataset("MUON",     data=mu_data[test_mask])
            labels_group.create_dataset("MUON",     data=mu_indices[test_mask])

        if do_combined_leptons:
            inputs_group.create_dataset("LEPTON",     data=lep_data[test_mask])
            labels_group.create_dataset("LEPTON",     data=lep_indices[test_mask])

        inputs_group.create_dataset("JET",          data=jet_data[test_mask])
        labels_group.create_dataset("JET",          data=jet_indices[test_mask])
        inputs_group.create_dataset("GLOBAL",       data=global_data[test_mask])
        
        metadata_group.create_dataset("FullyMatched",  data = np.array(FullyMatched[test_mask], dtype= np.int32))  
        
    print(f"File {test_file} written")
    

def create_directory_structure(parent_dir):
    # Create the parent directory with the name of 'outfile'
    raw_dir = os.path.join(parent_dir, "raw")
    processed_dir = os.path.join(parent_dir, "processed")
    
    # Check if the parent directory exists, if not create it
    if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)
        print(f"Created directory: {parent_dir}")
    
    # Create 'raw' and 'processed' subdirectories
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(processed_dir, exist_ok=True)
    
    # Print out confirmation
    print(f"Created subdirectories: {raw_dir} and {processed_dir}")
    
def parse_arguments():
    parser = argparse.ArgumentParser(description="Parse input data and handle output options.")
    
    # Positional argument for data_array (input data)
    parser.add_argument('infile', type=str, help='The output file where results will be saved.')
    parser.add_argument('outfile', type=str, help='The output file where results will be saved.')

    parser.add_argument('--individual_leptons', type=str, choices=['True', 'False'], default='False', 
                        help='Flag to include individual leptons (True or False).')
    parser.add_argument('--combined_leptons', type=str, choices=['True', 'False'], default='True', 
                        help='Flag to include combined leptons (True or False).')
    parser.add_argument('--neutrinos', type=str, choices=['True', 'False'], default='False', 
                        help='Flag to include neutrinos (True or False).')
    parser.add_argument('--split', type=int, default=10, 
                        help='Train-test split percentage (default is 10).')

    return parser.parse_args()

# Main logic to parse arguments and call function
if __name__ == '__main__':
    args = parse_arguments()
    
    from rich.console import Console
    from rich.table import Table

    # Create a console object to print output with rich formatting
    console = Console()
    
    # Convert 'True'/'False' strings to actual booleans
    do_individual_leptons = args.individual_leptons == 'True'
    do_combined_leptons = args.combined_leptons == 'True'
    do_neutrinos = args.neutrinos == 'True'
    
    outfile = args.outfile.replace(".h5","")
    
    console.print(
        f"Converting input file [bold green]4{args.infile}[/bold green] "
        f"to output [bold cyan]h5 file {args.outfile}[/bold cyan] "
        "with the following settings:",
        style="bold white")    

    # Create a table to display the information
    table = Table(title="Configuration Settings")

    # Add columns to the table
    table.add_column("Setting", style="bold green")
    table.add_column("Value", style="bold cyan")

    # Add rows to the table
    table.add_row("Individual charged leptons", str(do_individual_leptons))
    table.add_row("Combined charged leptons", str(do_combined_leptons))
    table.add_row("Neutrinos", str(do_neutrinos))
    table.add_row("Train-test split", str(args.split))

    # Print the table to the console
    console.print(table)
    
    train_file = f"{outfile}_train.h5"
    test_file  = f"{outfile}_test.h5"

    if os.path.exists(train_file) or os.path.exists(test_file):
        raise FileExistsError(f"Error: Files of type '{args.outfile}' already exist.")

    # Parse the inpute file
    tree = uproot.open(f"{args.infile}:Reco")
    all_branches = [br for br in tree.keys() if br[0]!="n"]
    event_array = tree.arrays(all_branches)
    
    # Create the directory
    console.print(
    f"Creating directory called [bold green]{outfile}[/bold green]",
    style="bold white")   
    create_directory_structure(outfile)

    # Call the function with parsed arguments
    parse(data_array=event_array,
        outfile=f"{outfile}/raw/{outfile}",
        do_individual_leptons=do_individual_leptons,
        do_combined_leptons=do_combined_leptons,
        do_neutrinos=do_neutrinos,
        train_test_split=args.split)
        