import sys
import h5py  
import numpy as np
from pathlib import Path

def process(input_file:str, output_file:str , number:str = "all", shuffle:bool = False):
    
    input_path = Path(input_file)
    if not input_path.is_file():
        raise FileNotFoundError(f"{input_file} does not exist")
    
    output_path = Path(output_file)
    if output_path.is_file():
        raise FileExistsError(f"{input_file} already exists")
    
    # Open the input file in read mode and the output file in write mode
    with h5py.File(input_file, 'r') as infile, h5py.File(output_file, 'w') as outfile:
        
        # Get total number of events
        first_group = list(infile.items())[0][1]  # Get the first group in the file
        first_dataset = list(first_group.items())[0][1]  # Get the first dataset in the first group
        total_entries = first_dataset.shape[0]  # Assuming the first dimension is the one to shuffle
        
        if number=="all":
            N = total_entries
        else:
            N = int(number)
        # else: 
        #     raise ValueError("Specified input 3 should either be all or an integer")
        
        shuffle_indices = np.random.permutation(N)
        
        # Iterate over all the groups in the input file
        for group_name, group in infile.items():
            # Create a corresponding group in the output file
            out_group = outfile.create_group(group_name)
            
            # Iterate over all datasets in the group
            for dataset_name, dataset in group.items():
                data = dataset[()]
                total_entries = data.shape[0]  # Assuming the first dimension is the one you want to sample from
                
                # Check if N exceeds the total entries in the dataset
                if N > total_entries:
                    raise ValueError(f"N ({N}) cannot be greater than the total number of entries ({total_entries}) in dataset '{dataset_name}'.")
                
                # Slice the first N entries of the dataset
                sliced_data = data[:N]
                if shuffle:
                    out_data = sliced_data[shuffle_indices]
                else: 
                    out_data = sliced_data
                
                # Create a new dataset in the output file with the sliced data
                out_group.create_dataset(dataset_name, data=out_data)
                
if __name__=="__main__":
    # infile  = sys.argv[1]
    # outfile = sys.argv[2]
    # number  = sys.argv[3]
    # shuffle = sys.argv[4]
    
    # process(infile, outfile, number, shuffle)
    process(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])