import subprocess as sub
from os import listdir, getcwd, makedirs
from os.path import isfile, join, getsize
from tqdm import tqdm
import argparse
from itertools import chain
import csv

# from mad4batch.ROOT_skim import skim_delphes
from mad4batch.delphes_branches import branches as DB


def write_job_script(output_file,branch_string):
    text = f"""#!/bin/bash
./skim_delphes $1 {output_file} {branch_string}
        """
    output_file = output_file.replace(".root","")
    job_name = f"job_{output_file}.sh"
    with open(job_name,"w") as file:
        file.write(text)

    
def write_submit_file(infile_template: str , outfile: str, Nfiles: int):
        
        output_file = outfile.replace(".root","")
        job_name = f"job_{output_file}.sh"
            
        text=f"""# Submit file for HTCondor
universe   = vanilla
executable = {job_name}
arguments = {infile_template}_$(Process).root
output     = $(ClusterId).$(Process).out
error      = $(ClusterId).$(Process).err
log        = $(ClusterId).$(Process).log
getenv = True
request_cpus = 12
request_memory = 10 GB
transfer_executable = True
should_transfer_files = YES
transfer_input_files    = ../Mad4Batch/src/mad4batch/skim_delphes, {infile_template}_$(Process).root
transfer_output_files = {outfile}
transfer_output_remaps = "{outfile} = {outfile.replace(".root","")}_$(Process).root"

when_to_transfer_output = ON_EXIT 

queue {Nfiles}"""

        name = f"skim_{outfile.replace('.root','')}.sub"
        with open(name,"w") as file:
            file.write(text)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Skim multiple Delphes ROOT file")
    parser.add_argument("-i","--input_file" , type=str, help="Input file template name",required=True)
    parser.add_argument("-o","--output_file", type=str, help="Path to the output ROOT file",required=True)
    parser.add_argument("-b","--branches", type=str, nargs='+', help="Branches to keep",required=False)
    parser.add_argument("--categories", type=str)
    parser.add_argument("--branch-file",type=str)
    parser.add_argument("-c","--config",type=str,required=False)
    parser.add_argument("--Nfiles",type=int,required=False)

    args = parser.parse_args()    
    
    #Check input files    
    mypath = "."
    files_to_skim = [f for f in listdir(mypath) if isfile(join(mypath, f)) and args.input_file in f and getsize(f)!=0][:args.Nfiles]
    print(f"Skimming {len(files_to_skim)} files")
    
    if args.categories is not None:
        if args.categories == ["High-Level"]:
            keys = ["Event","Jet","Muon","Electron","MissingET"]
        else:
            keys = args.branches
        branch_dict = {k: DB[k] for k in (DB.keys() & keys)}
        branches_to_keep  = list(chain(*list(branch_dict.values())))
        print("Retaining these branches")
        
    if args.branches is not None:
        branches_to_keep = args.branches    
        
    if args.config is not None:
        with open(args.config) as file:
            branches_to_keep = [line.rstrip() for line in file]
            
    branches_to_keep.remove('')
        
    # Print branches being kept    
    for br in branches_to_keep:
        print(f" - {br}")   

    branch_string = " ".join(branches_to_keep)

    write_job_script(args.output_file,branch_string)
    write_submit_file(args.input_file,args.output_file,len(files_to_skim))
