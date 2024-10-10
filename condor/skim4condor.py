import subprocess as sub
from os import listdir, getcwd, makedirs
from os.path import isfile, join, getsize
from tqdm import tqdm
import argparse
from itertools import chain
import csv

# from mad4batch.ROOT_skim import skim_delphes
from mad4condor.delphes_branches import branches as DB
from mad4condor.delphes_branches import light_branches

def write_job_script(path,output_file,branch_string):
    text = f"""#!/bin/bash
setupATLAS -q 
lsetup "root 6.30.02-x86_64-el9-gcc13-opt"
./skim_delphes $1 {output_file} {branch_string}
        """
    output_file = output_file.replace(".root","")
    job_name = f"job_{output_file}.sh"
    with open(f"{path}/{job_name}","w") as file:
        file.write(text)

    
def write_submit_file(path:str, infile_template: str , outfile: str, Nfiles: int):
        
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
request_cpus = 4
request_memory = 1 GB
transfer_executable = True
should_transfer_files = YES
transfer_input_files    = ../mad4condor/skim_delphes, {infile_template}_$(Process).root
transfer_output_files = {outfile}
transfer_output_remaps = "{outfile} = {outfile.replace(".root","")}_$(Process).root"

when_to_transfer_output = ON_EXIT 

queue {Nfiles}"""

        name = f"{path}/skim_{outfile.replace('.root','')}.sub"
        with open(name,"w") as file:
            file.write(text)

def run(args):
    
    #Check input files    
    mypath = args.directory 
    files_to_skim = []
    for fname in listdir(mypath):
        f = join(mypath,fname)
        if isfile(f) and "delphes" in f and getsize(f)!=0:
            files_to_skim.append(fname)
    if args.Nfiles:
        files_to_skim = files_to_skim[:args.Nfiles]
        
    infile_template  = args.infile or "_".join(files_to_skim[0].split("_")[:2])
    outfile_template = args.outfile or f"{infile_template}_skim.root"
    if outfile_template[-5:]!=".root":
        outfile_template = f"{outfile_template}.root"
    
    print(f"Skimming {len(files_to_skim)} files in directory {mypath} with name template {infile_template}")    

    ## Assigning branches
    if args.categories is not None:
        if args.categories=="HighLevel":
            keys = ["Event","Jet","Muon","Electron","MissingET"]
            branch_dict = {k: DB[k] for k in (DB.keys() & keys)}
            branches_to_keep  = list(chain(*list(branch_dict.values())))
        elif args.categories=="Light":
            branches_to_keep  = light_branches

    if args.branches is not None:
        branches_to_keep = args.branches    
        
    if args.config is not None:
        with open(args.config) as file:
            branches_to_keep = [line.rstrip() for line in file]
    
    print("Retaining branches")
    # Print branches being kept    
    for br in branches_to_keep:
        print(f" - {br}")   

    branch_string = " ".join(branches_to_keep)
    
    write_job_script(path=mypath,
                     output_file=outfile_template,
                     branch_string=branch_string)
    
    write_submit_file(path=mypath,
                      infile_template=infile_template,
                      outfile=outfile_template,
                      Nfiles=len(files_to_skim))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Skim multiple Delphes ROOT file")
    parser.add_argument("-d","--directory" , type=str, help="Directory",required=True)
    parser.add_argument("-f","--infile"  , type=str, help="Filename template" , required=False)
    parser.add_argument("-o","--outfile", type=str, help="Path to the output ROOT file",required=False)
    parser.add_argument("-b","--branches", type=str, nargs='+', help="Branches to keep",required=False)
    parser.add_argument("--categories", type=str)
    parser.add_argument("--branch-file",type=str)
    parser.add_argument("-c","--config",type=str,required=False)
    parser.add_argument("--Nfiles",type=int,required=False)

    args = parser.parse_args()    
    
    run(args)
