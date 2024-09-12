import subprocess as sub
from os import listdir, getcwd, makedirs
from os.path import isfile, join, getsize
from tqdm import tqdm
import argparse
from itertools import chain
import time

# from mad4batch.ROOT_skim import skim_delphes
from mad4batch.delphes_branches import branches as DB


def call_skimmer(infile, outname, branches_to_keep):
    executable = "/data/els285/MCproduction/Mad4Batch/src/mad4batch/skim_delphes"
    branch_string = " ".join(branches_to_keep)
    command = f"{executable} {infile} {outname} {branch_string}"
    sub.run(command,shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Skim multiple Delphes ROOT file")
    parser.add_argument("--output_file", type=str, help="Path to the output ROOT file",required=True)
    parser.add_argument("--branches", type=str, nargs='+', help="Branches to keep",required=False)
    parser.add_argument("--categories", type=str)
    parser.add_argument("--branch-file",type=str)
    parser.add_argument("--Nfiles",type=int,required=False)

    args = parser.parse_args()    
    
    
    if args.categories == ["High-Level"]:
        keys = ["Event","Jet","Muon","Electron","MissingET"]
    else:
        keys = args.branches
    branch_dict = {k: DB[k] for k in (DB.keys() & keys)}
    branches_to_keep  = list(chain(*list(branch_dict.values())))
    print("Retaining these branches")
        
    if args.branches is not None:
        branches_to_keep = args.branches    
        
    
    # Print branches being kept    
    for br in branches_to_keep:
        print(f" - {br}")
    
    mypath = "."
    files_to_hadd = [f for f in listdir(mypath) if isfile(join(mypath, f)) and "delphes" in f and getsize(f)!=0][:args.Nfiles]
    print(f"Skimming {args.Nfiles} files")

    t0 = time.time()
    for i,file in tqdm(enumerate(files_to_hadd)):
        outname = args.output_file.replace(".root",f"_{i}.root")
        call_skimmer(file, outname, branches_to_keep)
    print(time.time()-t0)
