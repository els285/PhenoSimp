import yaml
import sys
import os 
from os import listdir
from os.path import isfile, join
from datetime import datetime
import argparse
from warnings import warn 
from numpy.random import randint

madlad_path = "/data/els285/MCproduction/MadLAD"

class Mad4Condor(object):
    
    
    def __init__(self,madladpath:str, config_name:str, cfg:dict, Njobs:int, lhe:bool, hepmc:bool):
        self.madlad_path    = madlad_path
        self.config_name    = config_name
        self.cfg            = cfg
        self.Njobs          = Njobs
        self.name           = self.cfg["gen"]["block_model"]["save_dir"]
        self.outputs_string = ""
        self.remaps_string  = ""
        self.out_lhe   = lhe 
        self.out_hepmc = hepmc
        
        # image = self.cfg["run"]["image"]
        # assert os.path.isfile(f"/MadLAD/{image}.sif"), "Singularity image not
        # found"
        if "block_delphes" not in list(self.cfg["gen"].keys()):
            warn("Delphes is turned off, so no ROOT file will be produced")
            
        
        if self.cfg["run"]["auto-launch"]==False:
            raise ValueError("auto-launch must be set to true for running with Condor")
        
    
    def create_directory(self):
        
        dt = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        self.condor_directory_name = f"{self.name}_condorrun_{dt}"
        assert not os.path.isdir(self.condor_directory_name), f"{self.condor_directory_name} is already a directory and will not be over-written"
        os.mkdir(self.condor_directory_name)
        

    def write_job_script(self):
        
        text=f"""#!/bin/bash

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data/zihanzhang/pkgs/miniforge3/bin/conda' 'shell.bash' 'hook' 2>/dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    . "/data/zihanzhang/pkgs/miniforge3/etc/profile.d/conda.sh" || export PATH="/data/zihanzhang/pkgs/miniforge3/bin:$PATH"
fi
unset __conda_setup
# <<< conda initialize <<<
conda activate madlad

# Read the random number from the first argument (provided by HTCondor)
ISEED=$1  # Pass the number from the input file

echo $ISEED

cd MadLAD   # Execute in MadLAD folder
python -m madlad.generate --config-name={self.config_name} gen.block_run.iseed=$ISEED
cd -        # Return to condor work directory

        """
        
        with open(f"{self.condor_directory_name}/job.sh","w") as file:
            file.write(text)
                
    def specify_outputs(self):
        
        evt_dir  = 'run_01_decayed_1'  if 'block_madspin' in list(self.cfg['gen'].keys()) else 'run_01'
        lhe_file = 'unweighted_events' if self.cfg["gen"]["block_model"]["order"]=="lo"   else 'events'
        hepmc_file = 'tag_1_pythia8_events.hepmc.gz' if self.cfg["gen"]["block_model"]["order"]=="lo" else 'events_PYTHIA8_0.hepmc'
        
        if self.out_lhe:
            self.outputs_string += f'MadLAD/{self.name}/Events/{evt_dir}/{lhe_file}.lhe.gz, '
            self.remaps_string  += f'{lhe_file}.lhe.gz = {lhe_file}_$(Cluster)_$(Process).lhe.gz; '
        
        if self.out_hepmc:
            if self.cfg["run"]["shower"] == True:
                self.outputs_string += f'MadLAD/{self.name}/Events/{evt_dir}/{hepmc_file}, '
                self.remaps_string  += f'{hepmc_file} = $(Cluster)_$(Process)_{hepmc_file}; '
                
            else: 
                print("Shower is turned off, will not transfer HepMC file or Delphes.root file")
                
        if 'block_delphes' in list(self.cfg['gen'].keys()):
            self.outputs_string += f" MadLAD/{self.name}.root, "
            self.remaps_string  += f"{self.name}.root = $(ClusterId)/{self.name}_$(Process).root; "
            
                
    def write_submit_file(self):
        
        text=f"""# Submit file for HTCondor
universe   = vanilla
executable = job.sh
arguments  = $(iseed)
output     = $(ClusterId)/logs/$(Process).out
error      = $(ClusterId)/logs/$(Process).err
log        = $(ClusterId)/logs/$(Process).log
request_cpus = 12
request_memory = 50 GB
transfer_executable = True
should_transfer_files = YES
transfer_input_files = {self.madlad_path}
transfer_output_files = {self.outputs_string}
transfer_output_remaps = "{self.remaps_string}"

when_to_transfer_output = ON_EXIT

# Queue jobs, passing each line from iseeds.txt
queue iseed from iseeds.txt
"""

        with open(f"{self.condor_directory_name}/submit.sub","w") as file:
            file.write(text)
            
            
    def generate_iseeds(self):
        
        """
        Generate iseeds.txt file with Njobs different seeds
        Iseed information: 
        https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/IntroGrid
        https://answers.launchpad.net/mg5amcnlo/+question/254698
        Max MadGraph iseed is 30081**2
        """
    
        iseeds = []
        with open(f"{self.condor_directory_name}/iseeds.txt", 'w') as f:

            for i in range(self.Njobs):
                low = i*2000
                f.write(f"{randint(low,low+1000)} \n")
             
def main(args):
         
    config_filepath = f"{madlad_path}/processes/{args.config}"
   
    with open(config_filepath) as stream:
        try:
            cfg = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
          
    config_name = config_filepath.split("/")[-1]
    
    RUN = Mad4Condor(madlad_path,config_name,cfg,args.Njobs,args.lhe,args.hepmc)
       
    RUN.create_directory()
    RUN.write_job_script()
    RUN.specify_outputs()
    RUN.write_submit_file()
    RUN.generate_iseeds()
    
    print(f"Directory {RUN.condor_directory_name} created   ")        

    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Generate HTCondor jobs for generation")
    parser.add_argument("--config", type=str, help="Config name, no need to point to the directory")
    parser.add_argument("--Njobs", type=int, help="Number of jobs")
    parser.add_argument("--lhe",action="store_true",required=False)
    parser.add_argument("--hepmc",action="store_true",required=False)
    
    args = parser.parse_args()   
    main(args)

