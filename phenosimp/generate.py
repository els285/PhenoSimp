import yaml
import sys
import os 
from os import listdir
from os.path import isfile, join
from datetime import datetime
import argparse


class Mad4Condor(object):
    
    def __init__(self,config_name,cfg,Njobs,lhe,hepmc):
        self.config_name    = config_name
        self.cfg            = cfg
        self.Njobs          = Njobs
        self.name           = self.cfg["gen"]["block_model"]["save_dir"]
        self.outputs_string = ""
        self.remaps_string  = ""
        self.out_lhe   = lhe 
        self.out_hepmc = hepmc
        
        self.create_directory()
        self.write_job_script()
        self.specify_outputs()
        self.write_submit_file()        

    
    def create_directory(self):
        
        dt = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        self.condor_directory_name = f"{self.name}_condorrun_{dt}"
        assert not os.path.isdir(self.condor_directory_name), f"{self.condor_directory_name} is already a directory and will not be over-written"
        os.mkdir(self.condor_directory_name)
        

    def write_job_script(self):
        text = f"""#!/bin/bash

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data/zihanzhang/pkgs/miniforge3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data/zihanzhang/pkgs/miniforge3/etc/profile.d/conda.sh" ]; then
        . "/data/zihanzhang/pkgs/miniforge3/etc/profile.d/conda.sh"
    else
        export PATH="/data/zihanzhang/pkgs/miniforge3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
conda activate madlad

cd MadLAD   # Execute in MadLAD folder
python -m madlad.generate --config-name={self.config_name} gen.block_run.iseed=$RANDOM
cd -        # Return to condor work directory

        """
        
        with open(f"{self.condor_directory_name}/job.sh","w") as file:
            file.write(text)
                
    def specify_outputs(self):
        
        evt_dir  = 'run_01_decayed_1'  if 'block_madspin' in list(self.cfg['gen'].keys()) else 'run_01'
        lhe_file = 'unweighted_events' if self.cfg["gen"]["block_model"]["order"]=="lo"   else 'events'
        
        if self.out_lhe:
            self.outputs_string += f'MadLAD/{self.name}/Events/{evt_dir}/{lhe_file}.lhe.gz, '
            self.remaps_string  += f'{lhe_file}.lhe.gz = {lhe_file}_$(Cluster)_$(Process).lhe.gz; '
        
        if self.out_hepmc:
            if self.cfg["run"]["shower"] == True:
                self.outputs_string += f'MadLAD/{self.name}/Events/{evt_dir}/events_PYTHIA8_0.hepmc, '
                self.remaps_string  += f'events_PYTHIA8_0.hepmc = events_$(Cluster)_$(Process).hepmc; '
                
            else: 
                print("Shower is turned off, will not transfer HepMC file or Delphes.root file")
                
        if 'block_delphes' in list(self.cfg['gen'].keys()):
            self.outputs_string += f" MadLAD/{self.name}.root, "
            self.remaps_string  += f"{self.name}.root = delphes_$(Cluster)_$(Process).root; "
            

                

    def write_submit_file(self):
        
        text=f"""# Submit file for HTCondor
universe   = vanilla
executable = job.sh
arguments  = $(RandomNumber)
output     = $(ClusterId).$(Process).out
error      = $(ClusterId).$(Process).err
log        = $(ClusterId).$(Process).log
request_cpus = 12
request_memory = 50 GB
transfer_executable = True
should_transfer_files = YES
transfer_input_files    = ../MadLAD
transfer_output_files = {self.outputs_string}
transfer_output_remaps = "{self.remaps_string}"

when_to_transfer_output = ON_EXIT 

queue {self.Njobs}"""

        with open(f"{self.condor_directory_name}/submit.sub","w") as file:
            file.write(text)
            
       
def main():
    
    parser = argparse.ArgumentParser(description="Skim multiple Delphes ROOT file")
    parser.add_argument("output_file", type=str, help="Path to the output ROOT file")
    parser.add_argument("Njobs", type=int, help="Branches to keep")
    parser.add_argument("--lhe",action="store_true",required=False)
    parser.add_argument("--hepmc",action="store_true",required=False)
    
    args = parser.parse_args()   
         
    config_filepath = f"MadLAD/processes/{args.output_file}"
   
    with open(config_filepath) as stream:
        try:
            cfg = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
          
    config_name = config_filepath.split("/")[-1]
    
    RUN = Mad4Condor(config_name,cfg,args.Njobs,args.lhe,args.hepmc)
    
main()
