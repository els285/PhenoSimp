import sys


def build_txt(parse_script,infile, Njobs):

    infile_edit = infile.replace(".root", "_$(Process).root")
    outfile = infile.replace(".root", "_$(Process)_parsed.root")

    txt = f"""# HTCondor submission script
    Universe       = vanilla
    Executable     = /home/els285/miniforge3/envs/madlad/bin/python
    Arguments      = -m phenosimp.delphes.{parse_script} {infile_edit} {outfile}
    transfer_executable = True
    Transfer_Input_Files = {infile_edit}
    Should_Transfer_Files = YES
    When_To_Transfer_Output = ON_EXIT
    Output         = logsB/job_$(Process).out
    Error          = logsB/job_$(Process).err
    Log            = logsB/job_$(Process).log
    Request_CPUs   = 2
    Request_Memory = 5GB
    Request_Disk   = 5GB

    Queue {Njobs}
    """
    
    return txt

def main(parse_script,infile,Njobs):
    
    txt = build_txt(parse_script,infile, Njobs)
    with open("condor_parse.sub", "w") as f:
        f.write(txt)
        
        
if __name__ == "__main__":
    
    parse_script = sys.argv[1]
    infile = sys.argv[2] 
    Njobs  = sys.argv[3]
    main(parse_script,infile,Njobs)