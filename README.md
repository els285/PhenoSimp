# PhenoSimp
Wrapper for tools used for particle physics phenomenology.

`MadLAD` is used to generate MC simulation. The software run in custom singularity containers.

`condor` contains a set of wrappers for running MadLAD and
post-processessing MadLAD outputs via HTCondor.

`lhe`contains methods for parsing LHE files into awkward structures.

`tools` contains scripts which don't fit into other categories

## Generation

The script `condor.generate` builds a directory containing HTCondor submission
script for running `MadLAD` in via the batch system.
It is run as
```bash
python -m condor.generate --config=<config_name> --Njobs=X --lhe --hepmc
```
where the last two flags are optional and provide the corresponding file formats
as outputs when the job completes.

### Skim
From `PhenoSimp` directory, run:
```
python -m condor.skim --directory <dir> <options>
```
Options:
* `infile` - template for Delphes output ROOT file inside directory. 
* `outfile` - output ROOT file name
* `Nfiles` - loops over a subset of files in the directory

Branch options:
* `categories` - `Light` = baseline high-level branches; `HighLevel` = all
  high-level branches. (The high-level categories are Jet, Electron, Muon,
  MissingET and Event)

Skim requires ROOT 6.30.02 which is set up in each Condor job.
