# PhenoSimp
Wrapper for tools I use for particle physics phenomenology.

`MadLAD` is used to generate MC simulation. The software run in custom singularity containers.

`condor` contains a set of wrappers for running MadLAD and
post-processessing MadLAD outputs via HTCondor.

`tools` contains scripts for downstream analysis of LHE, HepMC and ROOT files.

## Condor 

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
