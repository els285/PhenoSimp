# Simulation Preparation Instructions

Chain for generating HyPER and Nu2-Flows `h5` inputs.

1. Generate DELPHES output ROOT files with `generate`
2. Parse and combine separate files into common format
3. Apply object-level and event-level selections, saving to new ROOT file
4. Parse with specific parsers