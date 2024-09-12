
import numpy as np
import pylhe
import awkward as ak

class LHEparse:

    """
    Class for parsing LHE files to structed awkward arrays
    Particles are identified by the pdgid number as defined in the class
    dictionary PDGID
    Args:
    - 
    """
      
    PDGID = {
    1.0:"down"                  ,   -1.0:"anti_down",
    2.0:"up"                    ,   -2.0:"anti_up",
    3.0:"strange"               ,   -3.0:"anti_strange",
    4.0:"charm"                 ,   -4.0:"anti_charm",
    5.0:"bottom"                ,   -5.0:"anti_bottom",
    6.0:"top"                  ,   -6.0:"anti_top",

    7.0:"bprime"                ,   -7.0:"bprimebar",
    8.0:"tprime"                ,   -8.0:"tprimebar",

    11.0:"electron"             ,   -11.0:"anti_electron",
    12.0:"electron_neutrino"    ,   -12.0:"anti_electron_neutrino",
    13.0:"muon"                 ,   -13.0:"anti_muon",
    14.0:"muon_nuetrino"        ,   -14.0:"anti_muon_neutrino",
    15.0:"tau"                  ,   -15.0:"anti_tau",
    16.0:"tau_neutrino"         ,   -16.0:"anti_tau_neutirno",

    17.0:"tauprimeP"            ,   -17.0:"tauprimeM",
    18.0:"v_tauprimeP"          ,   -18.0:"v_tauprimeM",

    21.0:"gluon"          , 
    22.0:"photon"         , 
    23.0:"Z"              , 
    24.0:"W_plus"               ,   -24:"W_minus",
    25.0:"higgs"      
    }



    def __init__(self , file_name:str):
        self.file_name      = file_name
        print(f"Parsing LHE file {self.file_name}")
        self.array            =  pylhe.to_awkward(pylhe.read_lhe_with_attributes(self.file_name))


    def build(self):
        
        """
        Builds a flat awkward event array indexed by the particles names, 
            as given by the values of the PDGID dict.
        Supplemental composite fields for combined quark- and lepton-types.
        """
        
        # Parse only the PDGIDs which exist in the imported array
        unique_keys    = np.unique(self.array.particles.id.to_numpy())
        PDGID_filtered = {k:v for k,v in self.PDGID.items() if k in unique_keys}
        
        columns = {}
        for k,v in PDGID_filtered.items():
            columns[v] = self.array.particles[self.array.particles.id==k]
            
        columns["up_type_quarks"]        = ak.concatenate([columns["up"],columns["charm"],columns["top"]],axis=1)
        columns["anti_up_type_quarks"]   = ak.concatenate([columns["anti_up"],columns["anti_charm"],columns["anti_top"]],axis=1)

        columns["down_type_quarks"]      = ak.concatenate([columns["down"],columns["strange"],columns["bottom"]],axis=1)
        columns["anti_down_type_quarks"] = ak.concatenate([columns["anti_down"],columns["anti_strange"],columns["anti_bottom"]],axis=1)

        columns["positive_leptons"]      = ak.concatenate([columns["anti_electron"],columns["anti_muon"]],axis=1)
        columns["negative_leptons"]      = ak.concatenate([columns["electron"],columns["muon"]],axis=1)
        
        return ak.zip(columns, depth_limit=1, with_name="Event")
