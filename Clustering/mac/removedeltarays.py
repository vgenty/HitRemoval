import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)


from larlite import larlite as fmwk

from ROOT import twodimtools
a = twodimtools.Linearity()
print a

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    fname = sys.argv[x+1]
    my_proc.add_input_file(fname)
    
# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

# Specify analysis output root file name
my_proc.set_ana_output_file("hitremoval_ana.root");

# Specify data output root file name
my_proc.set_output_file("hitremoval_dr.root")

# prepare the various hit removal stages

_hitproducer    = "gaushit"
_clusproducer   = "pandoraCosmic"
_vtxproducer    = "numuCC_vertex"  #"mcvertex"#"numuCC_vertex"
_clusproducer01 = "clus01"
_clusproducer02 = "clus02"
_clusproducer03 = "clus03"

#muon removal
_trkproducer    = "pandoraNu"
_hitassproducer = "pandoraCosmicKHitRemoval"

def hitremoval00():

    algo = fmwk.RemoveDeltaRays()
    algo.setClusterProducer(_clusproducer)
    algo.setVertexProducer(_vtxproducer)
    algo.setVerbose(False)
    algo.setDeltaRayDistMin(1.0);
    algo.setDeltaRayDistMax(4.5);

    return algo

my_proc.add_process( hitremoval00() )

my_proc.set_data_to_write(fmwk.data.kMCShower, "mcreco"     )
#my_proc.set_data_to_write(fmwk.data.kMCTrack,  "mcreco"     )
my_proc.set_data_to_write(fmwk.data.kMCTruth,  "generator"  )
my_proc.set_data_to_write(fmwk.data.kVertex,   _vtxproducer )

my_proc.set_data_to_write(fmwk.data.kHit,         _hitproducer    )
my_proc.set_data_to_write(fmwk.data.kCluster,     _clusproducer   )
my_proc.set_data_to_write(fmwk.data.kAssociation, _clusproducer   )

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

