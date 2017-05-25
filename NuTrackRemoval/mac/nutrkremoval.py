import sys

from algoconfig import loadAlgo

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
my_proc.set_ana_output_file("ana.root");

# Specify data output root file name
my_proc.set_output_file("hitremoval_nutrk.root")

# prepare the various hit removal stages

my_proc.add_process(loadAlgo("VertexProximityRemoval") )
my_proc.add_process(loadAlgo("PandoraLinearRemoval") )
my_proc.add_process(loadAlgo("TrackDeltaRayRemoval") )
my_proc.add_process(loadAlgo("SimpleClusterer") )
my_proc.add_process(loadAlgo("VertexTrackRemoval") )
my_proc.add_process(loadAlgo("LinearRemoval") )
my_proc.add_process(loadAlgo("RemoveHitsNearVtx") )

my_proc.set_data_to_write(fmwk.data.kMCTruth,     "generator"     )
my_proc.set_data_to_write(fmwk.data.kMCShower,    "mcreco"        )
my_proc.set_data_to_write(fmwk.data.kMCTrack,     "mcreco"        )
my_proc.set_data_to_write(fmwk.data.kVertex,      "mcvertex"  )
my_proc.set_data_to_write(fmwk.data.kCluster,     "sc"            )
my_proc.set_data_to_write(fmwk.data.kAssociation, "sc"            )
my_proc.set_data_to_write(fmwk.data.kCluster,     "pandoraCosmic" )
my_proc.set_data_to_write(fmwk.data.kAssociation, "pandoraCosmic" )
my_proc.set_data_to_write(fmwk.data.kHit,         "gaushit"       )


print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

