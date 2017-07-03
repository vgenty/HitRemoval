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
my_proc.set_ana_output_file("ana.root");

# Specify data output root file name
my_proc.set_output_file("simplecluster.root")

# prepare the various hit removal stages

algo = fmwk.SimpleClusterer()
algo.setHitProducer("gaushit")
algo.setVtxProducer("mcvertex")
algo.setOutClusProducer("scall")
algo.setRadius(0.4)
algo.setCellSize(1.0)
algo.setUseVertex(True)
algo.setVtxRadius(2.0)
algo.setVerbose(False)
algo.setIgnoreGoF(True)
algo.setROIRadius(200.)

my_proc.add_process( algo )

#my_proc.set_data_to_write(fmwk.data.kHit,         "gaushit" )
#my_proc.set_data_to_write(fmwk.data.kCluster,     "sc2"  )
#my_proc.set_data_to_write(fmwk.data.kAssociation, "sc2"  )

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run(0,10)

sys.exit()

