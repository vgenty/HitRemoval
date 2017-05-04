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
my_proc.set_output_file("hitremoval_pandora.root")

# prepare the various hit removal stages

algo = fmwk.PandoraLinearRemoval()
#algo.setDebug(True)
algo.setClusterProducer("pandoraCosmic")
algo.setVertexProducer("mcvertex")
algo.setSlopeMin(1.0)
algo.setLLT(0.055)
algo.addSlopePt(1,0.05)
algo.addSlopePt(3,0.25)
algo.addSlopePt(10,0.4)
algo.setMaxDVtx(5.)
algo.setROIRadius(100.)
algo.setProtonDMax(10.)
algo.setMaxSSV(0.06)
algo.setDebug(False)

my_proc.add_process( algo )

my_proc.set_data_to_write(fmwk.data.kVertex,  "mcvertex"      )
my_proc.set_data_to_write(fmwk.data.kCluster, "pandoraCosmic" )
my_proc.set_data_to_write(fmwk.data.kAssociation, "pandoraCosmic" )
my_proc.set_data_to_write(fmwk.data.kHit,     "gaushit" )

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run(0,100)

sys.exit()

