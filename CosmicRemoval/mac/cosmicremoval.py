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
my_proc.set_output_file("cosmicremoval.root")

# prepare the various hit removal stages

# ROI REMOVAL
algo = fmwk.ROIRemoval()
algo.setClusterProducer("pandoraCosmic")
algo.setVertexProducer("mcvertex")
algo.setVerbose(False)
algo.setROI(100.)

my_proc.add_process( algo )

# VTX SLOPE CORRELATION REMOVAL

algo = fmwk.VertexSlopeCorrelation()
algo.setClusterProducer("pandoraCosmic")
algo.setVertexProducer("mcvertex")
algo.setVerbose(False)
algo.setCutFunction(80,-10,30,5)
algo.setMinNHits(10)
algo.setROIRadius(100.)

my_proc.add_process( algo )

# VTX ANGLE CORRELATION REMOVAL

algo = fmwk.VertexAngleCorrelation()
algo.setClusterProducer("pandoraCosmic")
algo.setVertexProducer("mcvertex")
algo.setVerbose(False)
algo.setCutFunction(100,-5,15,5)
algo.setMaxAngle(160.)

my_proc.add_process( algo )

# DELTA RAY REMOVAL

algo = fmwk.RemoveDeltaRays()
algo.setClusterProducer("pandoraCosmic")
algo.setVertexProducer("mcvertex")
algo.setVerbose(False)
algo.setDeltaRayDistMin(1.0);
algo.setDeltaRayDistMax(10.0);
algo.setMaxDeltaHits(50);
algo.setROI(120.)

my_proc.add_process( algo )

my_proc.set_data_to_write(fmwk.data.kMCShower, "mcreco"     )
my_proc.set_data_to_write(fmwk.data.kMCTruth,  "generator"  )
my_proc.set_data_to_write(fmwk.data.kVertex,   "mcvertex" )

my_proc.set_data_to_write(fmwk.data.kHit,         "gaushit"    )
my_proc.set_data_to_write(fmwk.data.kCluster,     "pandoraCosmic"   )
my_proc.set_data_to_write(fmwk.data.kAssociation, "pandoraCosmic"   )

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

