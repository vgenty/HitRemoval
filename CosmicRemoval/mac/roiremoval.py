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
my_proc.set_ana_output_file("roiremoval_ana.root");

# Specify data output root file name
my_proc.set_output_file("hitremoval_roi.root")

# prepare the various hit removal stages

algo = fmwk.ROIRemoval()
algo.setClusterProducer("pandoraCosmic")
algo.setVertexProducer("mcvertex")
algo.setVerbose(True)
algo.setROI(100.)

my_proc.add_process( algo )

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run(0,100)

sys.exit()

