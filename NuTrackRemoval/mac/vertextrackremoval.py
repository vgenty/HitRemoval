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
my_proc.set_output_file("hitremoval_vtx.root")

# prepare the various hit removal stages

algo = fmwk.VertexTrackRemoval()
#algo.setDebug(True)
algo.setClusterProducer("sc")
algo.setVertexProducer("mcvertex")
for n in xrange(15):
    nhits = n * 20
    algo.setMinNHits( nhits )
    print 'nhits : %i',nhits
    if (nhits <= 20) :
        algo.setMaxLinearity(0.03)
        print 'max lin : 0.02'
    else :
        val = (0.3 / 120.) * (nhits - 20)
        algo.setMaxLinearity( val )
        print 'max lin : %.02f'%val
    print
        
algo.setVtxRadius(3.5)
algo.setMaxProtonDist(15.)
algo.setMaxProtonLin(0.2)
algo.setDebug(False)

my_proc.add_process( algo )

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

