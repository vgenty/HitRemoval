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
my_proc.set_output_file("hitremoval.root")

# prepare the various hit removal stages

algo = fmwk.PandoraLinearRemoval()
#algo.setDebug(True)
algo.setClusterProducer("pandoraCosmic")
algo.setVertexProducer("mcvertex")
n_hits = [  3, 10, 20, 30, 50, 70,100,130,180,230]
lin_v  = [.02,.10,.10,.13,.16,.19,.22,.30,.40,.45]
for i,n in enumerate(n_hits):
    algo.setMaxLinearity( lin_v[i] )
    algo.setMinNHits( n )

algo.setMaxDVtx(5.)
    
my_proc.add_process( algo )

my_proc.set_data_to_write(fmwk.data.kVertex,  "mcvertex"      )
my_proc.set_data_to_write(fmwk.data.kCluster, "pandoraCosmic" )
my_proc.set_data_to_write(fmwk.data.kAssociation, "pandoraCosmic" )
my_proc.set_data_to_write(fmwk.data.kHit,     "gaushit" )

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run(0,7)

sys.exit()

