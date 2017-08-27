import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)


from larlite import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-2):
    fname = sys.argv[x+1]
    my_proc.add_input_file(fname)
    
# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

# Specify analysis output root file name
my_proc.set_ana_output_file("linearclusterremoval.root");

# Specify data output root file name
my_proc.set_output_file(sys.argv[-1])

algo = fmwk.PhotonClusterer()
algo.setHitProducer("gaushit")
algo.setVtxProducer("numuCC_vertex")
#algo.setVtxProducer("sel2")
algo.setOutClusProducer("photon")
algo.setRadius(1.0)
algo.setCellSize(3.0)
algo.setUseVertex(True)
algo.setVerbose(False)
algo.setROIRadius(200.)
algo.setMaxClusSize(100)
algo.setMinLin(0.06)

my_proc.add_process(algo)

#my_proc.set_data_to_write(fmwk.data.kHit,hitproducer)
my_proc.set_data_to_write(fmwk.data.kCluster,'photon')
my_proc.set_data_to_write(fmwk.data.kAssociation,'photon')
#my_proc.set_data_to_write(fmwk.data.kHit,"shrhits")

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

