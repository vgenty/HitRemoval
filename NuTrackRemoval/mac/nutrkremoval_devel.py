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
my_proc.set_output_file("hitremoval_nutrk_devel.root")

# prepare the various hit removal stages

algo = fmwk.VertexProximityRemoval()
algo.setClusterProducer("pandoraCosmic")
algo.setVertexProducer("randomvertex")
algo.setVtxRadius(3.5)
algo.setMaxLin(0.2)
algo.setVerbose(False)

my_proc.add_process( algo )

# 0th step: removal of pandoraCosmic tracks
    
algo = fmwk.PandoraLinearRemoval()
algo.setClusterProducer("pandoraCosmic")
algo.setVertexProducer("randomvertex")
algo.setSlopeMin(1.0)
algo.setLLT(0.055)
algo.addSlopePt(1,0.05)
algo.addSlopePt(3,0.25)
algo.addSlopePt(10,0.4)
algo.setMaxDVtx(5.)
algo.setROIRadius(100.)
algo.setProtonDMax(10.)
algo.setMaxSSV(0.06)
algo.setVerbose(False)

my_proc.add_process( algo )

# 1st: remove delta-rays

algo = fmwk.TrackDeltaRayRemoval()
algo.setClusterProducer("pandoraCosmic")
algo.setVertexProducer("randomvertex")
algo.setVerbose(False)
algo.setDeltaRayDistMin(1.0);
algo.setDeltaRayDistMax(4.5);
algo.setMaxNHitsDelta(20);
algo.setROI(100.);

my_proc.add_process( algo )

# 2nd step: simple clustering

algo = fmwk.SimpleClusterer()
algo.setHitProducer("gaushit")
algo.setVtxProducer("randomvertex")
algo.setOutClusProducer("sc")
algo.setRadius(0.4)
algo.setCellSize(1.0)
algo.setUseVertex(True)
algo.setVtxRadius(2.0)
algo.setVerbose(False)
algo.setROIRadius(100.)

my_proc.add_process( algo )


algo = fmwk.VertexTrackRemoval()
algo.setClusterProducer("sc")
algo.setVertexProducer("randomvertex")
for n in xrange(15):
    nhits = n * 20
    algo.setMinNHits( nhits )
    print 'nhits : %i',nhits
    if (nhits <= 20) :
        algo.setMaxLinearity(0.03)
    else :
        val = (0.3 / 120.) * (nhits - 20)
        algo.setMaxLinearity( val )
algo.setVtxRadius(3.5)
algo.setMaxProtonDist(15.)
algo.setMaxProtonLin(0.2)
algo.setVerbose(False)

my_proc.add_process( algo )

# 4th: linear cluster removal
algo = fmwk.LinearRemoval()
algo.setClusterProducer("sc")
n_hits = [i*10 for i in xrange(1,30)]
for i,n in enumerate(n_hits):
    algo.setMinNHits( n )
    if (n < 20) :
        algo.setMaxLinearity(0.0001)
    else:
        algo.setMaxLinearity( (0.1 / 120.) * (n - 20) )
algo.setMaxSSV(0.1)
algo.setVerbose(False)

my_proc.add_process( algo )

# 5th: remove hits near vertex

algo = fmwk.RemoveHitsNearVtx()
algo.setHitProducer("gaushit")
algo.setVertexProducer("randomvertex")
algo.setVtxRadius(3.5)

my_proc.add_process( algo )

my_proc.set_data_to_write(fmwk.data.kMCTruth,     "generator"     )
my_proc.set_data_to_write(fmwk.data.kMCShower,    "mcreco"        )
my_proc.set_data_to_write(fmwk.data.kMCTrack,     "mcreco"        )
my_proc.set_data_to_write(fmwk.data.kVertex,      "randomvertex"      )
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

