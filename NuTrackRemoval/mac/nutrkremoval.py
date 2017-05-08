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
my_proc.set_output_file("hitremoval_nutrk.root")

# prepare the various hit removal stages

# 0th step: removal of pandoraCosmic tracks
def algo00():
    
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

    return algo

# 1st: remove delta-rays
def algo01():

    algo = fmwk.TrackDeltaRayRemoval()
    algo.setClusterProducer("pandoraCosmic")
    algo.setVerbose(False)
    algo.setDeltaRayDistMin(1.0);
    algo.setDeltaRayDistMax(4.5);
    algo.setMaxNHitsDelta(20);

    return algo

# 2nd step: simple clustering
def algo02():

    algo = fmwk.SimpleClusterer()
    algo.setHitProducer("gaushit")
    algo.setVtxProducer("mcvertex")
    algo.setOutClusProducer("sc")
    algo.setRadius(0.4)
    algo.setCellSize(1.0)
    #algo.setMaxHitRMS(19)
    algo.setUseVertex(True)
    algo.setVtxRadius(2.0)
    algo.setVerbose(False)
    #algo.setMinTick(800)
    #algo.setMaxTick(5445)
    return algo

# 3rd: linear clusters near vtx (simple cluster)
def algo03():

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
        else :
            val = (0.3 / 120.) * (nhits - 20)
            algo.setMaxLinearity( val )
    algo.setVtxRadius(3.5)
    algo.setMaxProtonDist(15.)
    algo.setMaxProtonLin(0.2)
    algo.setDebug(False)
    return algo

# 4th: linear cluster removal
def algo04():

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
    algo.setDebug(False)
    return algo

# 5th: remove hits near vertex
def algo05():
    algo = fmwk.RemoveHitsNearVtx()
    algo.setHitProducer("gaushit")
    algo.setVertexProducer("mcvertex")
    algo.setVtxRadius(3.5)

    return algo

my_proc.add_process( algo00() )
my_proc.add_process( algo01() )
my_proc.add_process( algo02() )
my_proc.add_process( algo03() )
my_proc.add_process( algo04() )
my_proc.add_process( algo05() )


#my_proc.set_data_to_write(fmwk.data.kMCTruth,     "generator"     )
#my_proc.set_data_to_write(fmwk.data.kVertex,      "mcvertex"      )
#my_proc.set_data_to_write(fmwk.data.kCluster,     "sc"            )
#my_proc.set_data_to_write(fmwk.data.kAssociation, "sc"            )
#my_proc.set_data_to_write(fmwk.data.kCluster,     "pandoraCosmic" )
#my_proc.set_data_to_write(fmwk.data.kAssociation, "pandoraCosmic" )
#my_proc.set_data_to_write(fmwk.data.kHit,         "gaushit"       )


print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

