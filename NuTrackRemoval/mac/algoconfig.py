import sys,os

import ROOT

from larlite import larlite as fmwk

VTXSMEAR = 0.0

def VPR():

    algo = fmwk.VertexProximityRemoval()
    algo.setClusterProducer("pandoraCosmic")
    algo.setVertexProducer("mcvertex")
    algo.setVtxRadius(3.5)
    algo.setMaxLin(0.2)
    algo.setVerbose(False)
    algo.setVertexSmearing(VTXSMEAR)

    return algo
    
# 0th step: removal of pandoraCosmic tracks
def PLR():
    
    algo = fmwk.PandoraLinearRemoval()
    algo.setClusterProducer("pandoraCosmic")
    algo.setVertexProducer("mcvertex")
    algo.setSlopeMin(1.0)
    algo.setLLT(0.055)
    algo.addSlopePt(1,0.05)
    algo.addSlopePt(3,0.25)
    algo.addSlopePt(10,0.4)
    algo.setROIRadius(100.)
    algo.setProtonDMax(10.)
    algo.setMaxSSV(0.06)
    algo.setVertexSmearing(VTXSMEAR)
    algo.setVerbose(False)

    return algo

# 1st: remove delta-rays
def TDR():
    
    algo = fmwk.TrackDeltaRayRemoval()
    algo.setClusterProducer("pandoraCosmic")
    algo.setVertexProducer("mcvertex")
    algo.setVerbose(False)
    algo.setDeltaRayDistMin(1.0);
    algo.setDeltaRayDistMax(4.5);
    algo.setMaxNHitsDelta(20);
    algo.setROI(100.);
    algo.setVertexSmearing(VTXSMEAR)
    
    return algo

# 2nd step: simple clustering
def SC():
    
    algo = fmwk.SimpleClusterer()
    algo.setHitProducer("gaushit")
    algo.setVtxProducer("mcvertex")
    algo.setOutClusProducer("sc")
    algo.setRadius(0.4)
    algo.setCellSize(1.0)
    algo.setUseVertex(True)
    algo.setVtxRadius(2.0)
    algo.setVerbose(False)
    algo.setROIRadius(100.)
    #algo.setVertexSmearing(VTXSMEAR)
    
    return algo

# Vertex Track Removal
def VTR():
    algo = fmwk.VertexTrackRemoval()
    algo.setClusterProducer("sc")
    algo.setVertexProducer("mcvertex")
    for n in xrange(15):
        nhits = n * 20
        algo.setMinNHits( nhits )
        if (nhits <= 20) :
            algo.setMaxLinearity(0.03)
        else :
            val = (0.3 / 120.) * (nhits - 20)
            algo.setMaxLinearity( val )
    algo.setVtxRadius(3.5)
    algo.setMaxProtonDist(15.)
    algo.setMaxProtonLin(0.2)
    algo.setVerbose(False)
    algo.setVertexSmearing(VTXSMEAR)
    
    return algo

# 4th: linear cluster removal
def LIN():
    
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
    algo.setVertexSmearing(VTXSMEAR)
    
    return algo

# 5th: remove hits near vertex
def RHNV():
    
    algo = fmwk.RemoveHitsNearVtx()
    algo.setHitProducer("gaushit")
    algo.setVertexProducer("mcvertex")
    algo.setVtxRadius(3.5)
    algo.setVerbose(False)
    algo.setVertexSmearing(VTXSMEAR)
    
    return algo




def loadAlgo(algoname):

    if (algoname == "VertexProximityRemoval"):
        return VPR()
    if (algoname == "PandoraLinearRemoval"):
        return PLR()
    if (algoname == "TrackDeltaRayRemoval"):
        return TDR()
    if (algoname == "SimpleClusterer"):
        return SC()
    if (algoname == "VertexTrackRemoval"):
        return VTR()
    if (algoname == "LinearRemoval"):
        return LIN()
    if (algoname == "RemoveHitsNearVtx"):
        return RHNV()

    print
    print 'ALGORITHM %s NOT FOUND. QUIT'%algoname
    print

    sys.exit()
