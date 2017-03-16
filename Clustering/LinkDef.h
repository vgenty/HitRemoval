//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class larlite::SimpleClusterer+;
#pragma link C++ class larlite::ProximityClusterer+;
#pragma link C++ class larlite::LinearHitRemoval+;
#pragma link C++ class larlite::LinearClusterRemoval+;
#pragma link C++ class larlite::CosmicRemoval+;
#pragma link C++ class larlite::LinearClusterLocalRemoval+;
#pragma link C++ class larlite::ClusterFilter+;
#pragma link C++ class larlite::PhotonClusterer+;
#pragma link C++ class larlite::RemoveMuonHits+;
#pragma link C++ class larlite::RemoveCosmicTracks+;
#pragma link C++ class larlite::RemoveDeltaRays+;
//ADD_NEW_CLASS ... do not change this line
#endif








