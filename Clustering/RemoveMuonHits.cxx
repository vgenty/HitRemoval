#ifndef LARLITE_REMOVEMUONHITS_CXX
#define LARLITE_REMOVEMUONHITS_CXX

#include "RemoveMuonHits.h"
#include "DataFormat/hit.h"
#include "DataFormat/vertex.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/track.h"
#include <math.h>

namespace larlite {

  bool RemoveMuonHits::initialize() {

   _event = 0 ;
   _vtx_producer="";
   _trk_producer="";
   _hit_ass_producer="";

    return true;
  }
  
  bool RemoveMuonHits::analyze(storage_manager* storage) {

    //std::cout<<"\nNew event! "<<_event<<std::endl ;
    _event++;

    auto ev_vtx = storage->get_data<event_vertex>(_vtx_producer);
    if ( !ev_vtx || !ev_vtx->size() ){ std::cout<<"No Vertex! "<<std::endl ; return false; }

    auto ev_trk = storage->get_data<event_track>(_trk_producer);
    if ( !ev_trk || !ev_trk->size() ) {std::cout<<"No Track!" <<std::endl ; return false; }

    auto vtx = ev_vtx->at(0); 
    std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };

    std::vector<int> trk_ids;
    trk_ids.reserve(ev_trk->size());
    
    //Map of lengths -> track id
    std::multimap<float,int> trk_map ;

    // Find closest + longest pandoraNu track to vertex
    for ( size_t ti = 0; ti < ev_trk->size(); ti++ ) { 

      auto t_vtx = ev_trk->at(ti).Vertex() ;
      auto t_end = ev_trk->at(ti).End() ;
       
      float dist_st = sqrt( pow(t_vtx.X() - vtxXYZ[0],2) + 
                            pow(t_vtx.Y() - vtxXYZ[1],2) + 
                            pow(t_vtx.Z() - vtxXYZ[2],2) ); 

      float dist_end = sqrt( pow(t_end.X() - vtxXYZ[0],2) + 
                             pow(t_end.Y() - vtxXYZ[1],2) + 
                             pow(t_end.Z() - vtxXYZ[2],2) ); 
       if ( dist_st < 2 || dist_end < 2 ){
          float len = ev_trk->at(ti).Length();
          trk_ids.emplace_back(ti);
          trk_map.emplace(1./len,ti);
           }
      }
    //std::cout<<"Longest track + track ID : "<<1./trk_map.begin()->first<<", "<<trk_map.begin()->second<<std::endl ;

    // Want to ignore longest track associated to vertex-- Get handle to association 
    auto ev_hit_cosRem = storage->get_data<event_hit>(_hit_ass_producer);
    auto ev_ass = storage->get_data<larlite::event_ass>(_trk_producer);

    if ( !ev_hit_cosRem || ev_hit_cosRem->size() == 0 ) {
      std::cout << "No such hits associated to track! " << std::endl;
      return false;
      }

    if ( !ev_ass || ev_ass->size() == 0 ) {
      std::cout << "No such association! " << std::endl;
      return false;
      }

    // Get association to trk => hit and hit => trk
    auto const& ass_hit_v = ev_ass->association(ev_trk->id(), ev_hit_cosRem->id());

    if ( ass_hit_v.size() == 0) {
      std::cout << "No ass from track => hit! " << std::endl;
      return false;
      }

    if( trk_map.size() ) {

    // Count shower hits in radius
    for(int i = 0; i < ass_hit_v.at(trk_map.begin()->second).size(); i++){

      auto h_i = ass_hit_v.at(trk_map.begin()->second).at(i) ;
      auto h_ass = ev_hit_cosRem->at(h_i);

      //if( h_ass.GoodnessOfFit() < 0 || h_ass.WireID().Plane != 2 ) continue;

       h_ass.set_goodness(-1) ;

	}
     }

    return true;
  }

  bool RemoveMuonHits::finalize() {

    //if(_fout) { _fout->cd(); _tree->Write(); }

  
    return true;
  }


}
#endif
