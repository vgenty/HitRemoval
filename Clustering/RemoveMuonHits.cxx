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
   e=0;

    return true;
  }
  
  bool RemoveMuonHits::analyze(storage_manager* storage) {

    //std::cout<<"\nNew event! "<<_event<<std::endl ;
    _event++;

    auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex");
    if ( !ev_vtx || !ev_vtx->size() ){ std::cout<<"No Vertex! "<<std::endl ; return false; }

    auto ev_trk = storage->get_data<event_track>("pandoraNu");
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

	  //std::cout<<"STUFF : "<<len<<", "<<dist_st<<", "<<dist_end<<std::endl ;
           }
      }
    // Want to ignore longest track associated to vertex-- Get handle to association 
    auto ev_hit_cosRem = storage->get_data<event_hit>("pandoraCosmicHitRemoval"); //gaushit");

    if ( !ev_hit_cosRem || ev_hit_cosRem->size() == 0 ) {
      std::cout << "No such hits associated to track! " << std::endl;
      return false;
      }

    auto ev_ass = storage->get_data<larlite::event_ass>("pandoraNu"); //pandoraCosmic");

    if ( !ev_ass || ev_ass->size() == 0 ) {
      std::cout << "No such association! " << std::endl;
      return false;
      }

    e=0;
    if( trk_map.size() ) {

      //std::cout<<"Longest track + track ID : "<<1./trk_map.begin()->first<<", "<<trk_map.begin()->second<<std::endl ;

      // Get association to trk => hit and hit => trk
      auto const& ass_hit_v = ev_ass->association(ev_trk->id(), ev_hit_cosRem->id());

      //std::cout<<"Number of Tracks! "<<ass_hit_v.size()<<std::endl;

      if ( ass_hit_v.size() == 0) {
        std::cout << "No ass from track => hit! " << std::endl;
        return false;
        }

      //std::cout<<"HITS ASSOIC WITH TRACK! "<< ass_hit_v.at(trk_map.begin()->second).size()<<std::endl ;

      // Count shower hits in radius
      for(int i = 0; i < ass_hit_v.at(trk_map.begin()->second).size(); i++){

        auto h_i = ass_hit_v.at(trk_map.begin()->second).at(i) ;
	//std::cout<<"GOOD: "<<ev_hit_cosRem->at(h_i).GoodnessOfFit()<<std::endl ;
        //if ( ev_hit_cosRem->at(h_i).WireID().Plane == 2 && 
	if ( ev_hit_cosRem->at(h_i).GoodnessOfFit() != -1 ){
          e++;
          ev_hit_cosRem->at(h_i).set_goodness(-1);
	  }
   
	  //std::cout<<"ID: "<<h_i<<", "<<ev_hit_cosRem->at(h_i).WireID().Plane<<std::endl ;
	}
        //std::cout<<"REMOVING HIT! "<<_event<<std::endl ;

       //if ( e != 0 ) std::cout<<"\nEvent: "<<_event-1 << " removed: "<<e<<" hits"<<std::endl; 

     }



    return true;
  }

  bool RemoveMuonHits::finalize() {

    //if(_fout) { _fout->cd(); _tree->Write(); }

  
    return true;
  }


}
#endif
