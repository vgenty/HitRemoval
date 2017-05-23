#ifndef LARLITE_REMOVEHITSNEARVERTEX_CXX
#define LARLITE_REMOVEHITSNEARVERTEX_CXX

#include "RemoveHitsNearVtx.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"

namespace larlite {

  RemoveHitsNearVtx::RemoveHitsNearVtx()
    : HitRemovalBase()
  {

    _name        = "RemoveHitsNearVtx";
    _fout        = 0;
    
  }

  bool RemoveHitsNearVtx::initialize() {

    return true;
  }
  
  bool RemoveHitsNearVtx::analyze(storage_manager* storage) {

    _event_watch.Start();

    if ( (_hitProducer == "") || (_vertexProducer == "") ) {
      print(larlite::msg::kERROR,__FUNCTION__,"did not specify producers");
      return false;
    }
    
    auto ev_hit = storage->get_data<event_hit>(_hitProducer);

    // load vertex
    auto ev_vtx  = storage->get_data<event_vertex>(_vertexProducer);

    //set event ID through storage manager
    storage->set_id(storage->run_id(),storage->subrun_id(),storage->event_id());

    if (!ev_hit){
      print(larlite::msg::kERROR,__FUNCTION__,"no hits");
      return false;
    }

    if (!ev_vtx){
      print(larlite::msg::kERROR,__FUNCTION__,"no vertex");
      return false;
    }
    
    if (loadVertex(ev_vtx) == false) {
      print(larlite::msg::kERROR,__FUNCTION__,"num. vertices != 1");
      return false;
    }

    // loop trhough each cluster and calculate linaerity
    // if above some thresdhold, remove cluster

    for (size_t i=0; i < ev_hit->size(); i++){

      auto hit = ev_hit->at(i);

      if (hit.GoodnessOfFit() < 0) continue;
      
      auto pl = hit.WireID().Plane;

      auto hw = hit.WireID().Wire * _wire2cm;
      auto ht = hit.PeakTime() * _time2cm;
      
      double d = sqrt ( ( ( (hw - _vtx_w_cm[pl]) * (hw - _vtx_w_cm[pl]) ) +
			  ( (ht - _vtx_t_cm[pl]) * (ht - _vtx_t_cm[pl]) ) ) );

      if (d < _vtx_rad)
	ev_hit->at(i).set_goodness(-1.0);
      
    }// for all hits

    _event_time += _event_watch.RealTime();
    _event_num  += 1;
  
  return true;
  }

}
#endif
