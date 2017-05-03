#ifndef LARLITE_REMOVEHITSNEARVERTEX_CXX
#define LARLITE_REMOVEHITSNEARVERTEX_CXX

#include "RemoveHitsNearVtx.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"

#include "DataFormat/cluster.h"
#include "DataFormat/vertex.h"

namespace larlite {

  RemoveHitsNearVtx::RemoveHitsNearVtx()  {

    _name        = "RemoveHitsNearVtx";
    _fout        = 0;
    
    _verbose     = false;
    _debug       = false;
    
    _hitProducer     = "";
    _vertexProducer  = "";
    
    _vtx_w_cm = {0,0,0};
    _vtx_t_cm = {0,0,0};

  }

  bool RemoveHitsNearVtx::initialize() {

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();

    std::cout << "********************************" << std::endl;
    std::cout << "Wire -> cm conversion : " << _wire2cm << std::endl;
    std::cout << "Time -> cm conversion : " << _time2cm << std::endl;
    std::cout << "********************************" << std::endl;

    return true;
  }
  
  bool RemoveHitsNearVtx::analyze(storage_manager* storage) {

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
  
  return true;
  }

  bool RemoveHitsNearVtx::finalize() {

    return true;
  }


  bool RemoveHitsNearVtx::loadVertex(event_vertex* ev_vtx) {

    if (ev_vtx->size() != 1) return false;
    
    // get vertex position on each plane
    if ( (ev_vtx->size() == 1) ){
      auto const& vtx = ev_vtx->at(0);
      auto geoH = larutil::GeometryHelper::GetME();
      auto geom = larutil::Geometry::GetME();
      std::vector<double> xyz = {vtx.X(), vtx.Y(), vtx.Z()};
      for (size_t pl = 0; pl < 3; pl++){
	double *origin;
	origin = new double[3];
	geom->PlaneOriginVtx(pl,origin);
	auto const& pt = geoH->Point_3Dto2D(xyz,pl);
	_vtx_w_cm[pl] = pt.w;
	_vtx_t_cm[pl] = pt.t + 800 * _time2cm - origin[0];
      }
    }    

    return true;
  }

}
#endif
