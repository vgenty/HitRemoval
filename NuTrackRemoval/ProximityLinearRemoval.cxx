#ifndef LARLITE_PROXIMITYLINEARREMOVAL_CXX
#define LARLITE_PROXIMITYLINEARREMOVAL_CXX

#include "ProximityLinearRemoval.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"

#include "DataFormat/hit.h"
#include "DataFormat/cluster.h"

namespace larlite {

  bool ProximityLinearRemoval::initialize() {

    _vtx_w_cm = std::vector<double> (3,0.);
    _vtx_t_cm = std::vector<double> (3,0.);

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();

    return true;
  }

  
  bool ProximityLinearRemoval::analyze(storage_manager* storage) {
  
    if ( (_clusterProducer == "") || (_vertexProducer == "") ) {
      print(larlite::msg::kERROR,__FUNCTION__,"did not specify producers");
      return false;
    }
    
    auto ev_clus = storage->get_data<event_cluster>(_clusterProducer);

    // load vertex
    auto ev_vtx  = storage->get_data<event_vertex>(_vertexProducer);

    larlite::event_hit* ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());

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

    std::cout << "scanning " << ass_cluster_hit_v.size() << " clusters" << std::endl;

    for (size_t i=0; i < ass_cluster_hit_v.size(); i++){

      auto hit_idx_v = ass_cluster_hit_v[i];

      bool remove = false;

      // determine the linearity threshold for this cluster
      double max_lin = _max_lin_v[0];
      for (size_t n=0; n < _min_n_hits_v.size(); n++){
	auto const& min_n_hits = _min_n_hits_v[n];
	if ( hit_idx_v.size() > min_n_hits )
	  max_lin = _max_lin_v[n];
      }

      if (max_lin <= 0) continue;
      
      // get coordinates of hits to calculate linearity
      std::vector<double> hit_w_v;
      std::vector<double> hit_t_v;

      auto pl = ev_hit->at(hit_idx_v[0]).WireID().Plane;
      
      for (auto const& hit_idx : hit_idx_v){
	hit_w_v.push_back( ev_hit->at(hit_idx).WireID().Wire  * _wire2cm );
	hit_t_v.push_back( ev_hit->at(hit_idx).PeakTime()     * _time2cm );
      }

      twodimtools::Linearity lin(hit_w_v,hit_t_v);

      std::cout << "nhits : " << hit_idx_v.size() << " \t lin : " << lin._local_lin_truncated_avg << std::endl;

      if (lin._local_lin_truncated_avg < max_lin){
	remove = true;
      }

      if (remove) {
	std::cout << "\t\tremoved!" << std::endl;
	for (auto const& hit_idx : hit_idx_v){
	  ev_hit->at(hit_idx).set_goodness(-1.0);
	}
      }
      
    }// for all clusters


    return true;
  }

  bool ProximityLinearRemoval::finalize() {

    return true;
  }


  bool ProximityLinearRemoval::loadVertex(event_vertex* ev_vtx) {

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
