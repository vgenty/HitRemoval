#ifndef LARLITE_ROIREMOVAL_CXX
#define LARLITE_ROIREMOVAL_CXX

#include "ROIRemoval.h"

namespace larlite {

  ROIRemoval::ROIRemoval() {

    _name = "ROIRemoval";
    _clusProducer = "";
    _vertexProducer = "";

  }

  bool ROIRemoval::initialize() {

    return true;
  }
  
  bool ROIRemoval::analyze(storage_manager* storage) {

    _event_watch.Start();

    auto ev_clus = storage->get_data<event_cluster>(_clusProducer);
    auto ev_vtx  = storage->get_data<event_vertex> (_vertexProducer);

    larlite::event_hit *ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());

    if (!ev_hit){
      print(larlite::msg::kWARNING,__FUNCTION__,"no hits");
      return false;
    }

    if (loadVertex(ev_vtx) == false) {
      print(larlite::msg::kERROR,__FUNCTION__,"num. vertices != 1");
      return false;
    }

    // loop through all clusters
    //if cluster has hits out of ROI -> remove cluster

    if (_verbose) std::cout << "looping through " << ass_cluster_hit_v.size() << " clusters" << std::endl;
    
    for (size_t i=0; i < ass_cluster_hit_v.size(); i++) {

      bool OutOfROI = false;
      
      auto hit_idx_v = ass_cluster_hit_v[i];

      if (_verbose) std::cout << "new cluster of size " << hit_idx_v.size() << std::endl;

      if (hit_idx_v.size() == 0) continue;

      if (hit_idx_v[0] >= ev_hit->size() )  {
	print(larlite::msg::kERROR,__FUNCTION__,"cluster -> hit ass vector points to hit index out of bounds! skip cluster...");
	continue;
      }

      int pl = ev_hit->at(hit_idx_v[0]).WireID().Plane;

      for (auto const& hit_idx : hit_idx_v) {

	auto const& hit = ev_hit->at(hit_idx);

	double wcm = fabs( hit.WireID().Wire * _wire2cm - _vtx_w_cm[pl] );
	double tcm = fabs( hit.PeakTime() * _time2cm - _vtx_t_cm[pl] );

	if (_verbose) std::cout << "wcm = " << wcm << "\t tcm = " << tcm << std::endl;

	if ( (wcm > _roi) || (tcm > _roi) ) {
	  OutOfROI = true;
	  if (_verbose) { std::cout << "ROI dist = " << wcm << ", " << tcm << " -> remove" << std::endl; }
	  break;
	}
	
      }// for all hits in cluster
      
      if (OutOfROI == true) {
	for (auto const& hit_idx : hit_idx_v)
	  ev_hit->at(hit_idx).set_goodness(-1.0);
      }// if out of ROI

    }// for all clusters

    _event_time += _event_watch.RealTime();
    _event_num  += 1;
      
    return true;
  }

}
#endif
