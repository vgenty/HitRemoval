#ifndef LARLITE_VERTEXPROXIMITYREMOVAL_CXX
#define LARLITE_VERTEXPROXIMITYREMOVAL_CXX

#include "VertexProximityRemoval.h"

namespace larlite {

  VertexProximityRemoval::VertexProximityRemoval()  {

    _name        = "VertexProximityRemoval";
    _fout        = 0;
    
    _clusterProducer = "";
    _vertexProducer  = "";
    
  }

  bool VertexProximityRemoval::initialize() {

    return true;
  }
  
  bool VertexProximityRemoval::analyze(storage_manager* storage) {

    _event_watch.Start();

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

    if (loadVertex(ev_vtx) == false) {
      print(larlite::msg::kERROR,__FUNCTION__,"num. vertices != 1");
      return false;
    }
    
    // loop trhough each cluster and calculate linaerity
    // if above some thresdhold, remove cluster

    // select only clusters not previously removed
    auto const& clus_idx_v = AvailableClusterIndices(ev_hit, ass_cluster_hit_v);

    std::cout << "examining " << clus_idx_v.size()
	      << " of " << ass_cluster_hit_v.size() << " total clusters" << std::endl; 

    for (auto const& i : clus_idx_v) {

      auto hit_idx_v = ass_cluster_hit_v[i];

      if (hit_idx_v.size() == 0) continue;

      bool remove = false;

      // get coordinates of hits to calculate linearity
      std::vector<double> hit_w_v;
      std::vector<double> hit_t_v;

      auto pl = ev_hit->at(hit_idx_v[0]).WireID().Plane;
      
      for (auto const& hit_idx : hit_idx_v){
	hit_w_v.push_back( ev_hit->at(hit_idx).WireID().Wire  * _wire2cm );
	hit_t_v.push_back( ev_hit->at(hit_idx).PeakTime()     * _time2cm );
      }

      // figure out the minimum distance between this cluster
      // and the neutrino vertex
      double dvtx_min = 1e6;
      double dvtx_max = 0.;

      for (size_t u=0; u < hit_w_v.size(); u++) {

	auto wpt = hit_w_v[u];
	auto tpt = hit_t_v[u];
	
	double dd = ( ( (wpt - _vtx_w_cm[pl]) * (wpt - _vtx_w_cm[pl]) ) +
		      ( (tpt - _vtx_t_cm[pl]) * (tpt - _vtx_t_cm[pl]) ) );
	
	if (dd < dvtx_min) { dvtx_min = dd; }
	if (dd > dvtx_max) { dvtx_max = dd; }
	
      }

      dvtx_min = sqrt(dvtx_min);
      dvtx_max = sqrt(dvtx_max);

      if (_verbose)
	std::cout << "Cluster size       : " << hit_w_v.size() << std::endl
		  << "\t plane           : " << pl << std::endl
		  << "\t vtx dist        : " << dvtx_min << std::endl;


      if (dvtx_min > _vtx_rad) continue;
      
      twodimtools::Linearity lin(hit_w_v,hit_t_v);

      if (lin._local_lin_truncated_avg < _lin_max) remove = true;
      
      if (remove) {
	for (auto const& hit_idx : hit_idx_v)
	  ev_hit->at(hit_idx).set_goodness(-1.0);
      }
      
    }// for all clusters

    _event_time += _event_watch.RealTime();
    _event_num  += 1;
  
  return true;
  }

}
#endif
