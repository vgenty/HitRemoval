#ifndef LARLITE_LINEARREMOVAL_CXX
#define LARLITE_LINEARREMOVAL_CXX

#include "LinearRemoval.h"

namespace larlite {

  bool LinearRemoval::initialize() {

    return true;
  }

  
  bool LinearRemoval::analyze(storage_manager* storage) {
  
    if ( (_clusterProducer == "") ) {
      print(larlite::msg::kERROR,__FUNCTION__,"did not specify producers");
      return false;
    }
    
    auto ev_clus = storage->get_data<event_cluster>(_clusterProducer);

    larlite::event_hit* ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());

    //set event ID through storage manager
    storage->set_id(storage->run_id(),storage->subrun_id(),storage->event_id());

    if (!ev_hit){
      print(larlite::msg::kERROR,__FUNCTION__,"no hits");
      return false;
    }

    // loop trhough each cluster and calculate linaerity
    // if above some thresdhold, remove cluster

    // select only clusters not previously removed
    auto const& clus_idx_v = AvailableClusterIndices(ev_hit, ass_cluster_hit_v);

    for (size_t i=0; i < clus_idx_v.size(); i++){

      auto hit_idx_v = ass_cluster_hit_v[i];

      bool remove = false;

      // determine the linearity threshold for this cluster
      double max_lin = _max_lin_v[0];
      for (size_t n=0; n < _min_n_hits_v.size(); n++){
	auto const& min_n_hits = _min_n_hits_v[n];
	if ( hit_idx_v.size() > min_n_hits )
	  max_lin = _max_lin_v[n];
      }

      if (max_lin < 0) continue;
      if (hit_idx_v.size() < 8) continue;
      
      // get coordinates of hits to calculate linearity
      std::vector<double> hit_w_v;
      std::vector<double> hit_t_v;

      auto pl = ev_hit->at(hit_idx_v[0]).WireID().Plane;
      
      for (auto const& hit_idx : hit_idx_v){
	hit_w_v.push_back( ev_hit->at(hit_idx).WireID().Wire  * _wire2cm );
	hit_t_v.push_back( ev_hit->at(hit_idx).PeakTime()     * _time2cm );
      }

      twodimtools::Linearity lin(hit_w_v,hit_t_v);

      if (lin._local_lin_truncated_avg < max_lin)
	remove = true;

      if (_verbose)
	std::cout << "Pl : " << pl << "\t nhit : " << hit_w_v.size() << "\t lin : " << lin._local_lin_truncated_avg
		  << "\t ssv : " << lin._summed_square_variance << std::endl;
      
      if ( lin._summed_square_variance < _ssv ) {
	if (_verbose) std::cout << "\tremoved SSV" << std::endl;
	remove = true;
      }

      if (remove) {
	for (auto const& hit_idx : hit_idx_v){
	  ev_hit->at(hit_idx).set_goodness(-1.0);
	}
      }
      
    }// for all clusters


    return true;
  }

  bool LinearRemoval::finalize() {

    return true;
  }

}
#endif
