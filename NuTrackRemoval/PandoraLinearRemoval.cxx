#ifndef LARLITE_PANDORALINEARREMOVAL_CXX
#define LARLITE_PANDORALINEARREMOVAL_CXX

#include "PandoraLinearRemoval.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"

#include "DataFormat/cluster.h"
#include "DataFormat/vertex.h"

namespace larlite {

  PandoraLinearRemoval::PandoraLinearRemoval()  {

    _name        = "PandoraLinearRemoval";
    _fout        = 0;
    
    _verbose     = false;
    
    _clusterProducer = "";
    _vertexProducer  = "";
    
  }

  bool PandoraLinearRemoval::initialize() {

    std::cout << "********************************" << std::endl;
    std::cout << "Wire -> cm conversion : " << _wire2cm << std::endl;
    std::cout << "Time -> cm conversion : " << _time2cm << std::endl;
    std::cout << "********************************" << std::endl;


    _intercept.clear();
    _slope.clear();

    for (size_t i=1; i < _pts_x_v.size(); i++) {

      auto const& x1 = _pts_x_v[i-1];
      auto const& x2 = _pts_x_v[i];
      auto const& y1 = _pts_y_v[i-1];
      auto const& y2 = _pts_y_v[i];
      
      double s = (y2-y1)/(x2-x1);
      double m = y1 - s * x1;

      if (_verbose) {
	std::cout << "pt0 -> [ " << x1 << ", " << y1 << " ]" << std::endl;
	std::cout << "pt1 -> [ " << x2 << ", " << y2 << " ]" << std::endl;
	std::cout << "slope = " << s << "\t intercept = " << m << std::endl;
      }

      _slope.push_back( s );
      _intercept.push_back( m );
      
    }

    return true;
  }
  
  bool PandoraLinearRemoval::analyze(storage_manager* storage) {
    
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

    for (size_t i=0; i < ass_cluster_hit_v.size(); i++){

      auto hit_idx_v = ass_cluster_hit_v[i];

      bool remove = false;

      // get coordinates of hits to calculate linearity
      std::vector<double> hit_w_v;
      std::vector<double> hit_t_v;

      // don't look at clusters that have already been removed
      if (ev_hit->at(hit_idx_v[0]).GoodnessOfFit() < 0) continue;

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

      // only look at clusters starting close to the vertex.
      if (dvtx_min > _dvtx_max) continue;

      // if the cluster exits the ROI, remove
      if (dvtx_max > _roi_rad) remove = true;

      twodimtools::Linearity lin(hit_w_v,hit_t_v);

      _nhits = hit_w_v.size();
      _lin                 = lin._lin;
      _local_lin_avg       = lin._local_lin_avg;
      _local_lin_truncated = lin._local_lin_truncated_avg;

      if (_verbose)
	std::cout << "Cluster size : " << hit_w_v.size()   << std::endl
		  << "\t slope           : " << lin._slope << std::endl
		  << "\t lin             : " << lin._lin   << std::endl
		  << "\t local lin avg   : " << lin._local_lin_avg           << std::endl
		  << "\t local lin trunc : " << lin._local_lin_truncated_avg << std::endl;

      // remove muons

      // 2 separate cases depending on slope of the cluster

      if (_nhits > 20) {
	
	if ( fabs(lin._slope) < _slope_min) {
	  if (lin._local_lin_truncated_avg < _llt_min) {
	    if (_verbose) std::cout << "\t REMOVE MUON/PION w/ small slope" << std::endl;
	    remove = true;
	  }
	}// small slope case
	else {
	  if (lineCut(fabs(lin._slope),lin._local_lin_truncated_avg) == true) {
	    if (_verbose) std::cout << "\t REMOVE MUON/PION w/ large slope" << std::endl;
	    remove = true;
	  }
	}// if slope is large
	
      }// require min num of hits
      
      // remove protons
      // require maximum SSV
      // and impose a maximum distance from vertex
      if ( (lin._summed_square_variance < _ssv_max) && (dvtx_max < _proton_dmax) ) {
	if (_verbose) std::cout << "\t REMOVE PROTON" << std::endl;
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

  bool PandoraLinearRemoval::finalize() {

    return true;
  }


  bool PandoraLinearRemoval::lineCut(const double& x, const double& y) {

    if (_verbose) { std::cout << "\t\t slope : " << x << "\t lin : " << y << std::endl; } 
    
    for (size_t pt=1; pt < _pts_x_v.size(); pt++) {

      auto const& ptx = _pts_x_v[pt];

      if (x < ptx) {
      
	auto const& s = _slope[pt-1];
	auto const& m = _intercept[pt-1];

	if (_verbose) { std::cout << "\t\t slope < " << ptx << std::endl; }

	if (_verbose) { std::cout << "\t\t cut value is : " << (s*x+m) << std::endl; } 

	if (y < (s*x+m)) return true;
	else return false;
	
      }
      
    }// for all points on the line

    std::cout << "ERROR" << std::endl;
    
    return false;
    
  }

}
#endif
