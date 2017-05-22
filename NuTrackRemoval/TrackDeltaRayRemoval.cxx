#ifndef LARLITE_TRACKDELTARAYREMOVAL_CXX
#define LARLITE_TRACKDELTARAYREMOVAL_CXX

#include "TrackDeltaRayRemoval.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"
#include "DataFormat/vertex.h"

namespace larlite {

  TrackDeltaRayRemoval::TrackDeltaRayRemoval() {

    _name        = "TrackDeltaRayRemoval";
    _fout        = 0;
    _verbose     = false;
    _clusProducer = "";

    _d_delta_min = _d_delta_max = 0.;
    _nhitmax = 0;
    
  }

  bool TrackDeltaRayRemoval::initialize() {

    if ( _clusProducer == "" ) {
      print(larlite::msg::kERROR,__FUNCTION__,"did not specify producers");
      return false;
    }

    if ( (_nhitmax == 0) || (_d_delta_min == 0) || (_d_delta_max == 0) ) {
      print(larlite::msg::kERROR,__FUNCTION__,"did not set algorithm input parameters");
      return false;
    }
      
    return true;
  }
  
  bool TrackDeltaRayRemoval::analyze(storage_manager* storage) {

    _event_watch.Start();

    auto ev_clus = storage->get_data<event_cluster>(_clusProducer);

    _ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), _ev_hit, ev_clus->name());

    if (!_ev_hit){
      print(larlite::msg::kERROR,__FUNCTION__,"no hits");
      return false;
    }

    // keep track of all clusters identified as delta-rays
    std::vector<size_t> delta_ray_v;
    
    // loop through clusters and identify those removed
    // by cosmic removal, per plane
    std::vector< std::vector<size_t> > cosmic_clus_v(3,std::vector<size_t>());

    for (size_t i=0; i < ass_cluster_hit_v.size(); i++) {

      auto hit_idx_v = ass_cluster_hit_v[i];
      
      int pl = _ev_hit->at(hit_idx_v[0]).WireID().Plane;

      if (hit_idx_v.size() == 0) continue;

      // require min number of hits to be a cosmic
      if (hit_idx_v.size() < _nhitmax) continue;

      // if negative GoF -> removed.
      if ( _ev_hit->at(hit_idx_v.at(0)).GoodnessOfFit() < 0 )
	cosmic_clus_v[pl].push_back( i );

    }// for all clusters

    // now for small clusters that have not been removed
    // identify if they are close to a removed cosmic
    // cluster and delta-ray like

    for (size_t i=0; i < ass_cluster_hit_v.size(); i++) {
      
      auto hit_idx_v = ass_cluster_hit_v[i];

      int pl = _ev_hit->at(hit_idx_v[0]).WireID().Plane;

      if (hit_idx_v.size() == 0) continue;

      // if negative GoF -> removed, ignore
      if ( _ev_hit->at(hit_idx_v.at(0)).GoodnessOfFit() < 0 ) continue;

      // if too many hits -> remove
      if (hit_idx_v.size() > _nhitmax) continue;
      
      // compare this delta-ray to all removed muons in the plane
      for (auto const& muidx : cosmic_clus_v[pl]){
	if (DeltaRay( ass_cluster_hit_v[muidx], hit_idx_v ) == true )
	delta_ray_v.push_back( i );
      }
      
    }// for all clusters

    for (auto const& delta_ray : delta_ray_v){
      auto hit_idx_v = ass_cluster_hit_v[delta_ray];
      for (auto const& hit_idx : hit_idx_v)
	_ev_hit->at(hit_idx).set_goodness(-1.0);
    }// for all delta-rays

    _event_time += _event_watch.RealTime();
    _event_num  += 1;

    return true;
  }

  bool TrackDeltaRayRemoval::DeltaRay(const std::vector<unsigned int>& muon,
				      const std::vector<unsigned int>& deltaray) {

    // min dist to muon
    // idx of muon hit with min dist
    double ddmin = 1000000.;
    size_t hidxmin;
    // max distance to muon
    double ddmax = 0.0;

    // find minimum distance between delta-ray candidate and cosmic muon
    for (auto mu_h_idx : muon) {
      for (auto dr_h_idx : deltaray) {

	auto mu_h = _ev_hit->at(mu_h_idx);
	auto dr_h = _ev_hit->at(dr_h_idx);

	double dd = _distSq_(mu_h,dr_h);

	//std::cout << "\t dist = " << dd << std::endl;

	if (dd < ddmin) { ddmin = dd; hidxmin = mu_h_idx; }
      }
    }

    //std::cout << "ddmin is " << ddmin << std::endl;

    if (  sqrt(ddmin) > _d_delta_min ) return false;

    // find max distance to this point
    auto const& mu_h_min = _ev_hit->at(hidxmin);
    
    for (auto const& dr_h_idx : deltaray) {
      
      auto const& dr_h = _ev_hit->at(dr_h_idx);
      
      double dd = _distSq_(dr_h, mu_h_min);

      if (dd > ddmax) { ddmax = dd; }

    }// for all delta-ray hits

    if ( sqrt(ddmax) > _d_delta_max ) return false;

    return true;
  }


  double TrackDeltaRayRemoval::_distSq_(const larlite::hit& h1, const larlite::hit& h2) {

    double t1 = h1.PeakTime() * _time2cm;
    double w1 = h1.WireID().Wire * _wire2cm;
    double t2 = h2.PeakTime() * _time2cm;
    double w2 = h2.WireID().Wire * _wire2cm;
    
    return (t1-t2)*(t1-t2) + (w1-w2)*(w1-w2);
  }
      

}
#endif
