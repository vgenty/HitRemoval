#ifndef LARLITE_REMOVEDELTARAYS_CXX
#define LARLITE_REMOVEDELTARAYS_CXX

#include "RemoveDeltaRays.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"
#include "DataFormat/vertex.h"

namespace larlite {

  RemoveDeltaRays::RemoveDeltaRays() {

    _name        = "RemoveDeltaRays";
    _fout        = 0;
    _verbose     = false;
    _clusProducer = "";
    _vertexProducer = "";

    _d_delta_min = _d_delta_max = 0.;

    _vtx_w_cm = {0,0,0};
    _vtx_t_cm = {0,0,0};
    
  }

  bool RemoveDeltaRays::initialize() {

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();
    
    return true;
  }
  
  bool RemoveDeltaRays::analyze(storage_manager* storage) {

    auto ev_clus = storage->get_data<event_cluster>(_clusProducer);
    auto ev_vtx  = storage->get_data<event_vertex> (_vertexProducer);

    _ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), _ev_hit, ev_clus->name());

    if (!_ev_hit){
      print(larlite::msg::kWARNING,__FUNCTION__,"no hits");
      return false;
    }

    if (loadVertex(ev_vtx) == false) {
      print(larlite::msg::kERROR,__FUNCTION__,"num. vertices != 1");
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

      // only consider as cosmics the tracks with more hits than the max allowed for a delta-ray
      if (hit_idx_v.size() < _max_delta_hits) continue;

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
      if (hit_idx_v.size() > _max_delta_hits) continue;

      // if out of ROI, ignore this delta-ray (OK not to remove)
      auto const& hit0 = _ev_hit->at(hit_idx_v.at(0));
      double wcm = hit0.WireID().Wire * _wire2cm;
      double tcm = hit0.PeakTime() * _time2cm;
      double dvtx = sqrt( (wcm - _vtx_w_cm[pl]) * (wcm - _vtx_w_cm[pl]) +
			  (tcm - _vtx_t_cm[pl]) * (tcm - _vtx_t_cm[pl]) );
      if (dvtx > _roi) continue;
      
      // compare this delta-ray to all removed muons in the plane
      for (auto const& muidx : cosmic_clus_v[pl]){
	if (DeltaRay( ass_cluster_hit_v[muidx], hit_idx_v ) == true )
	delta_ray_v.push_back( i );
      }
      
    }// for all clusters

    for (auto const& delta_ray : delta_ray_v){
      auto hit_idx_v = ass_cluster_hit_v[delta_ray];
      for (auto const& hit_idx : hit_idx_v){
	//std::cout << "removing hit @ " << _ev_hit->at(hit_idx).PeakTime() * _time2cm << ", " << _ev_hit->at(hit_idx).WireID().Wire * _wire2cm << "]" << std::endl;
	_ev_hit->at(hit_idx).set_goodness(-1.0);
      }
    }// for all delta-rays
      
    return true;
  }

  bool RemoveDeltaRays::finalize() {

    return true;
  }

  
  bool RemoveDeltaRays::DeltaRay(const std::vector<unsigned int>& muon,
				 const std::vector<unsigned int>& deltaray) {

    // min dist to muon
    // idx of muon hit with min dist
    double ddmin = 1000000.;
    size_t hidxmin;
    // max distance to muon
    double ddmax = 0.0;

    //std::cout << muon.size() << " muon hits" << std::endl;
    //std::cout << deltaray.size() << " delta-ray hits" << std::endl;
    //std::cout << "delta-ray @ [" << _ev_hit->at(deltaray[0]).PeakTime() * _time2cm << ", " << _ev_hit->at(deltaray[0]).WireID().Wire * _wire2cm << "]" << std::endl;
    
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

    //std::cout << "dmin = " << (int)(10.*sqrt(ddmin)) << ", dmax = " << (int)(10.*sqrt(ddmax)) << " mm" << std::endl;
	
    return true;
  }


  double RemoveDeltaRays::_distSq_(const larlite::hit& h1, const larlite::hit& h2) {

    double t1 = h1.PeakTime() * _time2cm;
    double w1 = h1.WireID().Wire * _wire2cm;
    double t2 = h2.PeakTime() * _time2cm;
    double w2 = h2.WireID().Wire * _wire2cm;
    
    return (t1-t2)*(t1-t2) + (w1-w2)*(w1-w2);
  }

  
  bool RemoveDeltaRays::loadVertex(event_vertex* ev_vtx) {
    
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
