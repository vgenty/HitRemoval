#ifndef LARLITE_DRREMOVAL_CXX
#define LARLITE_DRREMOVAL_CXX

#include "DRRemoval.h"
#include "GeoAlgo/GeoAlgo.h"

#include "DataFormat/track.h"

namespace larlite {

  DRRemoval::DRRemoval()
    : HitRemovalBase()
  {

    _name        = "DRRemoval";
    _trackProducer  = "";
    _clusProducer   = "";
    _vertexProducer = "";

    _dmin = _dmax = 0.;
    
  }

  bool DRRemoval::initialize() {

    return true;
  }

  /*
  // This module projects reconstructed pandoraCosmic 3D tracks onto each plane
  // clusters in the event are 
   */
  
  bool DRRemoval::analyze(storage_manager* storage) {

    _event_watch.Start();

    _trk2D_v.clear();
    _trk2D_v = std::vector< std::vector< Traj2D > >(3,std::vector<Traj2D>());

    // geometry helper tools
    auto geoHelper = larutil::GeometryHelper::GetME();

    auto ev_clus = storage->get_data<event_cluster>(_clusProducer   );
    auto ev_trk  = storage->get_data<event_track>  (_trackProducer  );
    auto ev_vtx  = storage->get_data<event_vertex> (_vertexProducer );

    auto ev_clus_out = storage->get_data<event_cluster>("EM");
    auto cluster_ass_v = storage->get_data<event_ass>(ev_clus_out->name());
    std::vector<std::vector<unsigned int> > ass_cluster_hit_v_out;

    //set event ID through storage manager
    storage->set_id(storage->run_id(),storage->subrun_id(),storage->event_id());

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

    // define AABox which encompasses the vertex up to ROI size
    auto const& vtx = ev_vtx->at(0);
    ::geoalgo::AABox ROI(vtx.X()-_roi,vtx.Y()-_roi,vtx.Z()-_roi,
			 vtx.X()+_roi,vtx.Y()+_roi,vtx.Z()+_roi);

    if (_verbose)
      std::cout << "AABox [ " << ROI.Min()[0] << ", " << ROI.Min()[1] << ", " << ROI.Min()[2] << " ] -> [ "
		<< ROI.Max()[0] << ", " << ROI.Max()[1] << ", " << ROI.Max()[2] << " ]" << std::endl;

    // loop through reconstructed tracks and project each onto all 3 planes
    for (size_t t=0; t < ev_trk->size(); t++) {
      auto const& trk = ev_trk->at(t);
      if (_verbose)
	std::cout << "Trk start : [ " << trk.Vertex().X() << ", " << trk.Vertex().Y() << ", " << trk.Vertex().Z() << " ] -> ["
		  << trk.End().X() << ", " << trk.End().Y() << ", " << trk.End().Z() << " ]" << std::endl;
      // the track must exit the ROI for it to be a candidate track for DeltaRay removal
      ::geoalgo::Trajectory traj;
      auto const& npt = trk.NumberTrajectoryPoints();
      for (size_t n=0; n < npt; n++){
	auto const& pt = trk.LocationAtPoint(n);
	traj.push_back( ::geoalgo::Vector(pt.X(),pt.Y(),pt.Z()) );
      }// for all points in the track
      // does the track intersect the ROI?
      auto const& intersection = _geoAlgo.Intersection(ROI,traj);
      if (intersection.size() > 0){
	if (_verbose) std::cout << "\t INTERSECTIon!" << std::endl;
	// project the track to 2D, for each plane
	for (size_t pl=0; pl < 3; pl++){
	  Traj2D trk2D;
	  for (size_t n=0; n < npt; n++){
	    auto const& point = geoHelper->Point_3Dto2D(trk.LocationAtPoint(n), pl);
	    trk2D.push_back(std::make_pair(point.w,point.t + 800. * _time2cm));
	  }// for all track trajectory points
	  _trk2D_v.at(pl).push_back( trk2D );
	}// for all planes
      }// if track intersects ROI
    }// for all tracks

    if (_verbose){
      for (size_t pl=0; pl < 3; pl++)
	std::cout << "@ plane " << pl << " there are " << _trk2D_v[pl].size() << " tracks passing the ROI" << std::endl;
    }

    // for all clusters, check if compatible with Delta Ray
    for (size_t c=0; c < ass_cluster_hit_v.size(); c++) {
      auto const& hit_idx_v = ass_cluster_hit_v[c];
      if (hit_idx_v.size() == 0) continue;
      // check plane
      auto const& pl = _ev_hit->at(hit_idx_v[0]).WireID().Plane;
      if (IsDeltaRay(hit_idx_v,pl) == false){
	ev_clus_out->emplace_back(ev_clus->at(c));
	ass_cluster_hit_v_out.push_back(hit_idx_v);
      }// add to list of output clusters
    }// for all clusters in the event
    
    cluster_ass_v->set_association(ev_clus_out->id(),product_id(data::kHit,_ev_hit->name()),ass_cluster_hit_v_out);    


    _event_time += _event_watch.RealTime();
    _event_num  += 1;
      
    return true;
  }

  // input vectors are the indices of the hits of the muon and delta-ray respectively
  bool DRRemoval::IsDeltaRay(const std::vector<unsigned int>& hitidx_v, const int& pl) {

    // create bounding polygon for this cluster
    double tmin = 10000.;
    double tmax = 0.;
    double wmin = 10000.;
    double wmax = 0.;
    for (auto const& idx : hitidx_v) {
      auto const& hit = _ev_hit->at(idx);
      auto const& t = hit.PeakTime() * _time2cm;
      auto const& w = hit.WireID().Wire * _wire2cm;
      if (t < tmin) { tmin = t; }
      if (t > tmax) { tmax = t; }
      if (w < wmin) { wmin = w; }
      if (w > wmax) { wmax = w; }
    }// find rectangle which bounds delta-ray

    if (_verbose)
      std::cout << "DR bbox : [ " << wmin << ", " << wmax << " ] -> [ " << tmin << ", " << tmax << " ]" << std::endl;

    // loop through all 2D projected tracks
    for (auto const& trj2D : _trk2D_v[pl] ) {

      // does the track intersect the binding box of the Delta Ray?
      bool intersect = false;
      for (auto const& pt2D : trj2D) {
	if (Contained(pt2D,tmin,tmax,wmin,wmax) == true){
	  intersect = true;
	  break;
	}// if point is contained
      }// for all 2D trajectory points

      // if track intersects cluster bounding box, calculate min/max distances
      if (intersect == false) continue;

      if (_verbose) std::cout << "\t contained! " << std::endl;

      double dmin = 10000.;

      // compare all 2D hits of cluster to all 2D projected points in track
      for (auto const& idx : hitidx_v) {
	auto const& hit = _ev_hit->at(idx);
	auto const& t = hit.PeakTime() * _time2cm;
	auto const& w = hit.WireID().Wire * _wire2cm;
	for (auto const& pt2D : trj2D) {
	  double dd = (pt2D.first-w)*(pt2D.first-w) + (pt2D.second-t)*(pt2D.second-t);
	  if (dd < dmin) dmin = dd;
	}// for all trajectory 2D points
      }// for all 2D cluster hits

      dmin = sqrt(dmin);

      if (_verbose) std::cout << "\t dmin : " << dmin << std::endl;

      // is this a delta-ray?
      if ( (dmin < _dmin) )
	return true;

    }// for all 2D trajectories on this plane

      
    return false;
  }
  
  bool DRRemoval::Contained(const std::pair<double,double>& pt2D,
			    const double& tmin, const double& tmax,
			    const double& wmin, const double& wmax) {

    // is the point conatined within the buffer amount?
    if ( (pt2D.first  > wmin-_dmax) && (pt2D.first  < wmax+_dmax) &&
	 (pt2D.second > tmin-_dmax) && (pt2D.second < tmax+_dmax) )
      return true;

    return false;
  }
    
  }
#endif
