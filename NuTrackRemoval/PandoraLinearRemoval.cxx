#ifndef LARLITE_PANDORALINEARREMOVAL_CXX
#define LARLITE_PANDORALINEARREMOVAL_CXX

#include "PandoraLinearRemoval.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"
#include "DataFormat/vertex.h"

namespace larlite {

  PandoraLinearRemoval::PandoraLinearRemoval()
    : _tree(nullptr)
  {

    _name        = "PandoraLinearRemoval";
    _fout        = 0;
    
    _verbose     = false;
    _debug       = false;
    
    _clusterProducer = "";
    _vertexProducer  = "";
    
    _max_lin_v = {0.0};
    _min_n_hits_v = {0};

    _vtx_w_cm = {0,0,0};
    _vtx_t_cm = {0,0,0};

  }

  bool PandoraLinearRemoval::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("linearclusterremoval","Linear Cluster Removal TTree");
    _tree->Branch("_nhits",&_nhits,"nhits/I");
    _tree->Branch("_lin",  &_lin  ,"lin/D"  );
    _tree->Branch("_local_lin_truncated",  &_local_lin_truncated  ,"local_lin_truncated/D"  );
    _tree->Branch("_local_lin_avg",  &_local_lin_avg  ,"local_lin_avg/D"  );

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();

    std::cout << "********************************" << std::endl;
    std::cout << "Wire -> cm conversion : " << _wire2cm << std::endl;
    std::cout << "Time -> cm conversion : " << _time2cm << std::endl;
    std::cout << "********************************" << std::endl;

    // make sure _max_lin_v and _min_n_his_v have the same size
    if (_max_lin_v.size() != _min_n_hits_v.size()){
      std::cout << "_max_lin_v and _min_n_hits_v do not have the same size! exit..." << std::endl;
      return false;
    }

    // make sure values are increasing for nhits requirement and decreasing for linearity requirement
    for (size_t i=0; i < _max_lin_v.size() - 1; i++){
      if (_max_lin_v[i+1] < _max_lin_v[i]){
	std::cout << "_max_lin_v values decreasing! quit..." << std::endl;
	return false;
      }
    }
    
    for (size_t i=0; i < _min_n_hits_v.size() - 1; i++){
      if (_min_n_hits_v[i+1] < _min_n_hits_v[i]){
	std::cout << "_min_n_hits_v values decreasing! quit..." << std::endl;
	return false;
      }
    }

    return true;
  }
  
  bool PandoraLinearRemoval::analyze(storage_manager* storage) {
    
    if ( (_clusterProducer == "") || (_vertexProducer == "") )
      { std::cout << "Producers not specified...exit" << std::endl; return false; }

    auto ev_clus = storage->get_data<event_cluster>(_clusterProducer);

    // load vertex
    auto ev_vtx  = storage->get_data<event_vertex>(_vertexProducer);

    larlite::event_hit* ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());

    //set event ID through storage manager
    storage->set_id(storage->run_id(),storage->subrun_id(),storage->event_id());

    if (!ev_hit){
      std::cout << "No hits!" << std::endl;
      return false;
    }

    if (!ev_vtx){
      std::cout << "No hits!" << std::endl;
      return false;
    }

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

    // loop trhough each cluster and calculate linaerity
    // if above some thresdhold, remove cluster

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
      
      // get coordinates of hits to calculate linearity
      std::vector<double> hit_w_v;
      std::vector<double> hit_t_v;

      auto pl = ev_hit->at(hit_idx_v[0]).WireID().Plane;
      
      for (auto const& hit_idx : hit_idx_v){
	hit_w_v.push_back( ev_hit->at(hit_idx).WireID().Wire  * _wire2cm );
	hit_t_v.push_back( ev_hit->at(hit_idx).PeakTime()     * _time2cm );
	//std::cout << "\t\t pt @ [ " << hit_w_v.at(hit_w_v.size()-1) << ", "
	//<< hit_t_v.at(hit_t_v.size()-1) << " ]" << std::endl;
      }

      // figure out the minimum distance between this cluster
      // and the neutrino vertex
      double dvtx = 1e6;
      for (size_t u=0; u < hit_w_v.size(); u++) {
	auto wpt = hit_w_v[u];
	auto tpt = hit_t_v[u];
	double dd = ( ( (wpt - _vtx_w_cm[pl]) * (wpt - _vtx_w_cm[pl]) ) +
		      ( (tpt - _vtx_t_cm[pl]) * (tpt - _vtx_t_cm[pl]) ) );
	if (dd < dvtx) { dvtx = dd; }
      }

      std::cout << "dvtx : " << sqrt(dvtx) << " @ plane " << pl << std::endl;

      if (sqrt(dvtx) > _dvtx_max) continue;
      

      twodimtools::Linearity lin(hit_w_v,hit_t_v);

      _nhits = hit_w_v.size();
      _lin                 = lin._lin;
      _local_lin_avg       = lin._local_lin_avg;
      _local_lin_truncated = lin._local_lin_truncated_avg;

      if (_debug)
	std::cout << "Cluster size : " << hit_w_v.size() << std::endl
		  << "\t lin             : " << lin._lin << std::endl
		  << "\t local lin avg   : " << lin._local_lin_avg << std::endl
		  << "\t local lin trunc : " << lin._local_lin_truncated_avg << std::endl
		  << "\t MAX LIN         : " << max_lin << std::endl;
	
	
      _tree->Fill();

      if (_local_lin_truncated < max_lin){
	if (_debug) std::cout << "\t REMOVE CLUSTER" << std::endl;
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

  bool PandoraLinearRemoval::finalize() {

    if (_fout){
      _fout->cd();
      _tree->Write();
    }
      

    return true;
  }

}
#endif
