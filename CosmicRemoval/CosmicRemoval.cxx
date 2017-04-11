#ifndef LARLITE_COSMICREMOVAL_CXX
#define LARLITE_COSMICREMOVAL_CXX

#include "CosmicRemoval.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"
#include "DataFormat/vertex.h"

namespace larlite {

  CosmicRemoval::CosmicRemoval()
    : _tree(nullptr)
  {

    _name        = "CosmicRemoval";
    _fout        = 0;
    _verbose     = false;
    _clusterProducer = "";
    _vtxProducer     = "";
    _out_hitProducer = "";
    _max_lin_v = {0.0};
    _min_n_hits_v = {0};

    _vtx_w_cm = {0,0,0};
    _vtx_t_cm = {0,0,0};

  }

  bool CosmicRemoval::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("linearclusterremoval","Linear Cluster Removal TTree");
    _tree->Branch("_nhits",&_nhits,"nhits/I");
    _tree->Branch("_lin",  &_lin  ,"lin/D"  );
    _tree->Branch("_angle",  &_angle  ,"angle/D"  );
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
  
  bool CosmicRemoval::analyze(storage_manager* storage) {

    auto ev_clus = storage->get_data<event_cluster>( _clusterProducer );
    auto evt_vtx = storage->get_data<event_vertex> ( _vtxProducer     );

    larlite::event_hit* ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());
    
    //auto out_hit = storage->get_data<event_hit>( _out_hitProducer );
    
    //set event ID through storage manager
    storage->set_id(storage->run_id(),storage->subrun_id(),storage->event_id());

    if (!ev_hit){
      std::cout << "No hits!" << std::endl;
      return false;
    }

    // get vertex position on each plane
    if ( (evt_vtx->size() == 1) ){
      auto const& vtx = evt_vtx->at(0);
      auto geoH = larutil::GeometryHelper::GetME();
      std::vector<double> xyz = {vtx.X(), vtx.Y(), vtx.Z()};
      for (size_t pl = 0; pl < 3; pl++){
	auto const& pt = geoH->Point_3Dto2D(xyz,pl);
	_vtx_w_cm[pl] = pt.w;
	_vtx_t_cm[pl] = pt.t + 800 * _time2cm;
      }
    }

    // loop trhough each cluster and calculate linaerity
    // if above some thresdhold, remove cluster

    for (size_t i=0; i < ass_cluster_hit_v.size(); i++){

      // store output cluster hit indices
      std::vector<unsigned int> out_cluster_hit_idx_v;
      larlite::cluster out_clus;

      auto hit_idx_v = ass_cluster_hit_v[i];

      bool remove = false;

      // determine the linearity threshold for this cluster
      double max_lin = _max_lin_v[0];
      for (size_t n=0; n < _min_n_hits_v.size(); n++){
	auto const& min_n_hits = _min_n_hits_v[n];
	if ( hit_idx_v.size() > min_n_hits )
	  max_lin = _max_lin_v[n];
      }

      int pl = ev_hit->at(hit_idx_v[0]).WireID().Plane;

      // get coordinates of hits to calculate linearity
      std::vector<double> hit_w_v;
      std::vector<double> hit_t_v;

      // minimum distance to vtx
      // and coordinates for closest point
      double dmin = 1000.;
      double wmin = 0;
      double tmin = 0;
      
      for (auto const& hit_idx : hit_idx_v){
	double w = ev_hit->at(hit_idx).WireID().Wire  * _wire2cm;
	double t = ev_hit->at(hit_idx).PeakTime() * _time2cm;
	hit_w_v.push_back( w );
	hit_t_v.push_back( t );
	double d = ( w - _vtx_w_cm[pl] ) * ( w - _vtx_w_cm[pl] ) + ( t - _vtx_t_cm[pl] ) * ( t - _vtx_t_cm[pl] );
	if (d < dmin) { dmin = d; wmin = w; tmin = t; }
      }

      dmin = sqrt(dmin);

      if (_verbose)
	std::cout << "Pl : " << pl << "\t NHits : " << hit_w_v.size() << std::endl;
      
      twodimtools::Linearity lin(hit_w_v,hit_t_v);

      _nhits = hit_w_v.size();
      _lin                 = lin._lin;
      _local_lin_avg       = lin._local_lin_avg;
      _local_lin_truncated = lin._local_lin_truncated_avg;

      // get dot product between slope line and mean to vtx line
      double _vtx_mean_w = _vtx_w_cm[pl] - lin._meanx;
      double _vtx_mean_t = _vtx_t_cm[pl] - lin._meany;
      double mag = sqrt( _vtx_mean_w * _vtx_mean_w + _vtx_mean_t * _vtx_mean_t );
      _vtx_mean_w /= mag;
      _vtx_mean_t /= mag;
      
      double slope     = lin._slope;
      double intercept = lin._intercept;

      // impact parameter to vertex:
      double x0 = _vtx_w_cm[pl];
      double y0 = _vtx_t_cm[pl];
      double IP = fabs( - slope * x0 + y0 - intercept ) / sqrt( slope * slope + 1 );

      // slope of vertex to minimum point
      double slopeminpt = (_vtx_t_cm[pl] - tmin) / (_vtx_w_cm[pl] - wmin);
      double tanangle = (slopeminpt - slope) / (1 + slopeminpt * slope);
      _angle = fabs( atan(tanangle) * 180 / 3.14 );

      if (_verbose)
	std::cout << "\t Linearity : " << _local_lin_truncated
		  << "\t IP : " << IP
		  << "\t DMIN : " << dmin
		  << "\t angle : " << _angle
		  << std::endl << std::endl;

      _tree->Fill();

      if (_local_lin_truncated < max_lin and ( _angle > 20. ) )
	remove = true;

      if ( (IP > 20) and (_nhits > 20) and (_local_lin_truncated < 0.3) )
	remove = true;

      if (remove){
	for (auto const& hit_idx : hit_idx_v)
	  ev_hit->at(hit_idx).set_goodness(-1.0);
      }
      
    }// loop through all clusters

    return true;
  }

  bool CosmicRemoval::finalize() {

    if (_fout){
      _fout->cd();
      _tree->Write();
    }
      

    return true;
  }

}
#endif
