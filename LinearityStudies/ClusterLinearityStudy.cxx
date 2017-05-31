#ifndef LARLITE_CLUSTERLINEARITYSTUDY_CXX
#define LARLITE_CLUSTERLINEARITYSTUDY_CXX

#include "ClusterLinearityStudy.h"

#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"

namespace larlite {

  ClusterLinearityStudy::ClusterLinearityStudy()
    : _tree(nullptr)
  {
    _name="ClusterLinearityStudy";
    _fout=0;
    _vtx_w_cm = {0.,0.,0.};
    _vtx_t_cm = {0.,0.,0.};

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();
    
  }

  bool ClusterLinearityStudy::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("_tree","tree");
    _tree->Branch("_pl",&_pl,"pl/I");
    _tree->Branch("_dvtx_max" ,&_dvtx_max ,"dvtx_max/D" );
    _tree->Branch("_dvtx_min" ,&_dvtx_min ,"dvtx_min/D" );  
    _tree->Branch("_nhits",&_nhits,"nhits/I");
    _tree->Branch("_lin"  ,&_lin  ,"lin/D"  );
    _tree->Branch("_ssv"  ,&_ssv  ,"ssv/D"  );
    _tree->Branch("_slope",&_slope,"slope/D");
    _tree->Branch("_local_lin_truncated",  &_local_lin_truncated  ,"local_lin_truncated/D"  );
    _tree->Branch("_local_lin_avg",  &_local_lin_avg  ,"local_lin_avg/D"  );

    return true;
  }
  
  bool ClusterLinearityStudy::analyze(storage_manager* storage) {
  
    auto ev_clus = storage->get_data<event_cluster>(_clusterProducer);
    auto ev_vtx  = storage->get_data<event_vertex> (_vertexProducer );

    larlite::event_hit* ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());

    if (!ev_hit){
      print(larlite::msg::kWARNING,__FUNCTION__,"no hits");
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

      // get coordinates of hits to calculate linearity
      std::vector<double> hit_w_v;
      std::vector<double> hit_t_v;

      _pl = ev_hit->at(hit_idx_v[0]).WireID().Plane;

      if (ev_hit->at(hit_idx_v[0]).GoodnessOfFit() < 0) continue;

      double dvtx_min_sq = 1e6;
      double dvtx_max_sq = 0.;
      for (auto const& hit_idx : hit_idx_v){
	auto const& hit = ev_hit->at(hit_idx);
	double wpt = hit.WireID().Wire  * _wire2cm;
	double tpt = hit.PeakTime()     * _time2cm;
	hit_w_v.push_back( wpt );
	hit_t_v.push_back( tpt );
	double dd = ( ( (wpt - _vtx_w_cm[_pl]) * (wpt - _vtx_w_cm[_pl]) ) +
		      ( (tpt - _vtx_t_cm[_pl]) * (tpt - _vtx_t_cm[_pl]) ) );
	if (dd < dvtx_min_sq) dvtx_min_sq = dd;
	if (dd > dvtx_max_sq) dvtx_max_sq = dd;
      }

      _dvtx_max = sqrt(dvtx_max_sq);
      _dvtx_min = sqrt(dvtx_min_sq);

      twodimtools::Linearity lin(hit_w_v,hit_t_v);

      _nhits               = hit_w_v.size();
      _lin                 = lin._lin;
      _slope               = lin._slope;
      _local_lin_avg       = lin._local_lin_avg;
      _local_lin_truncated = lin._local_lin_truncated_avg;
      _ssv                 = lin._summed_square_variance;

      _tree->Fill();

    }// for all clusters


    return true;
  }

  bool ClusterLinearityStudy::finalize() {

    _fout->cd();
    _tree->Write();

    return true;
  }
  
  bool ClusterLinearityStudy::loadVertex(event_vertex* ev_vtx) {
    
    if (ev_vtx->size() != 1) return false;
    
    // get vertex position on each plane
    if ( (ev_vtx->size() == 1) ){
      auto const& vtx = ev_vtx->at(0);

      std::vector<double> xyz = {vtx.X(), vtx.Y(), vtx.Z()};
      
      auto geoH = larutil::GeometryHelper::GetME();
      auto geom = larutil::Geometry::GetME();

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
