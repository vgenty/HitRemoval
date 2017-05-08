#ifndef LARLITE_CLUSTERLINEARITYSTUDY_CXX
#define LARLITE_CLUSTERLINEARITYSTUDY_CXX

#include "ClusterLinearityStudy.h"

#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"

namespace larlite {

  bool ClusterLinearityStudy::initialize() {

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();

    if (_tree) delete _tree;
    _tree = new TTree("_tree","tree");
    _tree->Branch("_pl",&_pl,"pl/I");
    _tree->Branch("_nhits",&_nhits,"nhits/I");
    _tree->Branch("_lin",  &_lin  ,"lin/D"  );
    _tree->Branch("_ssv",  &_ssv  ,"ssv/D"  );
    _tree->Branch("_slope",  &_slope  ,"slope/D"  );
    _tree->Branch("_local_lin_truncated",  &_local_lin_truncated  ,"local_lin_truncated/D"  );
    _tree->Branch("_local_lin_avg",  &_local_lin_avg  ,"local_lin_avg/D"  );

    return true;
  }
  
  bool ClusterLinearityStudy::analyze(storage_manager* storage) {
  
    auto ev_clus = storage->get_data<event_cluster>(_clusterProducer);

    larlite::event_hit* ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());

    // loop trhough each cluster and calculate linaerity
    // if above some thresdhold, remove cluster

    for (size_t i=0; i < ass_cluster_hit_v.size(); i++){

      auto hit_idx_v = ass_cluster_hit_v[i];

      // get coordinates of hits to calculate linearity
      std::vector<double> hit_w_v;
      std::vector<double> hit_t_v;

      _pl = ev_hit->at(hit_idx_v[0]).WireID().Plane;
      
      for (auto const& hit_idx : hit_idx_v){
	hit_w_v.push_back( ev_hit->at(hit_idx).WireID().Wire  * _wire2cm );
	hit_t_v.push_back( ev_hit->at(hit_idx).PeakTime()     * _time2cm );
      }

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

}
#endif
