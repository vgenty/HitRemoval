#ifndef LARLITE_CALCULATEHITREMOVALEFF_CXX
#define LARLITE_CALCULATEHITREMOVALEFF_CXX

#include "CalculateHitRemovalEff.h"
#include "DataFormat/hit.h"
#include "DataFormat/cluster.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/vertex.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"

namespace larlite {

  CalculateHitRemovalEff::CalculateHitRemovalEff()
    : _tree(nullptr)
  {
    _name = "CalculateHitRemovalEff";
    _fout = 0;
    _hit_amp = 0.;
    _entry = 0;
  }

  bool CalculateHitRemovalEff::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree("_tre","hit removal tree");

    _tree->Branch("_event",&_event,"event/I");
    _tree->Branch("_entry",&_entry,"entry/I");
    _tree->Branch("_run",&_run,"run/I");

    _tree->Branch("_nshr",&_nshr,"nshr/I");
    _tree->Branch("_edepshr",&_edepshr,"edepshr/D");
    
    _tree->Branch("_qtot",&_qtot,"qtot/D");
    _tree->Branch("_qremoved",&_qremoved,"qremoved/D");
    _tree->Branch("_ntot",&_ntot,"ntot/D");
    _tree->Branch("_nremoved",&_nremoved,"nremoved/D");
    _tree->Branch("_frac",&_frac,"frac/D");

    _tree->Branch("_qtot_0",&_qtot_0,"qtot_0/D");
    _tree->Branch("_qremoved_0",&_qremoved_0,"qremoved_0/D");
    _tree->Branch("_ntot_0",&_ntot_0,"ntot_0/I");
    _tree->Branch("_nremoved_0",&_nremoved_0,"nremoved_0/I");
    _tree->Branch("_frac_0",&_frac_0,"frac_0/D");

    _tree->Branch("_qtot_1",&_qtot_1,"qtot_1/D");
    _tree->Branch("_qremoved_1",&_qremoved_1,"qremoved_1/D");
    _tree->Branch("_ntot_1",&_ntot_1,"ntot_1/I");
    _tree->Branch("_nremoved_1",&_nremoved_1,"nremoved_1/I");
    _tree->Branch("_frac_1",&_frac_1,"frac_1/D");

    _tree->Branch("_qtot_2",&_qtot_2,"qtot_2/D");
    _tree->Branch("_qremoved_2",&_qremoved_2,"qremoved_2/D");
    _tree->Branch("_ntot_2",&_ntot_2,"ntot_2/I");
    _tree->Branch("_nremoved_2",&_nremoved_2,"nremoved_2/I");
    _tree->Branch("_frac_2",&_frac_2,"frac_2/D");

     _tree->Branch("_nclus_0_v","std::vector<int>",&_nclus_0_v);
    _tree->Branch("_nclus_1_v","std::vector<int>",&_nclus_1_v);
    _tree->Branch("_nclus_2_v","std::vector<int>",&_nclus_2_v);


    _tree->Branch("_nc" ,&_nc ,"nc/I");
    _tree->Branch("_pi0",&_pi0,"pi0/I");
    _tree->Branch("_pi0E",&_pi0E,"pi0E/D");

    _vtx_w_cm = std::vector<double> {0.,0.,0.};
    _vtx_t_cm = std::vector<double> {0.,0.,0.};

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();

    return true;
  }
  
  bool CalculateHitRemovalEff::analyze(storage_manager* storage) {

    _nclus_0_v.clear();
    _nclus_1_v.clear();
    _nclus_2_v.clear();

    // load vertex
    auto ev_vtx  = storage->get_data<event_vertex>(_vertexProducer);

    if (loadVertex(ev_vtx) == false) {
      print(larlite::msg::kERROR,__FUNCTION__,"num. vertices != 1");
      return false;
    }

    larlite::event_cluster *ev_clus = nullptr;
    larlite::event_hit     *ev_hit  = nullptr;
    std::vector< std::vector< unsigned int> > ass_cluster_hit_v;
    
    if (_use_cluster) {
      ev_clus = storage->get_data<event_cluster>("sc");
      ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());
      for (auto const& clus_ass : ass_cluster_hit_v)
	_nclus_2_v.push_back( clus_ass.size() );
    }
    else {
      ev_hit = storage->get_data<event_hit>("gaushit");
    }

    _nshr = 0;
    _edepshr = 0.;

    //auto ev_hit = storage->get_data<event_hit>("gaushit");
    
    auto ev_mcshr = storage->get_data<event_mcshower>("mcreco");


    auto ev_mctruth= storage->get_data<event_mctruth>("generator");

    _run   = storage->run_id();
    _event = storage->event_id();

    if (!ev_mctruth) { std::cout << "no truth" << std::endl; }

    if (_use_truth) {

      auto nu = ev_mctruth->at(0).GetNeutrino();
      auto parts = ev_mctruth->at(0).GetParticles();

      _nc = nu.CCNC();
      _pi0 = 0.;
      _pi0E = 0.;

      for (size_t i=0; i < parts.size(); i++){
	auto const& part = parts[i];
	if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) ){
	  _pi0E = part.Trajectory().at(0).E();
	  _pi0 += 1;
	}
      }// for all mcparticles
      
      
      _nshr = ev_mcshr->size();
      for (size_t i=0; i < _nshr; i++)
	_edepshr += ev_mcshr->at(i).DetProfile().E();
      
    }// if use truth

    _qtot = _qremoved = 0.;
    _qtot_0 = _qremoved_0 = 0.;
    _qtot_1 = _qremoved_1 = 0.;
    _qtot_2 = _qremoved_2 = 0.;

    _ntot = _nremoved = 0;
    _ntot_0 = _nremoved_0 = 0;
    _ntot_1 = _nremoved_1 = 0;
    _ntot_2 = _nremoved_2 = 0;

    for (size_t i=0; i < ev_hit->size(); i++) {

      auto const& hit = ev_hit->at(i);

      auto q   = hit.Integral();
      auto pl  = hit.WireID().Plane;
      auto GoF = hit.GoodnessOfFit();
      auto amp = hit.PeakAmplitude();
      
      auto const& wcm = fabs( (hit.WireID().Wire * _wire2cm) - _vtx_w_cm[pl]  );
      auto const& tcm = fabs( (hit.PeakTime()    * _time2cm) - _vtx_t_cm[pl]  );

      // out of ROI?
      if ( (wcm > _roi) or (tcm > _roi) ) continue;
      
      _qtot += q;
      _ntot += 1;
      
      if (pl == 0) { _qtot_0 += q; _ntot_0 += 1; }
      if (pl == 1) { _qtot_1 += q; _ntot_1 += 1; }
      if (pl == 2) { _qtot_2 += q; _ntot_2 += 1; }
      
      if ( (GoF < 0) || (amp < _hit_amp)) {
	_qremoved += q;
	_nremoved += 1;
	if (pl == 0) { _qremoved_0 += q; _nremoved_0 += 1; }
	if (pl == 1) { _qremoved_1 += q; _nremoved_1 += 1; }
	if (pl == 2) { _qremoved_2 += q; _nremoved_2 += 1; }
      }// if removed
      
    }// for all hits

    _frac = _qremoved/_qtot;
    _frac_0 = _qremoved_0/_qtot_0;
    _frac_1 = _qremoved_1/_qtot_1;
    _frac_2 = _qremoved_2/_qtot_2;

    _tree->Fill();

    _entry += 1;
  
    return true;
  }

  bool CalculateHitRemovalEff::finalize() {

    if (_fout) {
      _fout->cd();
      _tree->Write();
    }
  
    return true;
  }

  bool CalculateHitRemovalEff::loadVertex(event_vertex* ev_vtx) {
      
    if ( (!ev_vtx) || (ev_vtx->size() != 1) ) return false;
    
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
