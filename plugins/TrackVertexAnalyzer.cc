// -*- C++ -*-
//
// Package:    JetMetTesting/TrackVertexAnalyzer
// Class:      TrackVertexAnalyzer
// 
/**\class TrackVertexAnalyzer TrackVertexAnalyzer.cc JetMetTesting/TrackVertexAnalyzer/plugins/TrackVertexAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Young Ho Shin
//         Created:  Mon, 27 Jul 2015 13:21:10 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// data formats
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"
#include "CommonTools/Utils/interface/associationMapFilterValues.h"

#include "CommonTools/RecoUtils/interface/PFCand_AssoMapAlgos.h"

// tree/histogram writing
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TrackVertexAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit TrackVertexAnalyzer(const edm::ParameterSet&);
    ~TrackVertexAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    
    // ----------utility functions ---------------------------
    template <class T> bool vertexSelector(T t);

    // ----------member data ---------------------------
    // Config parameters
    bool UseAssociators;

    // TrackingParticle collections
    edm::EDGetTokenT<TrackingParticleCollection> label_tp_effic;
    edm::EDGetTokenT<TrackingParticleCollection> label_tp_fake;

    // reco::Track collections
    std::vector<edm::InputTag> label;
    std::vector<edm::EDGetTokenT<edm::View<reco::Track> > > labelToken;
    // std::vector<edm::EDGetTokenT<edm::View<TrajectorySeed> > > labelTokenSeed;

    // Track-TrackingParticle associators or associator maps
    std::vector<edm::InputTag> associators;
    std::vector<edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator>> associatorTokens;
    std::vector<edm::EDGetTokenT<reco::SimToRecoCollection>> associatormapStRs;
    std::vector<edm::EDGetTokenT<reco::RecoToSimCollection>> associatormapRtSs;

    // For vertex performance only
    edm::EDGetTokenT<TrackingVertexCollection> label_tv;
    edm::EDGetTokenT<edm::View<reco::Vertex> > recoVertexToken_;
    edm::EDGetTokenT<reco::VertexToTrackingVertexAssociator> vertexAssociatorToken_;

    edm::EDGetTokenT<PFCandToVertexAssMap> pfCandidateToVertexAssociationsToken_;

    // TTree objects
    TTree* outputTree;
    // TTree* tpTree;
    // TTree* rtTree;
    // Class containing all the objects to be saved
    struct outputClass {
      public:
        // per-TrackingVertex quantities
        std::vector<int>    tv_nMatch; // number of matched reco::Vertex
        std::vector<int>    tv_tag; // 0 = HS, 1 = PU
        std::vector<int>    tv_n_tracks;
        std::vector<double> tv_sumEt_tracks;
        std::vector<double> tv_sumPt2_tracks;
        std::vector<double> tv_qualitySum;
        std::vector<double> tv_x;
        std::vector<double> tv_y;
        std::vector<double> tv_z;
        // per-reco::Vertex quantities
        std::vector<int>    rv_nMatch; // number of matched TrackingVertex
        std::vector<int>    rv_n_tracks;
        std::vector<double> rv_sumEt_tracks;
        std::vector<double> rv_sumPt2_tracks;
        std::vector<double> rv_sumEt_tv; // track sumEt of matched TrackingVertex - Only when nMatch==1
        std::vector<double> rv_sumPt2_tv; // track sumPt2 of matched TrackingVertex - Only when nMatch==1
        std::vector<double> rv_x;
        std::vector<double> rv_y;
        std::vector<double> rv_z;
    };
    outputClass output;
    // const values for tv_tag
    static const int TV_TAG_UNDEFINED = -1; 
    static const int TV_TAG_HS = 0; // TrackingVertex originates from HS
    static const int TV_TAG_PU = 1; // TrackingVertex originates from PU
};

//
// constants, enums and typedefs
//
// from CommonTools/RecoUtils/plugins/PFCand_AssoMap.cc
typedef edm::AssociationMap<edm::OneToManyWithQuality< reco::VertexCollection, reco::PFCandidateCollection, int> > PFCandToVertexAssMap;

//
// static data member definitions
//

//
// constructors and destructor
//
TrackVertexAnalyzer::TrackVertexAnalyzer(const edm::ParameterSet& pset) :
  output( {std::vector<int>(), std::vector<int>() } )
{
  //now do what ever initialization is needed
  usesResource("TFileService");

  // TrackingParticle collections
  label_tp_effic = consumes<TrackingParticleCollection>(pset.getParameter< edm::InputTag >("label_tp_effic"));
  label_tp_fake = consumes<TrackingParticleCollection>(pset.getParameter< edm::InputTag >("label_tp_fake"));

  // reco::Track collections
  label = pset.getParameter< std::vector<edm::InputTag> >("label");
  bool isSeed = false;
  if (isSeed) {
    // for (auto itag : label) labelTokenSeed.push_back(iC.consumes<edm::View<TrajectorySeed> >(itag));
  } else {
    for (auto itag : label) labelToken.push_back(consumes<edm::View<reco::Track> >(itag));
  }

  // TrackingParticle-Track associators
  UseAssociators = pset.getParameter< bool >("UseAssociators");
  associators = pset.getUntrackedParameter< std::vector<edm::InputTag> >("associators");
  if(UseAssociators) {
    for (auto const& src: associators) {
      associatorTokens.push_back(consumes<reco::TrackToTrackingParticleAssociator>(src));
    }
  } else {   
    for (auto const& src: associators) {
      associatormapStRs.push_back(consumes<reco::SimToRecoCollection>(src));
      associatormapRtSs.push_back(consumes<reco::RecoToSimCollection>(src));
    }
  }

  // TrackingVertex and reco::Vertex collections
  label_tv = consumes<TrackingVertexCollection>(pset.getParameter< edm::InputTag >("label_tv"));
  recoVertexToken_ = consumes<edm::View<reco::Vertex> >(pset.getUntrackedParameter<edm::InputTag>("label_vertex"));
  // TrackingVertex-Vertex associator
  vertexAssociatorToken_ = consumes<reco::VertexToTrackingVertexAssociator>(pset.getUntrackedParameter<edm::InputTag>("vertexAssociator"));

  ////////////////////////////////////////////////////////////
  // Inputs for track-vertex algorithm performance monitoring
  ////////////////////////////////////////////////////////////
  pfCandidateToVertexAssociationsToken_ = consumes<PFCandToVertexAssMap>(pset.getParameter<edm::InputTag>("srcPFCandidateToVertexAssociations"));

  edm::Service<TFileService> fs;
  outputTree = fs->make<TTree>("TVATree", "TVATree");

  outputTree->Branch("tv_nMatch"        , &output.tv_nMatch        ) ;
  outputTree->Branch("tv_tag"           , &output.tv_tag           ) ;
  outputTree->Branch("tv_n_tracks"      , &output.tv_n_tracks      ) ;
  outputTree->Branch("tv_sumEt_tracks"  , &output.tv_sumEt_tracks  ) ;
  outputTree->Branch("tv_sumPt2_tracks" , &output.tv_sumPt2_tracks ) ;
  outputTree->Branch("tv_qualitySum"    , &output.tv_qualitySum    ) ;
  outputTree->Branch("tv_x"             , &output.tv_x             ) ;
  outputTree->Branch("tv_y"             , &output.tv_y             ) ;
  outputTree->Branch("tv_z"             , &output.tv_z             ) ;
  outputTree->Branch("rv_nMatch"        , &output.rv_nMatch        ) ;
  outputTree->Branch("rv_n_tracks"      , &output.rv_n_tracks      ) ;
  outputTree->Branch("rv_sumEt_tracks"  , &output.rv_sumEt_tracks  ) ;
  outputTree->Branch("rv_sumPt2_tracks" , &output.rv_sumPt2_tracks ) ;
  outputTree->Branch("rv_sumEt_tv"      , &output.rv_sumEt_tv      ) ;
  outputTree->Branch("rv_sumPt2_tv"     , &output.rv_sumPt2_tv     ) ;
  outputTree->Branch("rv_x"             , &output.rv_x             ) ;
  outputTree->Branch("rv_y"             , &output.rv_y             ) ;
  outputTree->Branch("rv_z"             , &output.rv_z             ) ;
}


TrackVertexAnalyzer::~TrackVertexAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
  void
TrackVertexAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  using namespace edm;
  using namespace reco;
  // Clear output trees
  output.tv_nMatch        .clear();
  output.tv_tag           .clear();
  output.tv_n_tracks      .clear();
  output.tv_sumEt_tracks  .clear();
  output.tv_sumPt2_tracks .clear();
  output.tv_qualitySum    .clear();
  output.tv_x             .clear();
  output.tv_y             .clear();
  output.tv_z             .clear();
  output.rv_nMatch        .clear();
  output.rv_n_tracks      .clear();
  output.rv_sumEt_tracks  .clear();
  output.rv_sumPt2_tracks .clear();
  output.rv_sumEt_tv      .clear();
  output.rv_sumPt2_tv     .clear();
  output.rv_x             .clear();
  output.rv_y             .clear();
  output.rv_z             .clear();

  LogTrace("TrackVertexAnalyzer") << "\n====================================================" << "\n"
    << "Analyzing new event" << "\n"
    << "====================================================\n" << "\n";
  ////////////////////////////////////////////////////////////
  // Begin: Track retrieval
  ////////////////////////////////////////////////////////////
  ////
  // TrackingParticle collections
  //
  edm::Handle<TrackingParticleCollection>  TPCollectionHeff ;
  event.getByToken(label_tp_effic,TPCollectionHeff);
  TrackingParticleCollection const & tPCeff = *(TPCollectionHeff.product());
  edm::Handle<TrackingParticleCollection>  TPCollectionHfake ;
  event.getByToken(label_tp_fake,TPCollectionHfake);
  ////
  // reco::Track collections and TrackingParticle-Track associators
  //
  bool ignoremissingtkcollection_ = false;
  int w=0; //counter counting the number of sets of histograms
  int ww = 0, www =0;
  //
  //get collections from the event
  //
  edm::Handle<View<Track> >  trackCollection;
  // if( !event.getByToken(labelToken[www], trackCollection) && ignoremissingtkcollection_ ) continue;
  if( !event.getByToken(labelToken[www], trackCollection) && ignoremissingtkcollection_ ) return;

  reco::RecoToSimCollection const * recSimCollP=nullptr;
  reco::SimToRecoCollection const * simRecCollP=nullptr;
  reco::RecoToSimCollection recSimCollL;
  reco::SimToRecoCollection simRecCollL;

  //associate tracks
  LogTrace("TrackVertexAnalyzer") << "Analyzing "
    << label[www] << " with "
    << associators[ww] <<"\n";
  if(UseAssociators){
    edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
    event.getByToken(associatorTokens[ww], theAssociator);

    LogTrace("TrackVertexAnalyzer") << "Calling associateRecoToSim method" << "\n";
    recSimCollL = std::move(theAssociator->associateRecoToSim(trackCollection,
          TPCollectionHfake));
    recSimCollP = &recSimCollL;
    LogTrace("TrackVertexAnalyzer") << "Calling associateSimToReco method" << "\n";
    simRecCollL = std::move(theAssociator->associateSimToReco(trackCollection,
          TPCollectionHeff));
    simRecCollP = &simRecCollL;
  }
  else{
    Handle<reco::SimToRecoCollection > simtorecoCollectionH;
    event.getByToken(associatormapStRs[ww], simtorecoCollectionH);
    simRecCollP = simtorecoCollectionH.product();

    // We need to filter the associations of the current track
    // collection only from SimToReco collection, otherwise the
    // SimToReco histograms get false entries
    simRecCollL = associationMapFilterValues(*simRecCollP, *trackCollection);
    simRecCollP = &simRecCollL;

    Handle<reco::RecoToSimCollection > recotosimCollectionH;
    event.getByToken(associatormapRtSs[ww],recotosimCollectionH);
    recSimCollP = recotosimCollectionH.product();

    // In general, we should filter also the RecoToSim collection.
    // But, that would require changing the input type of TPs to
    // View<TrackingParticle>, and either replace current
    // associator interfaces with (View<Track>, View<TP>) or
    // adding the View,View interface (same goes for
    // RefToBaseVector,RefToBaseVector). Since there is currently
    // no compelling-enough use-case, we do not filter the
    // RecoToSim collection here. If an association using a subset
    // of the Sim collection is needed, user has to produce such
    // an association explicitly.
  }

  reco::RecoToSimCollection const & recSimColl = *recSimCollP;
  reco::SimToRecoCollection const & simRecColl = *simRecCollP;
  w++;
  ////////////////////////////////////////////////////////////
  // End: Track retrieval
  ////////////////////////////////////////////////////////////

  // Tracking Vertex collection
  edm::Handle<TrackingVertexCollection> htv;
  event.getByToken(label_tv, htv);
  // reco::Vertex collection
  edm::Handle<edm::View<reco::Vertex> > hvertex;
  event.getByToken(recoVertexToken_, hvertex);

  edm::Handle<reco::VertexToTrackingVertexAssociator> hvassociator;
  event.getByToken(vertexAssociatorToken_, hvassociator);

  // Get RECO-to-SIM vertex association
  LogTrace("TrackVertexAnalyzer") << "Get RECO-to-SIM vertex association";
  auto v_r2s = hvassociator->associateRecoToSim(hvertex, htv);
  LogTrace("TrackVertexAnalyzer") << "hvertex->size(): " << hvertex->size() << std::endl;
  for ( unsigned int irv = 0; irv != hvertex->size(); irv++ ) {
    LogTrace("TrackVertexAnalyzer") << "Matching to RECO vertex: " << irv << std::endl;
    int nMatch = 0;
    auto rvRef = hvertex->refAt(irv);
    // if ( !vertexSelector<decltype(rvRef)>(rvRef) ) continue;
    if ( !vertexSelector(rvRef) ) continue;
    LogTrace("TrackVertexAnalyzer") << "rvRef == NULL: " << rvRef.isNull() << std::endl;
    // Skip unmatched RECO vertices
    TrackingVertexRef tv;
    if ( v_r2s.find(rvRef) != v_r2s.end() ) {
      auto vec_tv_quality = v_r2s[rvRef];
      for (const auto tv_quality : vec_tv_quality) {
        tv = tv_quality.first;
        LogTrace("TrackVertexAnalyzer") << tv->nDaughterTracks() << std::endl;
        LogTrace("TrackVertexAnalyzer") << tv->nSourceTracks() << std::endl;
        nMatch++;
      }
    }
    LogTrace("TrackVertexAnalyzer") << "Matched " << nMatch << " SIM vertices to RECO vertex: " << irv << std::endl;
    int n_tracks = 0;
    double sumEt = 0;
    double sumPt2 = 0;
    for ( auto track = rvRef->tracks_begin(); track != rvRef->tracks_end(); track++ ) {
      n_tracks++;
      sumEt += (*track)->pt();
      sumPt2 += TMath::Power( (*track)->pt(), 2 );
    }
    double sumEt_tv = 0;
    double sumPt2_tv = 0;
    if ( nMatch == 1 ) {
      for ( const auto track : tv->daughterTracks() ) {
        sumEt_tv += track->pt();
        sumPt2_tv += TMath::Power( track->pt(), 2 );
      }
    }
    output.rv_nMatch        .push_back ( nMatch     ) ;
    output.rv_n_tracks      .push_back ( n_tracks   ) ;
    output.rv_sumEt_tracks  .push_back ( sumEt      ) ;
    output.rv_sumPt2_tracks .push_back ( sumPt2     ) ;
    output.rv_sumEt_tv      .push_back ( sumEt_tv   ) ;
    output.rv_sumPt2_tv     .push_back ( sumPt2_tv  ) ;
    output.rv_x             .push_back ( rvRef->x() ) ;
    output.rv_y             .push_back ( rvRef->y() ) ;
    output.rv_z             .push_back ( rvRef->z() ) ;
  }

  // Get SIM-to-RECO vertex association
  auto v_s2r = hvassociator->associateSimToReco(hvertex, htv);
  int currentEvent = -1;
  LogTrace("TrackVertexAnalyzer") << "Number of SIM vertices in collection: " << htv->size() << std::endl;
  for ( unsigned int itv = 0; itv != htv->size(); itv++ ) {
    int nMatch = 0;
    int tag = TV_TAG_UNDEFINED;
    double sumEt = 0;
    double sumPt2 = 0;
    double qualitySum = 0;
    // auto tv = htv[itv];
    // edm::Ref<TrackingVertexCollection> tvRef(htv, itv);
    const TrackingVertexRef tvRef(htv, itv);
    // LogTrace("TrackVertexAnalyzer") << "BX: " << tvRef->eventId().bunchCrossing() << std::endl;
    // LogTrace("TrackVertexAnalyzer") << "event: " << tvRef->eventId().event() << std::endl;
    if(tvRef->eventId().bunchCrossing() != 0) continue;
    if(tvRef->eventId().event() != currentEvent) {
      currentEvent = tvRef->eventId().event();
      if ( currentEvent == 0 ) tag = TV_TAG_HS;
      else                     tag = TV_TAG_PU;
    }
    else {
      continue;
    }
    if ( v_s2r.find(tvRef) != v_s2r.end() ) {
      // Only run over primary vertices of the in-time pileup
      // events (BX=0, first vertex in each of the events)
      // LogTrace("TrackVertexAnalyzer") << "Matching to SIM vertex: " << itv << std::endl;
      // LogTrace("TrackVertexAnalyzer") << "event: " << tvRef->eventId().event() << std::endl;
      // LogTrace("TrackVertexAnalyzer") << "currentEvent: " << currentEvent << std::endl;
      const auto vec_rv_quality = v_s2r[tvRef];
      for (const auto rv_quality : vec_rv_quality) {
        auto rv = rv_quality.first;
        auto quality = rv_quality.second;
        if (!vertexSelector(rv)) continue;
        LogTrace("TrackVertexAnalyzer") << "SIM-RECO association quality: " << quality;
        qualitySum += quality;
        nMatch++;
      }
      for ( const auto track : tvRef->daughterTracks() ) {
        sumEt += track->pt();
        sumPt2 += TMath::Power( track->pt(), 2 );
      }
    }
    output . tv_nMatch        . push_back( nMatch );
    output . tv_tag           . push_back( tag );
    output . tv_n_tracks      . push_back( tvRef->nDaughterTracks() );
    output . tv_sumEt_tracks  . push_back( sumEt );
    output . tv_sumPt2_tracks . push_back( sumPt2 );
    output . tv_qualitySum    . push_back( qualitySum );
    output . tv_x             . push_back( tvRef->position().x() );
    output . tv_y             . push_back( tvRef->position().y() );
    output . tv_z             . push_back( tvRef->position().z() );
  }

  edm::Handle<PFCandToVertexAssMap> pfCandidateToVertexAssociations;
  event.getByToken(pfCandidateToVertexAssociationsToken_, pfCandidateToVertexAssociations);

  // Loop over vertex-PF association
  for ( PFCandToVertexAssMap::const_iterator vertexToCand = pfCandidateToVertexAssociations->begin();
    vertexToCand != pfCandidateToVertexAssociations->end(); ++vertexToCand ) {
    reco::VertexRef vertex = vertexToCand->key;
    const PFCandQualityPairVector& pfCandidates_vertex = vertexToCand->val;
    LogTrace("TrackVertexAnalyzer") << "Looking for pfCand->trackRef in recSimColl.";
    // float matchingEfficiency = 0;
    float nMatch = 0, nTracks = 0;
    for ( const auto pfCandidate_vertex : pfCandidates_vertex ) {
      const PFCandidateRef& pfRef = (pfCandidate_vertex).first;
      const reco::TrackRef& trackRef = pfRef->trackRef();
      auto trackBaseRef = RefToBase<reco::Track>(trackRef);
      if ( recSimColl.find(trackBaseRef) != recSimColl.end() ) {
        nMatch++;
        // LogTrace("TrackVertexAnalyzer") << "Retrieving vec_tp_quality.";
        auto vec_tp_quality = recSimColl[trackBaseRef];
      }
      nTracks++;
      // auto pfCandTrack = RefToBase<reco::Track>(*trackRef);
      // for ( const auto recSim = recSimColl.begin(); recSim != recSimColl.end(); recSim++ ) {
      //   Ref<>(recSim->key);
      // }
      // if ( recSimColl.find(pfCandTrack) != recSimColl.end() ) {
      //   LogTrace("TrackVertexAnalyzer") << "Found trackRef in recSimColl.\n";
      // }
      // else {
      //   LogTrace("TrackVertexAnalyzer") << "Found trackRef in recSimColl.\n";
      // }
    }
    LogTrace("TrackVertexAnalyzer") << "Found " << nMatch << " / " << nTracks << " tracks in recSimColl.";
  }


  // Write out outputTree
  outputTree->Fill();

  // // Get SIM-to-RECO vertex association
  // auto v_s2r = hvassociator->associateSimToReco(hvertex, htv);

  // // Get pv[0] from RECO. If it is not a good vertex, skip event.
  // auto pvPtr = hvertex->refAt(0);
  // if(pvPtr->isFake() || pvPtr->ndof() < 0) // skip junk vertices
  //   return;
  // // If pv[0] from RECO is NOT associated with a SIM vertex, skip event
  // auto pvFound = v_r2s.find(pvPtr);
  // if(pvFound == v_r2s.end())
  //   return;

  // // Get ParametersDefinerForTP, required to compute momentum/vertex-position for TrackingParticleRef
  // edm::ESHandle<ParametersDefinerForTP> parametersDefinerTPHandle;
  // setup.get<TrackAssociatorRecord>().get(parametersDefiner,parametersDefinerTPHandle);
  // //Since we modify the object, we must clone it
  // auto parametersDefinerTP = parametersDefinerTPHandle->clone();


  // edm::ESHandle<ParametersDefinerForTP> parametersDefinerTPHandle;
  // setup.get<TrackAssociatorRecord>().get(parametersDefiner,parametersDefinerTPHandle);
  // //Since we modify the object, we must clone it
  // auto parametersDefinerTP = parametersDefinerTPHandle->clone();
  //
  // edm::Handle<TrackingParticleCollection>  TPCollectionHeff ;
  // event.getByToken(label_tp_effic,TPCollectionHeff);
  // TrackingParticleCollection const & tPCeff = *(TPCollectionHeff.product());
  //
  // edm::Handle<TrackingParticleCollection>  TPCollectionHfake ;
  // event.getByToken(label_tp_fake,TPCollectionHfake);
  //
  //
  // if(parametersDefinerIsCosmic_) {
  //   edm::Handle<SimHitTPAssociationProducer::SimHitTPAssociationList> simHitsTPAssoc;
  //   //warning: make sure the TP collection used in the map is the same used in the MTV!
  //   event.getByToken(_simHitTpMapTag,simHitsTPAssoc);
  //   parametersDefinerTP->initEvent(simHitsTPAssoc);
  //   cosmictpSelector.initEvent(simHitsTPAssoc);
  // }
  //
  // if(doPlotsOnlyForTruePV_) {
  //   edm::Handle<TrackingVertexCollection> htv;
  //   event.getByToken(label_tv, htv);
  //
  //   edm::Handle<edm::View<reco::Vertex> > hvertex;
  //   event.getByToken(recoVertexToken_, hvertex);
  //
  //   edm::Handle<reco::VertexToTrackingVertexAssociator> hvassociator;
  //   event.getByToken(vertexAssociatorToken_, hvassociator);
  //
  //   auto v_r2s = hvassociator->associateRecoToSim(hvertex, htv);
  //   auto pvPtr = hvertex->refAt(0);
  //   if(pvPtr->isFake() || pvPtr->ndof() < 0) // skip junk vertices
  //     return;
  //   auto pvFound = v_r2s.find(pvPtr);
  //   if(pvFound == v_r2s.end())
  //     return;
  // }
  //
  // edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  // event.getByToken(bsSrc,recoBeamSpotHandle);
  // reco::BeamSpot const & bs = *recoBeamSpotHandle;
  //
  // edm::Handle< std::vector<PileupSummaryInfo> > puinfoH;
  // event.getByToken(label_pileupinfo,puinfoH);
  // PileupSummaryInfo puinfo;
  //
  // for (unsigned int puinfo_ite=0;puinfo_ite<(*puinfoH).size();++puinfo_ite){
  //   if ((*puinfoH)[puinfo_ite].getBunchCrossing()==0){
  //     puinfo=(*puinfoH)[puinfo_ite];
  //     break;
  //   }
  // }
  //
  // #<{(|
  // edm::Handle<TrackingVertexCollection> tvH;
  // event.getByToken(label_tv,tvH);
  // TrackingVertexCollection const & tv = *tvH;
  // |)}>#
  //
  // // Calculate the number of 3D layers for TPs
  // //
  // // I would have preferred to produce the ValueMap to Event and read
  // // it from there, but there would have been quite some number of
  // // knock-on effects, and again the fact that we take two TP
  // // collections do not support Ref<TP>'s would have complicated that.
  // //
  // // In principle we could use the SimHitTPAssociationList read above
  // // for parametersDefinerIsCosmic_, but since we don't currently
  // // support Ref<TP>s, we can't in general use it since eff/fake TP
  // // collections can, in general, be different.
  // TrackingParticleNumberOfLayers tpNumberOfLayersAlgo(event, simHitTokens_);
  // auto nlayers_tPCeff_ptrs = tpNumberOfLayersAlgo.calculate(TPCollectionHeff, setup);
  // const auto& nLayers_tPCeff = *(std::get<TrackingParticleNumberOfLayers::nTrackerLayers>(nlayers_tPCeff_ptrs));
  // const auto& nPixelLayers_tPCeff = *(std::get<TrackingParticleNumberOfLayers::nPixelLayers>(nlayers_tPCeff_ptrs));
  // const auto& nStripMonoAndStereoLayers_tPCeff = *(std::get<TrackingParticleNumberOfLayers::nStripMonoAndStereoLayers>(nlayers_tPCeff_ptrs));
  //
  // // Precalculate TP selection (for efficiency), and momentum and vertex wrt PCA
  // //
  // // TODO: ParametersDefinerForTP ESProduct needs to be changed to
  // // EDProduct because of consumes.
  // //
  // // In principle, we could just precalculate the momentum and vertex
  // // wrt PCA for all TPs for once and put that to the event. To avoid
  // // repetitive calculations those should be calculated only once for
  // // each TP. That would imply that we should access TPs via Refs
  // // (i.e. View) in here, since, in general, the eff and fake TP
  // // collections can be different (and at least HI seems to use that
  // // feature). This would further imply that the
  // // RecoToSimCollection/SimToRecoCollection should be changed to use
  // // View<TP> instead of vector<TP>, and migrate everything.
  // //
  // // Or we could take only one input TP collection, and do another
  // // TP-selection to obtain the "fake" collection like we already do
  // // for "efficiency" TPs.
  // std::vector<size_t> selected_tPCeff;
  // std::vector<std::tuple<TrackingParticle::Vector, TrackingParticle::Point>> momVert_tPCeff;
  // selected_tPCeff.reserve(tPCeff.size());
  // momVert_tPCeff.reserve(tPCeff.size());
  // int nIntimeTPs = 0;
  // if(parametersDefinerIsCosmic_) {
  //   for(size_t j=0; j<tPCeff.size(); ++j) {
  //     TrackingParticleRef tpr(TPCollectionHeff, j);
  //
  //     TrackingParticle::Vector momentum = parametersDefinerTP->momentum(event,setup,tpr);
  //     TrackingParticle::Point vertex = parametersDefinerTP->vertex(event,setup,tpr);
  //     if(doSimPlots_) {
  //       histoProducerAlgo_->fill_generic_simTrack_histos(momentum, vertex, tpr->eventId().bunchCrossing());
  //     }
  //     if(tpr->eventId().bunchCrossing() == 0)
  //       ++nIntimeTPs;
  //
  //     if(cosmictpSelector(tpr,&bs,event,setup)) {
  //       selected_tPCeff.push_back(j);
  //       momVert_tPCeff.emplace_back(momentum, vertex);
  //     }
  //   }
  // }
  // else {
  //   size_t j=0;
  //   for(auto const& tp: tPCeff) {
  //
  //     // TODO: do we want to fill these from all TPs that include IT
  //     // and OOT (as below), or limit to IT+OOT TPs passing tpSelector
  //     // (as it was before)? The latter would require another instance
  //     // of tpSelector with intimeOnly=False.
  //     if(doSimPlots_) {
  //       histoProducerAlgo_->fill_generic_simTrack_histos(tp.momentum(), tp.vertex(), tp.eventId().bunchCrossing());
  //     }
  //     if(tp.eventId().bunchCrossing() == 0)
  //       ++nIntimeTPs;
  //
  //     if(tpSelector(tp)) {
  //       selected_tPCeff.push_back(j);
  // TrackingParticleRef tpr(TPCollectionHeff, j);
  //       TrackingParticle::Vector momentum = parametersDefinerTP->momentum(event,setup,tpr);
  //       TrackingParticle::Point vertex = parametersDefinerTP->vertex(event,setup,tpr);
  //       momVert_tPCeff.emplace_back(momentum, vertex);
  //     }
  //     ++j;
  //   }
  // }
  // if(doSimPlots_) {
  //   histoProducerAlgo_->fill_simTrackBased_histos(nIntimeTPs);
  // }
  //
  // //calculate dR for TPs
  // float dR_tPCeff[tPCeff.size()];
  // {
  //   float etaL[tPCeff.size()], phiL[tPCeff.size()];
  //   for(size_t iTP: selected_tPCeff) {
  //     //calculare dR wrt inclusive collection (also with PU, low pT, displaced)
  //     auto const& tp2 = tPCeff[iTP];
  //     auto  && p = tp2.momentum();
  //     etaL[iTP] = etaFromXYZ(p.x(),p.y(),p.z());
  //     phiL[iTP] = atan2f(p.y(),p.x());
  //   }
  //   auto i=0U;
  //   for ( auto const & tp : tPCeff) {
  //     double dR = std::numeric_limits<double>::max();
  //     if(dRtpSelector(tp)) {//only for those needed for efficiency!
  //       auto  && p = tp.momentum();
  //       float eta = etaFromXYZ(p.x(),p.y(),p.z());
  //       float phi = atan2f(p.y(),p.x());
  //       for(size_t iTP: selected_tPCeff) {
  //         //calculare dR wrt inclusive collection (also with PU, low pT, displaced)
  //   if (i==iTP) {continue;}
  //         auto dR_tmp = reco::deltaR2(eta, phi, etaL[iTP], phiL[iTP]);
  //         if (dR_tmp<dR) dR=dR_tmp;
  //       }  // ttp2 (iTP)
  //     }
  //     dR_tPCeff[i++] = std::sqrt(dR);
  //   }  // tp
  // }
  //
  // edm::Handle<View<Track> >  trackCollectionForDrCalculation;
  // event.getByToken(labelTokenForDrCalculation, trackCollectionForDrCalculation);
  //
  // // dE/dx
  // // at some point this could be generalized, with a vector of tags and a corresponding vector of Handles
  // // I'm writing the interface such to take vectors of ValueMaps
  // std::vector<const edm::ValueMap<reco::DeDxData> *> v_dEdx;
  // if(dodEdxPlots_) {
  //   edm::Handle<edm::ValueMap<reco::DeDxData> > dEdx1Handle;
  //   edm::Handle<edm::ValueMap<reco::DeDxData> > dEdx2Handle;
  //   event.getByToken(m_dEdx1Tag, dEdx1Handle);
  //   event.getByToken(m_dEdx2Tag, dEdx2Handle);
  //   v_dEdx.push_back(dEdx1Handle.product());
  //   v_dEdx.push_back(dEdx2Handle.product());
  // }
  //
  // int w=0; //counter counting the number of sets of histograms
  // for (unsigned int ww=0;ww<associators.size();ww++){
  //   for (unsigned int www=0;www<label.size();www++){
  //     //
  //     //get collections from the event
  //     //
  //     edm::Handle<View<Track> >  trackCollection;
  //     if(!event.getByToken(labelToken[www], trackCollection)&&ignoremissingtkcollection_)continue;
  //
  //     reco::RecoToSimCollection const * recSimCollP=nullptr;
  //     reco::SimToRecoCollection const * simRecCollP=nullptr;
  //     reco::RecoToSimCollection recSimCollL;
  //     reco::SimToRecoCollection simRecCollL;
  //
  //     //associate tracks
  //     LogTrace("TrackVertexAnalyzer") << "Analyzing "
  //                                << label[www] << " with "
  //                                << associators[ww] <<"\n";
  //     if(UseAssociators){
  //       edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
  //       event.getByToken(associatorTokens[ww], theAssociator);
  //
  // LogTrace("TrackVertexAnalyzer") << "Calling associateRecoToSim method" << "\n";
  // recSimCollL = std::move(theAssociator->associateRecoToSim(trackCollection,
  //                                                                  TPCollectionHfake));
  //        recSimCollP = &recSimCollL;
  // LogTrace("TrackVertexAnalyzer") << "Calling associateSimToReco method" << "\n";
  // simRecCollL = std::move(theAssociator->associateSimToReco(trackCollection,
  //                                                                  TPCollectionHeff));
  //       simRecCollP = &simRecCollL;
  //     }
  //     else{
  // Handle<reco::SimToRecoCollection > simtorecoCollectionH;
  // event.getByToken(associatormapStRs[ww], simtorecoCollectionH);
  // simRecCollP = simtorecoCollectionH.product();
  //
  //       // We need to filter the associations of the current track
  //       // collection only from SimToReco collection, otherwise the
  //       // SimToReco histograms get false entries
  //       simRecCollL = associationMapFilterValues(*simRecCollP, *trackCollection);
  //       simRecCollP = &simRecCollL;
  //
  // Handle<reco::RecoToSimCollection > recotosimCollectionH;
  // event.getByToken(associatormapRtSs[ww],recotosimCollectionH);
  // recSimCollP = recotosimCollectionH.product();
  //
  //       // In general, we should filter also the RecoToSim collection.
  //       // But, that would require changing the input type of TPs to
  //       // View<TrackingParticle>, and either replace current
  //       // associator interfaces with (View<Track>, View<TP>) or
  //       // adding the View,View interface (same goes for
  //       // RefToBaseVector,RefToBaseVector). Since there is currently
  //       // no compelling-enough use-case, we do not filter the
  //       // RecoToSim collection here. If an association using a subset
  //       // of the Sim collection is needed, user has to produce such
  //       // an association explicitly.
  //     }
  //
  //     reco::RecoToSimCollection const & recSimColl = *recSimCollP;
  //     reco::SimToRecoCollection const & simRecColl = *simRecCollP;
  //
  //
  //
  //     // ########################################################
  //     // fill simulation histograms (LOOP OVER TRACKINGPARTICLES)
  //     // ########################################################
  //
  //     //compute number of tracks per eta interval
  //     //
  //     LogTrace("TrackVertexAnalyzer") << "\n# of TrackingParticles: " << tPCeff.size() << "\n";
  //     int ats(0);      //This counter counts the number of simTracks that are "associated" to recoTracks
  //     int st(0);       //This counter counts the number of simulated tracks passing the MTV selection (i.e. tpSelector(tp) )
  //     unsigned sts(0);   //This counter counts the number of simTracks surviving the bunchcrossing cut
  //     unsigned asts(0);  //This counter counts the number of simTracks that are "associated" to recoTracks surviving the bunchcrossing cut
  //
  //     //loop over already-selected TPs for tracking efficiency
  //     for(size_t i=0; i<selected_tPCeff.size(); ++i) {
  //       size_t iTP = selected_tPCeff[i];
  // TrackingParticleRef tpr(TPCollectionHeff, iTP);
  // const TrackingParticle& tp = tPCeff[iTP];
  //
  //       auto const& momVert = momVert_tPCeff[i];
  // TrackingParticle::Vector momentumTP;
  // TrackingParticle::Point vertexTP;
  //
  // double dxySim(0);
  // double dzSim(0);
  // double dR=dR_tPCeff[iTP];
  //
  // //---------- THIS PART HAS TO BE CLEANED UP. THE PARAMETER DEFINER WAS NOT MEANT TO BE USED IN THIS WAY ----------
  // //If the TrackingParticle is collison like, get the momentum and vertex at production state
  // if(!parametersDefinerIsCosmic_)
  //   {
  //     momentumTP = tp.momentum();
  //     vertexTP = tp.vertex();
  //     //Calcualte the impact parameters w.r.t. PCA
  //     const TrackingParticle::Vector& momentum = std::get<TrackingParticle::Vector>(momVert);
  //     const TrackingParticle::Point& vertex = std::get<TrackingParticle::Point>(momVert);
  //     dxySim = (-vertex.x()*sin(momentum.phi())+vertex.y()*cos(momentum.phi()));
  //     dzSim = vertex.z() - (vertex.x()*momentum.x()+vertex.y()*momentum.y())/sqrt(momentum.perp2())
  //       * momentum.z()/sqrt(momentum.perp2());
  //   }
  // //If the TrackingParticle is comics, get the momentum and vertex at PCA
  // else
  //   {
  //     momentumTP = std::get<TrackingParticle::Vector>(momVert);
  //     vertexTP = std::get<TrackingParticle::Point>(momVert);
  //     dxySim = (-vertexTP.x()*sin(momentumTP.phi())+vertexTP.y()*cos(momentumTP.phi()));
  //     dzSim = vertexTP.z() - (vertexTP.x()*momentumTP.x()+vertexTP.y()*momentumTP.y())/sqrt(momentumTP.perp2())
  //       * momentumTP.z()/sqrt(momentumTP.perp2());
  //   }
  // //---------- THE PART ABOVE HAS TO BE CLEANED UP. THE PARAMETER DEFINER WAS NOT MEANT TO BE USED IN THIS WAY ----------
  //
  //       //This counter counts the number of simulated tracks passing the MTV selection (i.e. tpSelector(tp) ), but only for in-time TPs
  //       if(tp.eventId().bunchCrossing() == 0) {
  //         st++;
  //       }
  //
  // // in the coming lines, histos are filled using as input
  // // - momentumTP
  // // - vertexTP
  // // - dxySim
  // // - dzSim
  //       if(!doSimTrackPlots_)
  //         continue;
  //
  // // ##############################################
  // // fill RecoAssociated SimTracks' histograms
  // // ##############################################
  // const reco::Track* matchedTrackPointer=0;
  // if(simRecColl.find(tpr) != simRecColl.end()){
  //   auto const & rt = simRecColl[tpr];
  //   if (rt.size()!=0) {
  //     ats++; //This counter counts the number of simTracks that have a recoTrack associated
  //     // isRecoMatched = true; // UNUSED
  //     matchedTrackPointer = rt.begin()->first.get();
  //     LogTrace("TrackVertexAnalyzer") << "TrackingParticle #" << st
  //                                      << " with pt=" << sqrt(momentumTP.perp2())
  //                                      << " associated with quality:" << rt.begin()->second <<"\n";
  //   }
  // }else{
  //   LogTrace("TrackVertexAnalyzer")
  //     << "TrackingParticle #" << st
  //     << " with pt,eta,phi: "
  //     << sqrt(momentumTP.perp2()) << " , "
  //     << momentumTP.eta() << " , "
  //     << momentumTP.phi() << " , "
  //     << " NOT associated to any reco::Track" << "\n";
  // }
  //
  //
  //
  //
  //       int nSimHits = tp.numberOfTrackerHits();
  //       int nSimLayers = nLayers_tPCeff[tpr];
  //       int nSimPixelLayers = nPixelLayers_tPCeff[tpr];
  //       int nSimStripMonoAndStereoLayers = nStripMonoAndStereoLayers_tPCeff[tpr];
  //       histoProducerAlgo_->fill_recoAssociated_simTrack_histos(w,tp,momentumTP,vertexTP,dxySim,dzSim,nSimHits,nSimLayers,nSimPixelLayers,nSimStripMonoAndStereoLayers,matchedTrackPointer,puinfo.getPU_NumInteractions(), dR);
  //         sts++;
  //         if(matchedTrackPointer)
  //           asts++;
  //         if(dRtpSelectorNoPtCut(tp)) {
  //           h_simul_coll_allPt[ww]->Fill(www);
  //           if (matchedTrackPointer) {
  //             h_assoc_coll_allPt[ww]->Fill(www);
  //           }
  //
  //           if(dRtpSelector(tp)) {
  //             h_simul_coll[ww]->Fill(www);
  //             if (matchedTrackPointer) {
  //               h_assoc_coll[ww]->Fill(www);
  //             }
  //           }
  //         }
  //
  //
  //
  //
  //     } // End  for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){
  //
  //     // ##############################################
  //     // fill recoTracks histograms (LOOP OVER TRACKS)
  //     // ##############################################
  //     if(!doRecoTrackPlots_)
  //       continue;
  //     LogTrace("TrackVertexAnalyzer") << "\n# of reco::Tracks with "
  //                                << label[www].process()<<":"
  //                                << label[www].label()<<":"
  //                                << label[www].instance()
  //                                << ": " << trackCollection->size() << "\n";
  //
  //     int sat(0); //This counter counts the number of recoTracks that are associated to SimTracks from Signal only
  //     int at(0); //This counter counts the number of recoTracks that are associated to SimTracks
  //     int rT(0); //This counter counts the number of recoTracks in general
  //
  //     //calculate dR for tracks
  //     float dR_trk[trackCollection->size()];
  //     int i=0;
  //     float etaL[trackCollectionForDrCalculation->size()];
  //     float phiL[trackCollectionForDrCalculation->size()];
  //     for (auto const & track2 : *trackCollectionForDrCalculation) {
  //        auto  && p = track2.momentum();
  //        etaL[i] = etaFromXYZ(p.x(),p.y(),p.z());
  //        phiL[i] = atan2f(p.y(),p.x());
  //        ++i;
  //     }
  //     for(View<Track>::size_type i=0; i<trackCollection->size(); ++i){
  // auto const &  track = (*trackCollection)[i];
  // auto dR = std::numeric_limits<float>::max();
  //       auto  && p = track.momentum();
  //       float eta = etaFromXYZ(p.x(),p.y(),p.z());
  //       float phi = atan2f(p.y(),p.x());
  // for(View<Track>::size_type j=0; j<trackCollectionForDrCalculation->size(); ++j){
  //   auto dR_tmp = reco::deltaR2(eta, phi, etaL[j], phiL[j]);
  //   if ( (dR_tmp<dR) & (dR_tmp>std::numeric_limits<float>::min())) dR=dR_tmp;
  // }
  // dR_trk[i] = std::sqrt(dR);
  //     }
  //
  //     for(View<Track>::size_type i=0; i<trackCollection->size(); ++i){
  //
  // RefToBase<Track> track(trackCollection, i);
  // rT++;
  //
  // bool isSigSimMatched(false);
  // bool isSimMatched(false);
  //       bool isChargeMatched(true);
  //       int numAssocRecoTracks = 0;
  // int nSimHits = 0;
  // double sharedFraction = 0.;
  //
  //       auto tpFound = recSimColl.find(track);
  //       isSimMatched = tpFound != recSimColl.end();
  //       if (isSimMatched) {
  //           const auto& tp = tpFound->val;
  //     nSimHits = tp[0].first->numberOfTrackerHits();
  //           sharedFraction = tp[0].second;
  //           if (tp[0].first->charge() != track->charge()) isChargeMatched = false;
  //           if(simRecColl.find(tp[0].first) != simRecColl.end()) numAssocRecoTracks = simRecColl[tp[0].first].size();
  //     at++;
  //     for (unsigned int tp_ite=0;tp_ite<tp.size();++tp_ite){
  //             TrackingParticle trackpart = *(tp[tp_ite].first);
  //       if ((trackpart.eventId().event() == 0) && (trackpart.eventId().bunchCrossing() == 0)){
  //              isSigSimMatched = true;
  //      sat++;
  //      break;
  //       }
  //           }
  //     LogTrace("TrackVertexAnalyzer") << "reco::Track #" << rT << " with pt=" << track->pt()
  //                                      << " associated with quality:" << tp.begin()->second <<"\n";
  // } else {
  //   LogTrace("TrackVertexAnalyzer") << "reco::Track #" << rT << " with pt=" << track->pt()
  //                                    << " NOT associated to any TrackingParticle" << "\n";
  // }
  //
  // double dR=dR_trk[i];
  // histoProducerAlgo_->fill_generic_recoTrack_histos(w,*track,bs.position(),isSimMatched,isSigSimMatched, isChargeMatched, numAssocRecoTracks, puinfo.getPU_NumInteractions(), nSimHits, sharedFraction,dR);
  //       h_reco_coll[ww]->Fill(www);
  //       if(isSimMatched) {
  //         h_assoc2_coll[ww]->Fill(www);
  //         if(numAssocRecoTracks>1) {
  //           h_looper_coll[ww]->Fill(www);
  //         }
  //         else if(!isSigSimMatched) {
  //           h_pileup_coll[ww]->Fill(www);
  //         }
  //       }
  //
  // // dE/dx
  // if (dodEdxPlots_) histoProducerAlgo_->fill_dedx_recoTrack_histos(w,track, v_dEdx);
  //
  //
  // //Fill other histos
  // if (!isSimMatched) continue;
  //
  // histoProducerAlgo_->fill_simAssociated_recoTrack_histos(w,*track);
  //
  // TrackingParticleRef tpr = tpFound->val.begin()->first;
  //
  // #<{(| TO BE FIXED LATER
  // if (associators[ww]=="trackAssociatorByChi2"){
  //   //association chi2
  //   double assocChi2 = -tp.begin()->second;//in association map is stored -chi2
  //   h_assochi2[www]->Fill(assocChi2);
  //   h_assochi2_prob[www]->Fill(TMath::Prob((assocChi2)*5,5));
  // }
  // else if (associators[ww]=="quickTrackAssociatorByHits"){
  //   double fraction = tp.begin()->second;
  //   h_assocFraction[www]->Fill(fraction);
  //   h_assocSharedHit[www]->Fill(fraction*track->numberOfValidHits());
  // }
  // |)}>#
  //
  //
  // //Get tracking particle parameters at point of closest approach to the beamline
  // TrackingParticle::Vector momentumTP = parametersDefinerTP->momentum(event,setup,tpr);
  // TrackingParticle::Point vertexTP = parametersDefinerTP->vertex(event,setup,tpr);
  // int chargeTP = tpr->charge();
  //
  // histoProducerAlgo_->fill_ResoAndPull_recoTrack_histos(w,momentumTP,vertexTP,chargeTP,
  //                                                   *track,bs.position());
  //
  //
  // //TO BE FIXED
  // //std::vector<PSimHit> simhits=tpr.get()->trackPSimHit(DetId::Tracker);
  // //nrecHit_vs_nsimHit_rec2sim[w]->Fill(track->numberOfValidHits(), (int)(simhits.end()-simhits.begin() ));
  //
  //     } // End of for(View<Track>::size_type i=0; i<trackCollection->size(); ++i){
  //
  //     histoProducerAlgo_->fill_trackBased_histos(w,at,rT,st);
  //
  //     LogTrace("TrackVertexAnalyzer") << "Total Simulated: " << st << "\n"
  //                                << "Total Associated (simToReco): " << ats << "\n"
  //                                << "Total Reconstructed: " << rT << "\n"
  //                                << "Total Associated (recoToSim): " << at << "\n"
  //                                << "Total Fakes: " << rT-at << "\n";
  //
  //     w++;
  //   } // End of  for (unsigned int www=0;www<label.size();www++){
  // } //END of for (unsigned int ww=0;ww<associators.size();ww++){
  //
  //
}


// ------------ method called once each job just before starting event loop  ------------
  void 
TrackVertexAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
  void 
TrackVertexAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackVertexAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

template <class T>
bool
TrackVertexAnalyzer::vertexSelector(T vertex) {
  bool vertexIsGood = (!vertex->isFake()) && (vertex->ndof()>4) && (fabs(vertex->z()<=24.0)) && (vertex->position().Rho()<=2.0);
  return vertexIsGood;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackVertexAnalyzer);
