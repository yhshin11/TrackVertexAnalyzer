#TrackVertexAnalyzer
Uses information from MC-truth (TrackingVertex, TrackingParticle, etc.) to analyze JetMET-relevant track/vertex performance. 

Modified from [MultiTrackValidator](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMultiTrackValidator).

#Instructions
```
# Checkout compatible CMSSW release and setup environment
cmsrel CMSSW_7_6_0_pre1
cd CMSSW_7_6_0_pre1/src
cmsenv
git cms-init
# Checkout MultiTrackValidator package
git cms-addpkg Validation/RecoTrack/
# Checkout TrackVertexAnalyzer in desired location
mkdir JetMetTesting; cd JetMetTesting
git clone https://github.com/yhshin11/TrackVertexAnalyzer.git
# Build and run TrackVertexAnalyzer
USER_CXXFLAGS="-DEDM_ML_DEBUG -g" scram b -v -j8
cmsRun TrackVertexAnalyzer/test/TrackVertexAnalyzer_cfg.py	
```
