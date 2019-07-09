//
// Simple class to test the concept of fittting noisy and lossy data to a square wave by
// an exhaustive search of templates over a finite set of bins.  This uses the bitset class
// to make an extremely fast comparison of data with templates (xor).
// testing makes use of the ROOT library (see https://root.cern.ch)
// Original author: David Brown (LBNL), 9 July 2019
//
#include <bitset>
#include <vector>
#include <iostream>
#include <chrono>
#include "TRandom3.h"
#include "TH1F.h"
using namespace chrono;

// binary square wave on a measurement space from 0 to 1
// lambda is the wavelength
// phase is the phase (= position of leading edge of first wave)
// width is the wave width

bool squarewave(float lambda, float phase, float width, float x){
  float val = fmod(x-phase,lambda);
  if(val<0) val += lambda;
  return (val >=0.0 && val < width);
}

// test class
static const size_t nbits(36);
class SquareWaveFitTest {
  public:
  typedef std::bitset<nbits> bset;
  typedef std::vector<bset> bsvec;

  
  TH1F* _hover, *_hmod, *_hlambda, *_hphase, *_hwidth, *_hdur;
  bsvec _bitmodels;
  std::vector<float> _lvec, _p0vec, _fvec;
  TRandom3 _rand;
  float _fnbits;
  
  SquareWaveFitTest();
  int test(float lambda, float phase, float width, float eff,float pur,unsigned ntrials);
  unsigned bestOverlap(bset test, unsigned& bover);
  void setRandom(bset const& model, bset& rbits,float eff,float pur);
  bool badParams(float lambda, float phase, float width);
  void setbits(bset& bits,float lambda, float phase, float width);
};

int SquareWaveFitTest::test(float lambda, float phase, float width, float eff,float pur,unsigned ntrials){
  if(badParams(lambda,phase,width)){
    std::cout << "Illegal values " << std::endl;
    return -1;
  }
  bset model;
  setbits(model,lambda,phase,width);
  std::cout << model << std::endl;
// find the best template match
  unsigned bover;
  unsigned imod = bestOverlap(model,bover);
  std::cout << "template matches model " << imod << " overlap " << bover
  << " bits " << _bitmodels[imod] << std::endl;
  // generate random profile according to the model
  bset rbits;
  setRandom(model,rbits,eff,pur);
  std::cout << "Random profile = " << rbits << std::endl;
  imod = bestOverlap(rbits,bover);
  std::cout << "random bits matches model " << imod << " overlap " << bover
  << " bits " << _bitmodels[imod] << std::endl;
  _hdur->Reset();
  _hover->Reset();
  _hmod->Reset();
  _hlambda->Reset();
  _hphase->Reset();
  _hwidth->Reset();
  // now, high-stats tests
  for(unsigned itrial=0;itrial<ntrials;itrial++){
    setRandom(model,rbits,eff,pur);
    auto start = high_resolution_clock::now();
    imod = bestOverlap(rbits,bover);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    _hdur->Fill(duration.count());
    _hover->Fill(bover);
    _hmod->Fill(imod);
    _hlambda->Fill(_lvec[imod]);
    _hphase->Fill(_p0vec[imod]);
    _hwidth->Fill(_fvec[imod]);
  }
  return 0;
}

SquareWaveFitTest::SquareWaveFitTest() : _fnbits(nbits) {
// book histograms
  _hover = new TH1F("hover","Best Overlap",nbits+1,-0.5,nbits+0.5);
  _hmod = new TH1F("hmod","Best Model",_bitmodels.size()+1,-0.5,_bitmodels.size()+0.5);
  _hlambda = new TH1F("hlambda","Wavelength;#Lambda",100,0.0,0.6);
  _hphase = new TH1F("hphase","Phase",100,0.0,1.0);
  _hwidth = new TH1F("hwidth","Width",100,0.0,1.0);
  _hdur = new TH1F("hdur","Search Time;#mu seconds",51,-0.5,50.5);
  // random engine
  _rand = TRandom3(1238123);
  // build the models
  unsigned nval =nbits*(nbits+1)*(2*nbits+1)/64;
  _bitmodels.reserve(nval);
  _lvec.reserve(nval);
  _p0vec.reserve(nval);
  _fvec.reserve(nval);
  // exhaustive loop over reasonable models
  for(size_t ilambda=2;ilambda<nbits/2;ilambda++){
    float lambda = ilambda/_fnbits;
    for(int iphase=0; iphase <= (int)ilambda; iphase++){
      float phase = iphase/_fnbits;
      for(size_t iwidth=1; iwidth < ilambda-1; iwidth++){
	float width = iwidth/_fnbits;
	bset bits;
	setbits(bits,lambda,phase,width);
	_bitmodels.push_back(bits);
	_lvec.push_back(lambda);
	_p0vec.push_back(phase);
	_fvec.push_back(width);
      }
    }
  }
  std::cout << "N models = " << _bitmodels.size() << std::endl;
}

unsigned SquareWaveFitTest::bestOverlap(bset test, unsigned& bover) {
  unsigned retval(0); 
  bover = nbits;
  for(size_t imodel=0;imodel < _bitmodels.size();imodel++){
    size_t overlap = (_bitmodels[imodel]^test).count();
    if(overlap < bover){
      bover = overlap;
      retval = imodel;
    }
  }
  return retval;
}

void SquareWaveFitTest::setRandom(bset const& model, bset& rbits,float eff,float pur) {
  for(size_t ibit=0;ibit<nbits; ibit++){
    double rval = _rand.Uniform();
    if(model.test(ibit))
      rbits.set(ibit,rval < eff);
    else
      rbits.set(ibit,rval > pur);
  }
}
bool SquareWaveFitTest::badParams(float lambda, float phase, float width) {
  return lambda < 2.0/_fnbits || lambda > 0.5 ||
    phase < 0.0 || phase > lambda ||
    width < 1.0/_fnbits || width > lambda-1.0/_fnbits;
}

void SquareWaveFitTest::setbits(bset& bits,float lambda, float phase, float width) {
  // set bits for original model
  for(size_t ibit = 0; ibit < nbits; ibit++){
    float x = (float(ibit)+0.5)/float(nbits);
    bits.set(ibit,squarewave(lambda,phase,width,x));
  }
}


