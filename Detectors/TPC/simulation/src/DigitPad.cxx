/// \file DigitPad.cxx
/// \author Andi Mathis, andreas.mathis@ph.tum.de

#include "TPCSimulation/DigitPad.h"
#include "TPCSimulation/Digitizer.h"
#include "TPCSimulation/Digit.h"

using namespace AliceO2::TPC;

DigitPad::DigitPad(Int_t pad)
  : mPad(pad),
    mTotalChargePad(0.),
    mADCCounts()
{}

DigitPad::~DigitPad()
{
  mADCCounts.resize(0);
  mTotalChargePad = 0;
}

void DigitPad::fillOutputContainer(TClonesArray *output, Int_t cru, Int_t timeBin, Int_t row, Int_t pad)
{  
  for(auto &aADCCounts : mADCCounts) {
    mTotalChargePad += aADCCounts.getADC();
  }
  
  const Float_t mADC = Digitizer::ADCvalue(mTotalChargePad);
  
  if(mADC > 0) {
    TClonesArray &clref = *output;
    new(clref[clref.GetEntriesFast()]) Digit(cru, mADC, row, pad, timeBin);
  }
}

void DigitPad::fillOutputContainer(TClonesArray *output, Int_t cru, Int_t timeBin, Int_t row, Int_t pad, Float_t commonMode)
{  
  for(auto &aADCCounts : mADCCounts) {
    mTotalChargePad += aADCCounts.getADC();
  }
  
  const Float_t mADC = Digitizer::ADCvalue(mTotalChargePad);
  
  if(mADC > 0) {
    TClonesArray &clref = *output;
    new(clref[clref.GetEntriesFast()]) Digit(cru, mADC, row, pad, timeBin, commonMode);
  }
}

void DigitPad::processCommonMode(Int_t cru, Int_t timeBin, Int_t row, Int_t pad)
{  
  for(auto &aADCCounts : mADCCounts) {
    mTotalChargePad += aADCCounts.getADC();
  }
}
