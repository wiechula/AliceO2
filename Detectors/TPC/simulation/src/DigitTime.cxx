/// \file DigitTime.cxx
/// \author Andi Mathis, andreas.mathis@ph.tum.de

#include "TPCSimulation/DigitTime.h"
#include "TPCSimulation/DigitRow.h"
#include "TPCBase/Mapper.h"

using namespace AliceO2::TPC;

DigitTime::DigitTime(Int_t timeBin, Int_t nrows)
  : mTotalChargeTimeBin(0.),
    mTimeBin(timeBin),
    mRows(nrows)
{}

DigitTime::~DigitTime()
{
  for(auto &aRow : mRows) {
    if(aRow == nullptr) continue;
    delete aRow;
  }
}

void DigitTime::setDigit(Int_t eventID, Int_t trackID, Int_t cru, Int_t row, Int_t pad, Float_t charge)
{
  DigitRow *result = mRows[row];
  if(result != nullptr) {
    mRows[row]->setDigit(eventID, trackID, pad, charge);
  }
  else{
    const Mapper& mapper = Mapper::instance();
    mRows[row] = new DigitRow(row, mapper.getPadRegionInfo(CRU(cru).region()).getPadsInRowRegion(row));
    mRows[row]->setDigit(eventID, trackID, pad, charge);
  }
}

void DigitTime::fillOutputContainer(TClonesArray *output, Int_t cru, Int_t timeBin)
{
  for(auto &aRow : mRows) {
    if(aRow == nullptr) continue;
    aRow->fillOutputContainer(output, cru, timeBin, aRow->getRow());
  }
}

void DigitTime::fillOutputContainer(TClonesArray *output, Int_t cru, Int_t timeBin, std::vector<CommonMode> &commonModeContainer)
{
  Float_t commonMode =0;
  for (auto &aCommonMode :commonModeContainer){
    if(aCommonMode.getCRU() == cru && aCommonMode.getTimeBin() == timeBin) {
      commonMode = aCommonMode.getCommonMode();
      break;
    }
  }

  for(auto &aRow : mRows) {
    if(aRow == nullptr) continue;
    aRow->fillOutputContainer(output, cru, timeBin, aRow->getRow(), commonMode);
  }
}

void DigitTime::processCommonMode(Int_t cru, Int_t timeBin)
{
  for(auto &aRow : mRows) {
    if(aRow == nullptr) continue;
    aRow->processCommonMode(cru, timeBin, aRow->getRow());
    mTotalChargeTimeBin += aRow->getTotalChargeRow();
  }
}
