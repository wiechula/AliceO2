// Copyright CERN and copyright holders of ALICE O3. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file   CalibPedestal.cxx
/// \author Jens Wiechula, Jens.Wiechula@ikf.uni-frankfurt.de

#include "TFile.h"
#include "TObjArray.h"
#include "TH2S.h"
#include "TH1F.h"
#include "TF1.h"

#include "TPCBase/ROC.h"
#include "MathUtils/MathBase.h"
#include "TPCCalibration/CalibPedestal.h"

using namespace o2::TPC;
using o2::mathUtils::mathBase::fitGaus;

CalibPedestal::CalibPedestal(PadSubset padSubset)
  : CalibRawBase(padSubset),
    mFirstTimeBin(0),
    mLastTimeBin(500),
    mADCMin(0),
    mADCMax(120),
    mNumberOfADCs(mADCMax - mADCMin + 1),
    mUseSlowFit(false),
    mWriteDebugHistos(false),
    mPedestal(padSubset),
    mNoise(padSubset),
    mADCdata()

{
  mADCdata.resize(ROC::MaxROC);
  mPedestal.setName("Pedestals");
  mNoise.setName("Noise");
}

//______________________________________________________________________________
Int_t CalibPedestal::updateROC(const Int_t roc, const Int_t row, const Int_t pad, const Int_t timeBin,
                               const Float_t signal)
{
  Int_t adcValue = Int_t(signal);
  if (timeBin < mFirstTimeBin || timeBin > mLastTimeBin) {
    return 0;
  }

  if (adcValue < mADCMin || adcValue > mADCMax) {
    return 0;
  }

  const GlobalPadNumber padInROC = mMapper.getPadNumberInROC(PadROCPos(roc, row, pad));
  Int_t bin = padInROC * mNumberOfADCs + (adcValue - mADCMin);
  vectorType& adcVec = *getVector(ROC(roc), kTRUE);
  ++(adcVec[bin]);

  // printf("bin: %5d, val: %.2f\n", bin, adcVec[bin]);

  return 0;
}

//______________________________________________________________________________
CalibPedestal::vectorType* CalibPedestal::getVector(ROC roc, bool create /*=kFALSE*/)
{
  vectorType* vec = mADCdata[roc].get();
  if (vec || !create) {
    return vec;
  }

  const size_t numberOfPads = (roc.rocType() == RocType::IROC) ? mMapper.getPadsInIROC() : mMapper.getPadsInOROC();

  vec = new vectorType;
  vec->resize(numberOfPads * mNumberOfADCs);

  mADCdata[roc] = std::unique_ptr<vectorType>(vec);

  return vec;
}

//______________________________________________________________________________
void CalibPedestal::analyse()
{
  ROC roc;

  std::vector<float> fitValues;

  TH1F hDummy("hDummy", "hDummy", mNumberOfADCs - 2, mADCMin + 1, mADCMax - 1);
  TF1 fgaus("fgaus", "gaus");

  for (auto& vecPtr : mADCdata) {
    auto vec = vecPtr.get();
    if (!vec) {
      ++roc;
      continue;
    }

    CalROC& calROCPedestal = mPedestal.getCalArray(roc);
    CalROC& calROCNoise = mNoise.getCalArray(roc);

    float* array = vec->data();

    const size_t numberOfPads = (roc.rocType() == RocType::IROC) ? mMapper.getPadsInIROC() : mMapper.getPadsInOROC();

    for (Int_t ichannel = 0; ichannel < numberOfPads; ++ichannel) {
      size_t offset = ichannel * mNumberOfADCs;

      if (mUseSlowFit) {
        hDummy.Set(mNumberOfADCs, array + offset);
        hDummy.Fit(&fgaus, "QN");
        const float pedestal = fgaus.GetParameter(1);
        const float noise = fgaus.GetParameter(2);
        calROCPedestal.setValue(ichannel, pedestal);
        calROCNoise.setValue(ichannel, noise);
        fgaus.SetParameters(0, 0, 0);
      } else {
        fitGaus(mNumberOfADCs, array + offset, float(mADCMin), float(mADCMax + 1), fitValues);

        const float pedestal = fitValues[1];
        const float noise = fitValues[2];
        // printf("roc: %2d, channel: %4d, pedestal: %.2f, noise: %.2f\n", roc.getRoc(), ichannel, pedestal, noise);
      }
    }

    ++roc;
  }
}

//______________________________________________________________________________
void CalibPedestal::resetData()
{
  for (auto& vecPtr : mADCdata) {
    auto vec = vecPtr.get();
    if (!vec) {
      continue;
    }
    vec->clear();
  }
}

//______________________________________________________________________________
void CalibPedestal::dumpToFile(const std::string filename)
{
  auto f = std::unique_ptr<TFile>(TFile::Open(filename.data(), "recreate"));
  f->WriteObject(&mPedestal, "Pedestals");
  f->WriteObject(&mNoise, "Noise");

  // debug data
  if (mWriteDebugHistos) {
    TObjArray arrHist;
    arrHist.SetOwner();

    for (int iroc = 0; iroc < mADCdata.size(); ++iroc) {
      auto vec = mADCdata[iroc].get();
      if (!vec)
        continue;

      const size_t numberOfPads =
        (ROC(iroc).rocType() == RocType::IROC) ? mMapper.getPadsInIROC() : mMapper.getPadsInOROC();
      TH2S* h = new TH2S(Form("hROC_02d", iroc), Form("Baseline distribution ROC %2d;ADC value;channel", iroc),
                         mNumberOfADCs, 0, mNumberOfADCs, numberOfPads, 0, numberOfPads);
      arrHist.AddAt(h, iroc);
      for (int ichannel = 0; ichannel < numberOfPads; ++ichannel) {
        for (int iadc = 0; iadc < mNumberOfADCs; ++iadc) {
          size_t offset = ichannel * mNumberOfADCs;
          h->SetBinContent(iadc + 1, ichannel + 1, vec->data()[offset + iadc]);
        }
      }
    }

    arrHist.Write("debug_histos", TObject::kSingleKey);
  }

  f->Close();
}
