// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file   CommonModeCorrection.cxx
/// \brief  Calculate the common mode correction factor

#include <random>
#include <algorithm>
#include <thread>
#include <mutex>
#include "CCDB/CcdbApi.h"
#include "TPCBase/CommonModeCorrection.h"
#include "TPCBase/Mapper.h"
#include "TPCBase/Utils.h"
#include "TPCBase/CRUCalibHelpers.h"
#include "TChain.h"
#include "TFile.h"
#include "MathUtils/RandomRing.h"
#include "CommonUtils/TreeStreamRedirector.h"

using namespace o2::tpc;
using namespace o2::tpc::cru_calib_helpers;
CommonModeCorrection::CMInfo CommonModeCorrection::getCommonMode(gsl::span<const float> values, gsl::span<const float> cmKValues, gsl::span<const float> pedestals) const
{
  if (values.size() == 0) {
    return CMInfo{};
  }
  // sanity check
  if (values.size() != cmKValues.size() || values.size() != pedestals.size()) {
    LOGP(error, "vector sizes of input values, cmKValues and pedestals don't match: {}, {}, {}", values.size(), cmKValues.size(), pedestals.size());
    return CMInfo{};
  }
  std::mt19937 rng(std::time(nullptr));
  std::uniform_int_distribution<std::mt19937::result_type> dist(0, values.size() - 1);

  static math_utils::RandomRing random(math_utils::RandomRing<>::RandomType::Flat);
  std::vector<std::vector<float>> adcCM(sNThreads); //< ADC values used for common mode calculation

  auto thread = [&](int iThread) {
    for (size_t iPad = iThread; iPad < values.size(); iPad += sNThreads) {
      const float kCM = mLimitKFactor ? fixedSizeToFloat<6>(floatToFixedSize<8, 6>(cmKValues[iPad])) : cmKValues[iPad];
      const float pedestal = mLimitPedestal ? fixedSizeToFloat(floatToFixedSize(pedestals[iPad])) : pedestals[iPad];
      const float adcPad = values[iPad] - pedestal;
      const float adcPadNorm = adcPad / kCM;

      if (adcPad > mQEmpty) {
        continue;
      }

      int nPadsOK = 0;

      for (int iRnd = 0; iRnd < mNPadsCompRamdom; ++iRnd) {
        // int padRnd = dist(rng);
        int padRnd = 0;
        do {
          // padRnd = dist(rng);
          // static std::mutex rand_mutex;
          // std::lock_guard lock{rand_mutex};
          padRnd = int(random.getNextValue() * (values.size() - 1));
        } while (padRnd == iPad);
        const float kCMRnd = mLimitKFactor ? fixedSizeToFloat<6>(floatToFixedSize<8, 6>(cmKValues[padRnd])) : cmKValues[padRnd];
        const float pedestalRnd = mLimitPedestal ? fixedSizeToFloat(floatToFixedSize(pedestals[padRnd])) : pedestals[padRnd];
        const float adcPadRnd = values[padRnd] - pedestalRnd;
        if (std::abs(adcPadNorm - adcPadRnd / kCMRnd) < mQComp) {
          ++nPadsOK;
        }
      }

      if (nPadsOK >= mNPadsCompMin) {
        adcCM[iThread].emplace_back(adcPadNorm);
      }
    }
  };

  // start threads
  std::vector<std::thread> threads(sNThreads);

  for (int i = 0; i < sNThreads; i++) {
    threads[i] = std::thread(thread, i);
  }

  // wait for the threads to finish
  for (auto& th : threads) {
    th.join();
  }

  int entriesCM = 0;
  float commonMode = 0.f;
  for (const auto& adcs : adcCM) {
    entriesCM += int(adcs.size());
    commonMode += std::accumulate(adcs.begin(), adcs.end(), 0.f);
  }

  if (entriesCM > 0) {
    commonMode /= float(entriesCM);
  }
  return CMInfo{commonMode, entriesCM};
}

void CommonModeCorrection::loadDefaultPadMaps(FEEConfig::Tags tag)
{
  o2::ccdb::CcdbApi cdbApi;
  cdbApi.init("http://alice-ccdb.cern.ch");
  const auto feeConfig = cdbApi.retrieveFromTFileAny<FEEConfig>("TPC/Config/FEE", {}, long(tag));
  if (!feeConfig) {
    LOGP(error, "Could not retrieve pad maps");
    return;
  }
  mPadMaps = feeConfig->padMaps;
  delete feeConfig;
}

CommonModeCorrection::CMdata CommonModeCorrection::collectCMdata(const std::vector<Digit>& digits, int cru, int timeBin)
{

  CMdata data;
  if (!padMapExists("CMkValues") || padMapExists("Pedestals")) {
    return data;
  }

  for (const auto& digit : digits) {
    if (digit.getTimeStamp() < timeBin) {
      continue;
    }

    if (digit.getTimeStamp() > timeBin) {
      break;
    }

    if (digit.getCRU() < cru) {
      continue;
    }

    if (digit.getCRU() > cru) {
      break;
    }

    const auto sector = CRU(digit.getCRU()).sector();
    data.adcValues.emplace_back(digit.getChargeFloat());
    data.cmKValues.emplace_back(mPadMaps["CMkValues"].getValue(sector, digit.getRow(), digit.getPad()));
    data.pedestals.emplace_back(mPadMaps["Pedestals"].getValue(sector, digit.getRow(), digit.getPad()));
  }
  return data;
}

int CommonModeCorrection::correctDigits(std::vector<Digit>& digits, std::vector<std::vector<CMInfo>>& cmValues, bool negativeOnly) const
{
  // calculation common mode values
  int maxTimeBin = -1;
  int lastCRU = -1;
  int lastTimeBin = -1;
  CMdata data;
  const auto& cmkValues = mPadMaps.at("CMkValues");
  const auto& pedestals = mPadMaps.at("Pedestals");

  for (size_t iDigit = 0; iDigit < digits.size(); ++iDigit) {
    auto& digit = digits[iDigit];
    if (((lastCRU > 0) && (digit.getCRU() != lastCRU)) || ((lastTimeBin > 0) && (digit.getTimeStamp() != lastTimeBin))) {
      auto& cmValuesCRU = cmValues[lastCRU];
      if (cmValuesCRU.size() <= lastTimeBin) {
        cmValuesCRU.resize(cmValuesCRU.size() + 500);
      }
      cmValuesCRU[lastTimeBin] = getCommonMode(data.adcValues, data.cmKValues, data.pedestals);
      // LOGP(info, "processing CRU {}, timeBin {}, CM = {}", lastCRU, lastTimeBin, cmVal);

      data.clear();
    }
    const auto sector = CRU(digit.getCRU()).sector();
    data.adcValues.emplace_back(digit.getChargeFloat());
    data.cmKValues.emplace_back(cmkValues.getValue(sector, digit.getRow(), digit.getPad()));
    data.pedestals.emplace_back(pedestals.getValue(sector, digit.getRow(), digit.getPad()));

    lastCRU = digit.getCRU();
    lastTimeBin = digit.getTimeStamp();
    maxTimeBin = std::max(lastTimeBin, maxTimeBin);
  }
  {
    auto& cmValuesCRU = cmValues[lastCRU];
    if (cmValuesCRU.size() <= lastTimeBin) {
      cmValuesCRU.resize(cmValuesCRU.size() + 500);
    }
    cmValuesCRU[lastTimeBin] = getCommonMode(data.adcValues, data.cmKValues, data.pedestals);
    // LOGP(info, "processing CRU {}, timeBin {}, CM = {}", lastCRU, lastTimeBin, cmVal);
    data.clear();
  }

  // apply correction
  for (auto& digit : digits) {
    const auto sector = CRU(digit.getCRU()).sector();
    const auto cmKValue = cmkValues.getValue(sector, digit.getRow(), digit.getPad());
    const auto cmValue = cmValues[digit.getCRU()][digit.getTimeStamp()].cmValue;
    if (!negativeOnly || cmValue < 0) {
      digit.setCharge(digit.getCharge() - cmValue * cmKValue);
    }
  }

  return maxTimeBin;
}

void CommonModeCorrection::correctDigits(std::string_view digiFileIn, Long64_t maxEntries, std::string_view digitFileOut, bool negativeOnly)
{

  TChain* tree = o2::tpc::utils::buildChain(fmt::format("ls {}", digiFileIn), "o2sim", "o2sim");
  Long64_t nEntries = tree->GetEntries();
  if (maxEntries > 0) {
    nEntries = std::min(nEntries, maxEntries);
  }

  std::unique_ptr<TFile> fOut{TFile::Open(digitFileOut.data(), "RECREATE")};
  TTree tOut("o2sim", "o2sim");

  std::array<std::vector<o2::tpc::Digit>*, 36> digitizedSignal;
  for (size_t iSec = 0; iSec < digitizedSignal.size(); ++iSec) {
    digitizedSignal[iSec] = nullptr;
    tree->SetBranchAddress(Form("TPCDigit_%zu", iSec), &digitizedSignal[iSec]);
    tOut.Branch(Form("TPCDigit_%zu", iSec), &digitizedSignal[iSec]);
  }

  o2::utils::TreeStreamRedirector pcstream("CommonModeValues.root", "recreate");

  for (Long64_t iTF = 0; iTF < nEntries; ++iTF) {
    tree->GetEntry(iTF);
    LOGP(info, "Processing entry {}/{}", iTF + 1, nEntries);

    std::vector<std::vector<CMInfo>> cmValues; // CRU * timeBin
    cmValues.resize(CRU::MaxCRU);
    int maxTimeBin = -1;

    for (size_t iSector = 0; iSector < 36; ++iSector) {
      auto digits = digitizedSignal[iSector];
      if (!digits || (digits->size() == 0)) {
        continue;
      }
      const int maxTimeBinSector = correctDigits(*digits, cmValues, negativeOnly);
      maxTimeBin = std::max(maxTimeBin, maxTimeBinSector);
    }
    tOut.Fill();

    for (int iCRU = 0; iCRU < cmValues.size(); ++iCRU) {
      int maxTBCRU = std::min(maxTimeBin, int(cmValues[iCRU].size()));
      for (int iTimeBin = 0; iTimeBin < maxTBCRU; ++iTimeBin) {
        pcstream << "cm"
                 << "iTF=" << iTF
                 << "iCRU=" << iCRU
                 << "iTimeBin=" << iTimeBin
                 << "cm=" << cmValues[iCRU][iTimeBin].cmValue
                 << "cmEntries=" << cmValues[iCRU][iTimeBin].nPadsUsed
                 << "\n";
      }
    }
  }

  pcstream.Close();
  fOut->cd();
  tOut.Write();
  fOut->Close();
}

float CommonModeCorrection::getCalPadValue(const std::string calibName, int icru, int pad) const
{
  if (mPadMaps.find(calibName) == mPadMaps.end()) {
    LOGP(error, "{} not set, cannot be used", calibName);
    return 0;
  }
  const auto& calPad = mPadMaps.at(calibName);
  const CRU cru(icru);
  const int roc = cru.roc();
  const int padOffset = (cru.isIROC()) ? Mapper::GLOBALPADOFFSET[cru.region()] : Mapper::GLOBALPADOFFSET[cru.region()] - Mapper::GLOBALPADOFFSET[4];

  const auto& calArray = calPad.getCalArray(roc);

  return calArray.getValue(padOffset + pad);
}

bool CommonModeCorrection::padMapExists(const std::string& calibName)
{
  if (mPadMaps.find(calibName) == mPadMaps.end()) {
    LOGP(error, "{} not in mPadMaps", calibName);
    return false;
  }
  return true;
}

void CommonModeCorrection::loadCalPad(std::string_view fileName, std::string_view nameInFile, std::string_view namePadMap)
{
  auto pads = o2::tpc::utils::readCalPads(fileName, nameInFile);
  if (pads.size() == 0) {
    LOGP(error, "Could not load object {} from file {}", nameInFile, fileName);
    return;
  }

  if (namePadMap.size() == 0) {
    namePadMap = nameInFile;
  }

  mPadMaps[namePadMap.data()] = *pads[0];
}
