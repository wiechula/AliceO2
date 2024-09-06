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
#include "CCDB/CcdbApi.h"
#include "TPCBase/CommonModeCorrection.h"
#include "TPCBase/Mapper.h"

using namespace o2::tpc;
float CommonModeCorrection::getCommonMode(gsl::span<const float> values, int cru) const
{
  std::mt19937 rng(std::time(nullptr));
  std::uniform_int_distribution<std::mt19937::result_type> dist(0, values.size() - 1);

  std::vector<float> adcCM; //< ADC values used for common mode calculation

  for (size_t iPad = 0; iPad < values.size(); ++iPad) {
    const float kPulser = getCalPadValue("CMkValues", cru, iPad);
    const float pedestal = getCalPadValue("Pedestals", cru, iPad);
    const float adcPad = values[iPad] - pedestal;
    const float adcPadNorm = adcPad / kPulser;

    if (adcPad > mQEmpty) {
      continue;
    }

    int nPadsOK = 0;

    for (int iRnd = 0; iRnd < mNPadsCompRamdom; ++iRnd) {
      int padRnd = dist(rng);
      while (padRnd == iPad) {
        padRnd = dist(rng);
      }
      const float kPulserRnd = getCalPadValue("CMkValues", cru, padRnd);
      const float pedestalRnd = getCalPadValue("Pedestals", cru, padRnd);
      const float adcPadRnd = values[padRnd] - pedestalRnd;
      if (std::abs(adcPadNorm - adcPadRnd / kPulserRnd) < mQComp) {
        ++nPadsOK;
      }
    }
    if (nPadsOK >= mNPadsCompMin) {
      adcCM.emplace_back(adcPadNorm);
    }
  }
  const float entriesCM = adcCM.size();
  const float commonMode = (entriesCM > 0) ? std::accumulate(adcCM.begin(), adcCM.end(), 0.f) / entriesCM : 0.f;
  return commonMode;
}

float CommonModeCorrection::getCommonMode(const std::vector<float>& values, int cru) const
{
  return getCommonMode(gsl::span(values), cru);
}

std::vector<float> CommonModeCorrection::getCommonMode(const std::vector<Digit>& digits, int timeBin)
{
  // we assume that have one value per pad and per time bin, ordered pad-by-pad, row-by-row, timebin-by-timebin
  // as it comes out of the raw to digits conversion

  std::vector<float> cmValues(10);
  const int nPadsSector = Mapper::getPadsInSector();

  // sanity checks
  const size_t firstDigitOffset = nPadsSector * timeBin;
  const size_t lastDigitOffset = nPadsSector * (timeBin + 1) - 1;
  if (firstDigitOffset >= digits.size()) {
    LOGP(error, "Time bin {} out of range: digits size: {}", timeBin, digits.size());
    return cmValues;
  }

  const auto& firstDigit = digits.at(firstDigitOffset);
  const auto firstDigitTB = firstDigit.getTimeStamp();

  if (firstDigitTB != timeBin) {
    LOGP(error, "expected digit time stamp {} does not match requested {}", firstDigitTB, timeBin);
    return cmValues;
  }

  if (lastDigitOffset < digits.size()) {
    const auto lastDigitTB = digits.at(lastDigitOffset).getTimeStamp();
    const auto nextDigitTB = digits.at(lastDigitOffset + 1).getTimeStamp();
    if (lastDigitTB != timeBin || nextDigitTB != (timeBin + 1)) {
      LOGP(error, "Time stamp mismatch. Requested {}, last digit {}, next digit {}", timeBin, lastDigitTB, nextDigitTB);
      return cmValues;
    }
  } else {
    if ((lastDigitOffset + 1) != digits.size()) {
      LOGP(error, "End of digits expected, but next digit does not exactly point to the end: {} != {}", lastDigitOffset + 1, digits.size());
      return cmValues;
    }
  }

  for (int icru = 0; icru < CRU::MaxCRU; ++icru) {
    const int nPads = Mapper::PADSPERREGION[icru];
    const int padOffset = Mapper::GLOBALPADOFFSET[icru];
    const int cruSector = digits[firstDigitOffset + padOffset].getCRU();

    std::vector<float> values(nPads);

    for (int iPad = 0; iPad < nPads; ++iPad) {
      values[iPad] = digits[firstDigitOffset + padOffset + iPad].getChargeFloat();
    }

    const auto cm = getCommonMode(values, cruSector);
    cmValues[icru] = cm;
  }
  return cmValues;
}

void CommonModeCorrection::applyCommonModeCorrection(std::vector<Digit>& digits, bool subtractPedestal)
{
  const int lastTimeBin = digits[digits.size() - 1].getTimeStamp();

  size_t iDigit = 0;
  for (int iTimeBin = 0; iTimeBin <= lastTimeBin; ++iTimeBin) {
    const auto cmValues = getCommonMode(digits, iTimeBin);
    for (; iDigit < digits.size(); ++iDigit) {
      auto& digit = digits[iDigit];
      if (digit.getTimeStamp() < iTimeBin) {
        continue;
      }
      if (digit.getTimeStamp() > iTimeBin) {
        break;
      }
      const CRU cru(digit.getCRU());
      float charge = digit.getChargeFloat() - cmValues[cru.region()];
      if (subtractPedestal) {
        const float pedestal = mPadMaps["Pedestals"].getValue(cru.sector(), digit.getRow(), digit.getPad());
        charge -= pedestal;
      }
      digit.setCharge(charge);
    }
  }
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
