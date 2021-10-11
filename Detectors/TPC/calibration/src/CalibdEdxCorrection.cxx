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

#include "TPCCalibration/CalibdEdxCorrection.h"

#include <algorithm>
#include <array>
#include <cstddef>

// o2 includes
#include "DataFormatsTPC/Defs.h"

using namespace o2::tpc;

size_t CalibdEdxCorrection::stackIndex(const StackID& stack, ChargeType charge)
{
  auto index = static_cast<size_t>(stack.sector);
  index += static_cast<size_t>(stack.side * SECTORSPERSIDE);
  index += static_cast<size_t>(stack.type * SECTORSPERSIDE * SIDES);
  index += static_cast<size_t>(charge * SECTORSPERSIDE * SIDES * GEMSTACKSPERSECTOR);
  return index;
}

float CalibdEdxCorrection::correction(const StackID& stack, ChargeType charge, float tgl, float snp) const
{
  const auto& p = getParams(stack, charge);
  float corr = p[0];
  if (mDims > 0) {
    corr += p[1] * tgl + p[2] * tgl * tgl;
    if (mDims > 1) {
      corr += p[3] * snp + p[4] * tgl * snp + p[5] * snp * snp;
    }
  }

  return corr;
}

void CalibdEdxCorrection::clear()
{
  for (auto& x : mParams) {
    std::fill(x.begin(), x.end(), 0.0);
  }
  std::fill(mChi2.begin(), mChi2.end(), 0.0);
  mDims = 0;
}
