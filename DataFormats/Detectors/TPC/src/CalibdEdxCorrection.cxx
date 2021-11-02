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

#include "DataFormatsTPC/CalibdEdxCorrection.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <string_view>

// o2 includes
#include "DataFormatsTPC/Defs.h"

// root includes
#include "TFile.h"

using namespace o2::tpc;

GPUd() float CalibdEdxCorrection::getCorrection(const StackID& stack, ChargeType charge, float z, float tgl) const
{
  const auto& p = getParams(stack, charge);
  float corr = p[0];
  if (mDims > 0) {
    corr += p[1] * z + p[2] * z * z;
    if (mDims > 1) {
      corr += p[3] * tgl + p[4] * z * tgl + p[5] * tgl * tgl;
    }
  }

  return corr;
}

#if !defined(GPUCA_GPUCODE)
void CalibdEdxCorrection::clear()
{
  for (auto& x : mParams) {
    std::fill(x.begin(), x.end(), 0.0);
  }
  std::fill(mChi2.begin(), mChi2.end(), 0.0);
  mDims = std::numeric_limits<int>::min();
}

void CalibdEdxCorrection::saveFile(std::string_view fileName) const
{
  std::unique_ptr<TFile> file(TFile::Open(fileName.data(), "recreate"));
  file->WriteObject(this, "CalibdEdxCorrection");
}

void CalibdEdxCorrection::loadFile(std::string_view fileName)
{
  std::unique_ptr<TFile> file(TFile::Open(fileName.data()));
  auto tmp = file->Get<CalibdEdxCorrection>("CalibdEdxCorrection");
  if (tmp != nullptr) {
    *this = *tmp;
  }
}
#endif
