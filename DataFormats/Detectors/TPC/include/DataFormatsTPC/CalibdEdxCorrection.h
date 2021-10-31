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

/// \file CalibdEdxCorrection.h
/// \author Thiago Badaró <thiago.saramela@usp.br>

#ifndef ALICEO2_TPC_CALIBDEDXCORRECTION_H_
#define ALICEO2_TPC_CALIBDEDXCORRECTION_H_

#include <cstddef>
#include <array>
#include <string_view>

// o2 includes
#include "DataFormatsTPC/Defs.h"

namespace o2::tpc
{

class CalibdEdxCorrection
{
 public:
  using Params = std::array<float, 6>;
  CalibdEdxCorrection() { clear(); }
  CalibdEdxCorrection(std::string_view fileName) { loadFile(fileName); }

  float getCorrection(const StackID&, ChargeType, float z = 0, float tgl = 0) const;

  const Params& getParams(const StackID& stack, ChargeType charge) const { return mParams[stackIndex(stack, charge)]; }
  float getChi2(const StackID& stack, ChargeType charge) const { return mChi2[stackIndex(stack, charge)]; }
  int getDims() const { return mDims; }

  void setParams(const StackID& stack, ChargeType charge, const Params& params) { mParams[stackIndex(stack, charge)] = params; }
  void setChi2(const StackID& stack, ChargeType charge, float chi2) { mChi2[stackIndex(stack, charge)] = chi2; }
  void setDims(int dims) { mDims = dims; }

  void clear();

  void saveFile(std::string_view fileName) const;
  void loadFile(std::string_view fileName);

 private:
  static size_t stackIndex(const StackID& stack, ChargeType charge)
  {
    return static_cast<size_t>(stack.index() + charge * SECTORSPERSIDE * SIDES * GEMSTACKSPERSECTOR);
  }

  std::array<Params, 288> mParams{};
  std::array<float, 288> mChi2{};
  int mDims{}; ///< Fit dimension
};

} // namespace o2::tpc

#endif
