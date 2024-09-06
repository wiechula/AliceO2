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

/// \file CommonModeCorrection.h
/// \brief  Calculate the common mode correction factor
/// \author Jens Wiechula, Jens.Wiechula@ikf.uni-frankfurt.de

#ifndef AliceO2_TPC_CommonModeCorrection_H_
#define AliceO2_TPC_CommonModeCorrection_H_

#include <gsl/span>
#include <vector>

#include <DataFormatsTPC/Digit.h>
#include <TPCBase/FEEConfig.h>
#include <TPCBase/CRU.h>

namespace o2::tpc
{

/// Class to calculate the common mode correction
///
/// Calculation of the common mode correction, based on the algorithm propsed by Marian Ivanov
/// The calculation is done for one single CRU and time bin
class CommonModeCorrection
{
 public:
  using CalPadMapType = std::unordered_map<std::string, CalPad>;

  /// Calculation of the common mode value
  ///
  /// \param value Charge values ordered per row and pad, the order is mendatory
  /// \param cru cru the charge values belong to
  float getCommonMode(gsl::span<const float> values, int cru) const;
  float getCommonMode(const std::vector<float>& values, int cru) const;

  /// Calculate common mode for all CRUs in one sector
  /// \param digits raw digits of one sector, not zero suppressed
  /// \param timeBin time bin to analyse
  /// \return vector of size 10 containing the CM of all CRUs
  std::vector<float> getCommonMode(const std::vector<Digit>& digits, int timeBin);

  /// Apply the common mode correction to all values of one sector
  /// \param digits raw digits of one sector, not zero suppressed
  void applyCommonModeCorrection(std::vector<Digit>& digits, bool subtractPedestal = false);

  void setNPadsCompRandom(int n) { mNPadsCompRamdom = n; }
  int getNPadsCompRandom() const { return mNPadsCompRamdom; }

  void setNPadsCompMin(int n) { mNPadsCompRamdom = n; }
  int getNPadsCompMin() const { return mNPadsCompRamdom; }

  void setQEmpty(float q) { mQEmpty = q; }
  float getQEmpty() const { return mQEmpty; }

  void setQComp(float q) { mQComp = q; }
  float getQComp() const { return mQComp; }

  void setITCorr(float corr) { mITCorr = corr; }
  float getITCorr() const { return mITCorr; }

  /// Pad maps loaded from FEEConfig
  void setPadMaps(CalPadMapType& padMaps) { mPadMaps = padMaps; }

  /// Custom pedestals, overwriting what was set in the padMaps
  void setPedestals(const CalPad& pedestals) { mPadMaps["Pedestals"] = pedestals; }

  void loadDefaultPadMaps(FEEConfig::Tags feeTag = FEEConfig::Tags::Physics30sigma);

 private:
  int mNPadsCompRamdom{10}; /// Number of random pads to compare with to check if the present pad is empty
  int mNPadsCompMin{7};     /// Minimum number of neighbouring pads with q close to present pad to define this as empty
  float mQEmpty{2};         ///< Threshold to enter check for empty pad
  float mQComp{1};          ///< Threshold for comparison with random pads
  float mITCorr{1.f};       ///< Ion tail scaling value
  CalPadMapType mPadMaps;   ///< Pad-by-pad CRU configuration values (Pedestal, Noise, ITF + CM parameters)

  /// Return the value stored in mPadMaps["calibName"]
  /// \param calibName name of calibraion in mPadMaps
  /// \param cru CRU number
  /// \param pad Pad number within the CRU
  float getCalPadValue(const std::string calibName, int icru, int pad) const;

  ClassDefNV(CommonModeCorrection, 0);
};

} // namespace o2::tpc
#endif
