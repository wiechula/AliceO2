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
#include <string_view>
#include <vector>

#include "DataFormatsTPC/Digit.h"
#include "TPCBase/FEEConfig.h"

namespace o2::tpc
{

/// Class to calculate the common mode correction
///
/// Calculation of the common mode correction, based on the algorithm propsed by Marian Ivanov
/// The calculation is done for one single CRU and time bin
class CommonModeCorrection
{
 public:
  struct CMdata {
    std::vector<float> adcValues;
    std::vector<float> cmKValues;
    std::vector<float> pedestals;

    void clear()
    {
      adcValues.clear();
      cmKValues.clear();
      pedestals.clear();
    }
  };

  struct CMInfo {
    float cmValue{};
    int nPadsUsed{};
  };

  using CalPadMapType = std::unordered_map<std::string, CalPad>;

  /// Calculation of the common mode value
  ///
  /// \param value pad-by-pad charge values
  /// \param cmKValues corresponding pad-by-pad common mode k-factors
  /// \param pedestals corresponding pad-by-pad pedestals
  /// \param
  CMInfo getCommonMode(gsl::span<const float> values, gsl::span<const float> cmKValues, gsl::span<const float> pedestals) const;
  CMInfo getCommonMode(const std::vector<float>& values, const std::vector<float>& cmKValues, const std::vector<float>& pedestals) const { return getCommonMode(gsl::span(values), gsl::span(cmKValues), gsl::span(pedestals)); }

  CMInfo getCommonMode(const CMdata& cmData) const { return getCommonMode(std::span(cmData.adcValues), std::span(cmData.cmKValues), std::span(cmData.pedestals)); }

  void setNPadsCompRandom(int n) { mNPadsCompRamdom = n; }
  int getNPadsCompRandom() const { return mNPadsCompRamdom; }

  void setNPadsCompMin(int n) { mNPadsCompRamdom = n; }
  int getNPadsCompMin() const { return mNPadsCompRamdom; }

  void setQEmpty(float q) { mQEmpty = q; }
  float getQEmpty() const { return mQEmpty; }

  void setQComp(float q) { mQComp = q; }
  float getQComp() const { return mQComp; }

  /// Pad maps loaded from FEEConfig
  void setPadMaps(CalPadMapType& padMaps) { mPadMaps = padMaps; }

  /// load a CalPad from file and add it to the local mPadMaps
  /// \param fileName input file name
  /// \param nameInFile name of the CalPad object in the file
  /// \param namePadMap name under which to store the object in the mPadMaps, if empty use the same as nameInFile
  void loadCalPad(std::string_view fileName, std::string_view nameInFile, std::string_view namePadMap = "");

  /// load CMkValues from file, assuming it is stored under the name "CMkValues
  void loadCMkValues(std::string_view fileName) { loadCalPad(fileName, "CMkValues"); }

  /// load Pedestals from file, assuming it is stored under the name "Pedestals
  void loadPedestals(std::string_view fileName) { loadCalPad(fileName, "Pedestals"); }

  /// Custom setting of CalPad, overwriting what was set in mPadMaps
  void setCalPad(const CalPad& calPad, std::string_view name) { mPadMaps[name.data()] = calPad; }

  /// load the Pad maps from CCDB
  void loadDefaultPadMaps(FEEConfig::Tags feeTag = FEEConfig::Tags::Physics30sigma);

  CMdata collectCMdata(const std::vector<Digit>& digits, int cru, int timeBin);

  /// corret digits for common mode
  /// \param cmValues will contain CM information for each CRU and time bin
  /// \param negativeOnly only correct negative common mode signals
  /// \return maximum
  int correctDigits(std::vector<Digit>& digits, std::vector<std::vector<CMInfo>>& cmValues, bool negativeOnly = false) const;

  void correctDigits(std::string_view digiFileIn, Long64_t maxEntries = -1, std::string_view digitFileOut = "tpcdigit_cmcorr.root", std::string_view cmFileOut = "CommonModeValues.root", bool negativeOnly = false, int nThreads = 1, bool writeOnlyCM = false);

  void limitKFactorPrecision(bool limit = true) { mLimitKFactor = limit; }
  void limitPedestalPrecision(bool limit = true) { mLimitPedestal = limit; }

  /// set the number of threads used for CM calculation
  /// \param nThreads number of threads
  static void setNThreads(const int nThreads) { sNThreads = nThreads; }

  /// \return returns the number of threads used for decoding
  static int getNThreads() { return sNThreads; }

 private:
  inline static int sNThreads{1}; /// Number of parallel threads for the CM calculation
  int mNPadsCompRamdom{10};       /// Number of random pads to compare with to check if the present pad is empty
  int mNPadsCompMin{7};           /// Minimum number of neighbouring pads with q close to present pad to define this as empty
  float mQEmpty{2};               ///< Threshold to enter check for empty pad
  float mQComp{1};                ///< Threshold for comparison with random pads
  bool mLimitKFactor{false};      ///< Limit the k-factor precision to 2I6F
  bool mLimitPedestal{false};     ///< Limit the preestal precision to 10I2F

  CalPadMapType mPadMaps; ///< Pad-by-pad CRU configuration values (Pedestal, Noise, ITF + CM parameters)

  /// Return the value stored in mPadMaps["calibName"]
  /// \param calibName name of calibraion in mPadMaps
  /// \param cru CRU number
  /// \param pad Pad number within the CRU
  float getCalPadValue(const std::string calibName, int icru, int pad) const;

  bool padMapExists(const std::string& calibName);

  ClassDefNV(CommonModeCorrection, 0);
};

} // namespace o2::tpc
#endif
