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

/// \file CalibdEdx.h
/// \brief This file provides the container used for time based dE/dx calibration.
/// \author Thiago Badaró <thiago.saramela@usp.br>

#ifndef ALICEO2_TPC_CALIBDEDX_H_
#define ALICEO2_TPC_CALIBDEDX_H_

#include <cstddef>
#include <gsl/span>
#include <string_view>
#include <tuple>

// o2 includes
#include "DataFormatsTPC/TrackCuts.h"
#include "DataFormatsTPC/Defs.h"

// boost includes
#include <boost/histogram.hpp>

// root includes
#include "TH2F.h"

namespace o2::tpc
{

// forward declaration
class TrackTPC;

struct CalibdEdxData {
  float mean{};
  float stdm{};
  int sector{};
  Side side{};
  GEMstack stack{};
  float tgl{};
  float snp{};
  ChargeType charge{};
};

/// Class that creates dE/dx histograms from a sequence of tracks objects
class CalibdEdx
{
 public:
  enum HistAxis {
    dEdx = 0,
    Sector = 1,
    Side = 2,
    Stack = 3,
    Tgl = 4,
    Snp = 5,
    Charge = 6,
    Size = 7 ///< Number of axes
  };

  // Interger histogram axis identifying the GEM stacks, without under and overflow bins.
  using HistIntAxis = boost::histogram::axis::integer<int, boost::histogram::use_default, boost::histogram::axis::option::none_t>;
  // Float axis to store data, without under and overflow bins.
  using HistFloatAxis = boost::histogram::axis::regular<float, boost::histogram::use_default, boost::histogram::use_default, boost::histogram::axis::option::none_t>;
  //  using HistFloatAxis = boost::histogram::axis::regular<>;

  // Define histogram axes types
  using HistAxesType = std::tuple<
    HistFloatAxis, // dEdx
    HistIntAxis,   // sector
    HistIntAxis,   // side
    HistIntAxis,   // stack type
    HistFloatAxis, // Tgl
    HistFloatAxis, // Snp
    HistIntAxis    // Charge
    >;

  using Hist = boost::histogram::histogram<HistAxesType>;
  using CalibContainer = std::vector<CalibdEdxData>;

  /// Default constructor
  CalibdEdx() = default;

  /// Constructor that enable tracks cuts.
  CalibdEdx(int nBins, float mindEdx, float maxdEdx, const TrackCuts& cuts);

  /// Constructor that enable tracks cuts, and creates a TrackCuts internally.
  CalibdEdx(int nBins, float mindEdx = 5, float maxdEdx = 70, float minP = 0.4, float maxP = 0.6, int minClusters = 60)
    : CalibdEdx(nBins, mindEdx, maxdEdx, {minP, maxP, static_cast<float>(minClusters)}) {}

  /// Fill histograms using tracks data.
  void fill(const gsl::span<const TrackTPC>);
  void fill(const std::vector<TrackTPC>&);
  void fill(const TrackTPC&);

  /// Add counts from other container.
  void merge(const CalibdEdx* other);

  /// Compute MIP position from dEdx histograms, and save result in the calib container.
  void finalize();

  /// Return the full, unprojected, histogram.
  const Hist& getHist() const { return mHist; }
  /// Return the projected histogram, the projected axes are summed over.
  auto getHist(const std::vector<int>& projected_axis) const
  {
    return boost::histogram::algorithm::project(mHist, projected_axis);
  }

  /// Return the projected histogram as a TH2F, the unkept axis are summed over.
  TH2F getRootHist(const std::vector<int>& projected_axis) const;
  /// Keep all axes
  TH2F getRootHist() const;

  const CalibContainer& getCalib() const { return mCalib; }

  /// \brief Check if there are enough data to compute the calibration.
  /// \param minEntries in each histogram
  /// \return false if any of the histograms has less entries than minEntries
  bool hasEnoughData(float minEntries) const;

  void setApplyCuts(bool apply) { mApplyCuts = apply; }
  bool getApplyCuts() const { return mApplyCuts; }
  void setCuts(const TrackCuts& cuts) { mCuts = cuts; }

  /// Find the sector of a track at a given GEM stack type.
  static int findTrackSector(const TrackTPC& track, GEMstack, bool& ok);

  /// Print the number of entries in each histogram.
  void print() const;

  /// Save the histograms to a file.
  void dumpToFile(std::string_view fileName) const;

  void writeTTree(std::string_view fileName) const;

 private:
  bool mApplyCuts{true}; ///< Whether or not to apply tracks cuts
  int mNBins;            ///< Number of dEdx bins
  TrackCuts mCuts;       ///< Cut class

  Hist mHist;              ///< TotdEdx multidimensional histogram
  CalibContainer mCalib{}; ///< Calibration output container

  ClassDefNV(CalibdEdx, 1);
};

} // namespace o2::tpc

#endif
