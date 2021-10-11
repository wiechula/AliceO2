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

#include "TPCCalibration/CalibdEdx.h"

#include <algorithm>
#include <array>
#include <boost/histogram/indexed.hpp>
#include <cmath>
#include <cstddef>
#include <gsl/span>
#include <limits>
#include <numeric>
#include <string_view>
#include <utility>

// o2 includes
#include "DataFormatsTPC/TrackTPC.h"
#include "DataFormatsTPC/TrackCuts.h"
#include "DataFormatsTPC/Defs.h"
#include "Framework/Logger.h"
#include "MathUtils/fit.h"
#include "MathUtils/Utils.h"

// root includes
#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"

// boost includes
#include <boost/histogram.hpp>

using namespace o2::tpc;
namespace bh = boost::histogram;

CalibdEdx::CalibdEdx(int nBins, float mindEdx, float maxdEdx, const TrackCuts& cuts)
  : mCuts{cuts}, mNBins{nBins}
{
  mHist = bh::make_histogram(
    HistFloatAxis(nBins, mindEdx, maxdEdx, "dEdx"),
    HistFloatAxis(18, -1, 1, "Tgl"),
    HistFloatAxis(18, -1, 1, "Snp"),
    HistIntAxis(0, SECTORSPERSIDE, "sector"),
    HistIntAxis(0, SIDES, "side"),
    HistIntAxis(0, GEMSTACKSPERSECTOR, "stacktype"),
    HistIntAxis(0, CHARGETYPES, "charge"));
}

void CalibdEdx::fill(const TrackTPC& track)
{
  // applying cuts
  if (track.hasBothSidesClusters() || (mApplyCuts && !mCuts.goodTrack(track))) {
    return;
  }

  const auto& dEdx = track.getdEdx();
  const auto side = track.hasASideClustersOnly() ? Side::A : Side::C;
  const std::array<float, 4> dEdxTot{dEdx.dEdxTotIROC, dEdx.dEdxTotOROC1, dEdx.dEdxTotOROC2, dEdx.dEdxTotOROC3};
  const std::array<float, 4> dEdxMax{dEdx.dEdxMaxIROC, dEdx.dEdxMaxOROC1, dEdx.dEdxMaxOROC2, dEdx.dEdxMaxOROC3};

  for (const GEMstack stack : {IROCgem, OROC1gem, OROC2gem, OROC3gem}) {
    // These are the x value in cm of the center of the stacks (IROC, OROC1, ...) in the local frame.
    // FIXME: these values are only approximations.
    constexpr std::array<float, 4> xks{109.625f, 129.275f, 148.775f, 169.725f};
    constexpr float b = 0.5;

    // We need a copy of the track to perform propagations
    auto cpTrack = track;
    bool ok = cpTrack.propagateTo(xks[stack], b);

    // Ignore stack if we are not able to find its sector
    if (!ok) {
      continue;
    }

    float phi = cpTrack.getXYZGlo().Phi();
    float tgl = cpTrack.getTgl();
    float snp = cpTrack.getSnp();

    o2::math_utils::bringTo02PiGen(phi);
    const auto sector = static_cast<int>(phi / SECPHIWIDTH);

    mHist(dEdxMax[stack], tgl, snp, sector, side, stack, ChargeType::Max);
    mHist(dEdxTot[stack], tgl, snp, sector, side, stack, ChargeType::Tot);
  }
}

void CalibdEdx::fill(const gsl::span<const TrackTPC> tracks)
{
  for (const auto& track : tracks) {
    fill(track);
  }
}

void CalibdEdx::fill(const std::vector<TrackTPC>& tracks)
{
  for (const auto& track : tracks) {
    fill(track);
  }
}

void CalibdEdx::merge(const CalibdEdx* other)
{
  mHist += other->getHist();
}

void CalibdEdx::finalize(std::array<float, 2> minEntries)
{
  const float entries = minStackEntries();
  mCalib.clear();

  TLinearFitter fit(2);

  // Choose the fit dimension based on the available statistics
  if (entries >= minEntries[1]) {
    fit.SetFormula("1 ++ x ++ x*x ++ y ++ x*y ++ y*y");
    mCalib.setDims(2);
  } else if (entries >= minEntries[0]) {
    fit.SetFormula("1 ++ x ++ x*x");
    mCalib.setDims(1);
  } else {
    fit.SetFormula("1");
    mCalib.setDims(0);
  }

  const int stack_bins = mHist.axis(0).size() * mHist.axis(1).size() * mHist.axis(2).size();
  auto entry = bh::indexed(mHist).begin();
  for (size_t stack = 0; stack < 288; ++stack) {
    CalibdEdxCorrection::StackID stackID{
      static_cast<int>(entry->bin(HistAxis::Sector).lower()),
      static_cast<enum Side>(entry->bin(HistAxis::Side).lower()),
      static_cast<GEMstack>(entry->bin(HistAxis::Stack).lower())};
    const auto charge = static_cast<ChargeType>(entry->bin(HistAxis::Charge).lower());

    for (int bin = 0; bin < stack_bins; ++bin, ++entry) {
      const float counts = *entry;
      if (counts == 0) {
        continue;
      }
      const double dedx_val = entry->bin(HistAxis::dEdx).center();
      // take account of dedx resolution ~ 5%
      const double perr = dedx_val / sqrt(counts) * 0.05;
      std::array<double, 2> angles{entry->bin(HistAxis::Tgl).center(),
                                   entry->bin(HistAxis::Snp).center()};
      fit.AddPoint(angles.data(), dedx_val, perr);
    }
    fit.Eval();

    CalibdEdxCorrection::Params params{0};
    for (int param = 0; param < fit.GetNumberFreeParameters(); ++param) {
      params[param] = fit.GetParameter(param);
    }
    mCalib.setParams(stackID, charge, params);
    mCalib.setChi2(stackID, charge, fit.GetChisquare());
  }
}

float CalibdEdx::minStackEntries() const
{
  // sum over the dEdx bins to get the number of entries per stack
  auto projection = bh::algorithm::project(mHist, std::vector<int>{HistAxis::Sector, HistAxis::Side, HistAxis::Stack});
  auto dEdxCounts = bh::indexed(projection);
  // find the stack with the least number of entries
  auto min_it = std::min_element(dEdxCounts.begin(), dEdxCounts.end());
  // the count is doubled since we sum qMax and qTot entries
  return static_cast<float>(*min_it / 2);
}

bool CalibdEdx::hasEnoughData(float minEntries) const
{
  return minStackEntries() >= minEntries;
}

TH2F CalibdEdx::getRootHist(const std::vector<int>& projected_axis) const
{
  const float lower = mHist.axis(0).begin()->lower();
  const float upper = mHist.axis(0).end()->lower();

  auto projectedHist = getHist(projected_axis);

  const int nHists = projectedHist.size() / mNBins;

  TH2F rootHist("", "", nHists, 0, nHists, mNBins, lower, upper);

  int stack = 0;
  float last_center = -1;
  // fill TH2
  for (auto&& x : bh::indexed(projectedHist)) {
    const auto y = x.bin(0).center(); // current bin interval along dEdx axis
    const auto w = *x;                // "dereference" to get the bin value
    rootHist.Fill(stack, y, w);

    if (y < last_center) {
      stack++;
    }
    last_center = y;
  }

  return rootHist;
}

TH2F CalibdEdx::getRootHist() const
{
  std::vector<int> keep_all(HistAxis::Size);
  std::iota(keep_all.begin(), keep_all.end(), 0);
  return getRootHist(keep_all);
}

void CalibdEdx::print() const
{
  const int unique_entries = std::accumulate(mHist.begin(), mHist.end(), 0.0) / GEMSTACKSPERSECTOR / 2;
  LOGP(info, "Total number of track entries: {}", unique_entries);
}

void CalibdEdx::dumpToFile(std::string_view fileName) const
{
  TFile file(fileName.data(), "recreate");
  const auto rootHist = getRootHist();
  file.WriteObject(&rootHist, "CalibHists");
  file.Close();
}

void CalibdEdx::writeTTree(std::string_view fileName) const
{
  TFile f(fileName.data(), "recreate");

  TTree tree("hist", "Saving boost histogram to TTree");

  std::vector<float> row(mHist.rank());
  for (int i = 0; i < mHist.rank(); ++i) {
    // FIXME: infer axis type and remove the hardcoded float
    tree.Branch(mHist.axis(i).metadata().c_str(), &row[i]);
  }
  float count = 0;
  tree.Branch("counts", &count);

  for (const auto& x : indexed(mHist)) {
    for (int i = 0; i < mHist.rank(); ++i) {
      row[i] = x.bin(i).center();
    }
    count = *x;
    tree.Fill();
  }

  f.Write();
  f.Close();
}
