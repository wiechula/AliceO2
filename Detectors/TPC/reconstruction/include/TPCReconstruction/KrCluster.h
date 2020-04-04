// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file KrCluster.h
/// \brief Struct for Krypton and X-ray clusters
/// \author Philip Hauer <hauer@hiskp.uni-bonn.de>

#ifndef ALICEO2_TPC_KrCluster_H_
#define ALICEO2_TPC_KrCluster_H_

namespace o2
{
namespace tpc
{

struct KrCluster {
 public:
  float totCharge = 0;
  float maxCharge = 0;
  int size = 0;
  // means and sigmas are weighted averages/sigmas
  float meanPad = 0;
  float meanRow = 0;
  float meanTime = 0;
  float sigmaPad = 0;
  float sigmaRow = 0;
  float sigmaTime = 0;
  int sector = 0;
};

} // namespace tpc
} // namespace o2

#endif
