// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file KrBoxClusterFinder3D.cpp
/// \brief Class source code for Krypton and X-ray events
/// \author Philip Hauer <hauer@hiskp.uni-bonn.de>

#include "TPCReconstruction/KrBoxClusterFinder.h"
#include "Framework/Logger.h"

using namespace o2::tpc;

// Constructor:
// This function creates a three dimensional vector (Pad,Row,Time) which is
// filled with the charges
KrBoxClusterFinder::KrBoxClusterFinder(std::vector<o2::tpc::Digit>& eventSector)
{
  std::vector<o2::tpc::Digit>::iterator cl = eventSector.begin();
  std::vector<int> Rows, Pads, Times;

  if (eventSector.size() == 0) {
    // prevents a segementation fault if the envent contains no data
    // segementation fault would occure later when trying to dereference the
    // max/min pointer (which are nullptr if data is empty) empty events should
    // be catched in the main function, hence, this "if" is no longer necessary
    LOGP(warning, "Sector size amount of data points) in current run is 0!");
    LOGP(warning, "mMapOfAllDigits with 0's is generated in order to prevent a segementation fault.");

    mMapOfAllDigits = std::vector<std::vector<std::vector<int>>>(10, std::vector<std::vector<int>>(10, std::vector<int>(10, 0)));
  } else {
    for (const auto& cli : eventSector) {
      Times.emplace_back(cli.getTimeStamp());
      Pads.emplace_back(cli.getPad());
      Rows.emplace_back(cli.getRow());
    }
    // here a "mMapOfAllDigits" large enough to contain all data points of the envents is
    // initialized and then filled
    mMapOfAllDigits = std::vector<std::vector<std::vector<int>>>(mMaxTimes, std::vector<std::vector<int>>(mMaxRows, std::vector<int>(mMaxPads)));
    for (int i = 0; (cl + i) < eventSector.end(); i++) {
      mMapOfAllDigits[Times[i]][Rows[i]][Pads[i]] = (cl + i)->getChargeFloat();
    }
  }
}

//#################################################

// Function to update the temporal cluster
void KrBoxClusterFinder::updateTempClusterFinal()
{
  double one_over_q_tot = 1. / mTempCluster.totCharge;
  mTempCluster.meanPad *= one_over_q_tot;
  mTempCluster.sigmaPad *= one_over_q_tot;
  mTempCluster.meanRow *= one_over_q_tot;
  mTempCluster.sigmaRow *= one_over_q_tot;
  mTempCluster.meanTime *= one_over_q_tot;
  mTempCluster.sigmaTime *= one_over_q_tot;
  mTempCluster.sigmaPad = std::sqrt(std::abs(mTempCluster.sigmaPad - mTempCluster.meanPad * mTempCluster.meanPad));
  mTempCluster.sigmaRow = std::sqrt(std::abs(mTempCluster.sigmaRow - mTempCluster.meanRow * mTempCluster.meanRow));
  mTempCluster.sigmaTime = std::sqrt(std::abs(mTempCluster.sigmaTime - mTempCluster.meanTime * mTempCluster.meanTime));
}

// Function to update the temporal cluster.
void KrBoxClusterFinder::updateTempCluster(int temp_Charge, int temp_Pad, int temp_Row, int temp_Time)
{
  if (temp_Charge >= mQThreshold) {
    mTempCluster.size += 1;
    mTempCluster.totCharge += temp_Charge;

    mTempCluster.meanPad += temp_Pad * temp_Charge;
    mTempCluster.sigmaPad += temp_Pad * temp_Pad * temp_Charge;

    mTempCluster.meanRow += temp_Row * temp_Charge;
    mTempCluster.sigmaRow += temp_Row * temp_Row * temp_Charge;

    mTempCluster.meanTime += temp_Time * temp_Charge;
    mTempCluster.sigmaTime += temp_Time * temp_Time * temp_Charge;
    if (temp_Charge > mTempCluster.maxCharge) {
      mTempCluster.maxCharge = temp_Charge;
    }
  } else {
    LOGP(warning, "Update cluster was called but current charge is below mQThreshold");
  }
}

// This function finds and evaluates all clusters in a 3D mMapOfAllDigits generated by the
// mMapOfAllDigitsCreator function, this function also updates the cluster tree
std::vector<std::tuple<int, int, int>> KrBoxClusterFinder::findLocalMaxima()
{
  std::vector<std::tuple<int, int, int>> localMaximaCoords;
  // loop over whole mMapOfAllDigits the find clusers
  for (int iRow = 0; iRow < mMaxRows; iRow++) {
    for (int iPad = 0; iPad < mMaxPads; iPad++) {
      for (int iTime = 0; iTime < mMaxTimes; iTime++) {

        mTempCluster.Reset();
        const int qMax = mMapOfAllDigits.at(iTime).at(iRow).at(iPad);

        // cluster Maximum must at least be larger than Threshold
        if (qMax <= mQThresholdMax) {
          continue;
        }

        // Acceptance condition: Require at least mMinNumberOfNeighbours neigbours
        // with signal in any direction!
        int noNeighbours = 0;
        if ((iPad + 1 < mMaxPads) && (mMapOfAllDigits[iTime][iRow][iPad + 1] > mQThreshold)) {
          noNeighbours++;
        }
        if ((iPad - 1 >= 0) && (mMapOfAllDigits[iTime][iRow][iPad - 1] > mQThreshold)) {
          noNeighbours++;
        }
        if ((iTime + 1 < mMaxTimes) && (mMapOfAllDigits[iTime + 1][iRow][iPad] > mQThreshold)) {
          noNeighbours++;
        }
        if ((iTime - 1 >= 0) && (mMapOfAllDigits[iTime - 1][iRow][iPad] > mQThreshold)) {
          noNeighbours++;
        }
        if ((iRow + 1 < mMaxRows) && (mMapOfAllDigits[iTime][iRow + 1][iPad] > mQThreshold)) {
          noNeighbours++;
        }
        if ((iRow - 1 >= 0) && (mMapOfAllDigits[iTime][iRow - 1][iPad] > mQThreshold)) {
          noNeighbours++;
        }
        if (noNeighbours < mMinNumberOfNeighbours) {
          continue;
        }

        // Check that this is a local maximum
        // Note that the checking is done so that if 2 charges have the same
        // qMax then only 1 cluster is generated
        // (that is why there is BOTH > and >=)
        // -> only the maximum with the smalest indices will be accepted
        bool thisIsMax = true;

        for (int i = -mMaxClusterSizePad; (i <= mMaxClusterSizePad) && thisIsMax;
             i++) {
          for (int k = -mMaxClusterSizeRow; (k <= mMaxClusterSizeRow) && thisIsMax;
               k++) {
            for (int j = -mMaxClusterSizeTime;
                 (j <= mMaxClusterSizeTime) && thisIsMax; j++) {
              if ((iPad + i < mMaxPads) && (iTime + j < mMaxTimes) &&
                  (iRow + k < mMaxRows) && (iPad + i >= 0) && (iTime + j) >= 0 &&
                  (iRow + k) >= 0 &&
                  mMapOfAllDigits[iTime + j][iRow + k][iPad + i] > qMax) {
                thisIsMax = false;
              }
            }
          }
        }

        if (!thisIsMax) {
          continue;
        } else {
          localMaximaCoords.push_back(std::make_tuple(iPad, iRow, iTime));
        }
      }
    }
  }
  return localMaximaCoords;
}

// Calculate the total charge as the sum over the region:
//
//    o o o o o
//    o i i i o
//    o i C i o
//    o i i i o
//    o o o o o
//
// with qmax at the center C.
//
// The inner charge (i) we always add, but we only add the outer
// charge (o) if the neighboring inner bin (i) has a signal.
//

// for loop over whole cluster, to determine if a charge should be added
// conditions are extrapolation of the 5x5 cluster case to arbitrary
// cluster sizes in 3 dimensions
KrCluster KrBoxClusterFinder::buildCluster(int clusterCenterPad, int clusterCenterRow, int clusterCenterTime)
{
  mTempCluster = KrCluster();

  for (int iTime = -mMaxClusterSizeTime; iTime <= mMaxClusterSizeTime; iTime++) {
    for (int iRow = -mMaxClusterSizeRow; iRow <= mMaxClusterSizeRow; iRow++) {
      for (int iPad = -mMaxClusterSizePad; iPad <= mMaxClusterSizePad; iPad++) {
        // First: Check if we might check outside of map:
        if (clusterCenterTime + iTime < 0 || clusterCenterTime + iTime >= mMaxTimes || clusterCenterPad + iPad < 0 || clusterCenterPad + iPad >= mMaxPads || clusterCenterRow + iRow < 0 || clusterCenterRow + iRow >= mMaxRows || mMapOfAllDigits.at(clusterCenterTime + iTime).at(clusterCenterRow + iRow).at(clusterCenterPad + iPad) <= mQThreshold) {
          continue;
        }
        // If not, there are several cases which were explained (for 2D) in the header of the code.
        // The first one is for the diagonal. So, the digit we are investigating here is on the diagonal:
        else if (std::abs(iTime) == std::abs(iPad) && std::abs(iTime) == std::abs(iRow)) {
          // Now we check, if the next inner digit has a signal above threshold:
          if (mMapOfAllDigits[clusterCenterTime + iTime - signnum(iTime)][clusterCenterRow + iRow - signnum(iRow)][clusterCenterPad + iPad - signnum(iPad)] > mQThreshold) {
            // If yes, the cluster gets updated with the digit on the diagonal.
            updateTempCluster(mMapOfAllDigits[clusterCenterTime + iTime][clusterCenterRow + iRow][clusterCenterPad + iPad], clusterCenterPad + iPad, clusterCenterRow + iRow, clusterCenterTime + iTime);
          }
        }
        // Basically, we go through every possible case in the next few if-else conditions:
        else if (std::abs(iTime) == std::abs(iPad)) {
          if (mMapOfAllDigits[clusterCenterTime + iTime - signnum(iTime)][clusterCenterRow + iRow][clusterCenterPad + iPad - signnum(iPad)] > mQThreshold) {
            updateTempCluster(mMapOfAllDigits[clusterCenterTime + iTime][clusterCenterRow + iRow][clusterCenterPad + iPad], clusterCenterPad + iPad, clusterCenterRow + iRow, clusterCenterTime + iTime);
          }
        } else if (std::abs(iTime) == std::abs(iRow)) {
          if (mMapOfAllDigits[clusterCenterTime + iTime - signnum(iTime)][clusterCenterRow + iRow - signnum(iRow)][clusterCenterPad + iPad] > mQThreshold) {
            updateTempCluster(mMapOfAllDigits[clusterCenterTime + iTime][clusterCenterRow + iRow][clusterCenterPad + iPad], clusterCenterPad + iPad, clusterCenterRow + iRow, clusterCenterTime + iTime);
          }
        } else if (std::abs(iPad) == std::abs(iRow)) {
          if (mMapOfAllDigits[clusterCenterTime + iTime][clusterCenterRow + iRow - signnum(iRow)][clusterCenterPad + iPad - signnum(iPad)] > mQThreshold) {
            updateTempCluster(mMapOfAllDigits[clusterCenterTime + iTime][clusterCenterRow + iRow][clusterCenterPad + iPad], clusterCenterPad + iPad, clusterCenterRow + iRow, clusterCenterTime + iTime);
          }
        } else if (std::abs(iTime) > std::abs(iPad) && std::abs(iTime) > std::abs(iRow)) {
          if (mMapOfAllDigits[clusterCenterTime + iTime - signnum(iTime)][clusterCenterRow + iRow][clusterCenterPad + iPad] > mQThreshold) {
            updateTempCluster(mMapOfAllDigits[clusterCenterTime + iTime][clusterCenterRow + iRow][clusterCenterPad + iPad], clusterCenterPad + iPad, clusterCenterRow + iRow, clusterCenterTime + iTime);
          }
        } else if (std::abs(iTime) < std::abs(iPad) && std::abs(iPad) > std::abs(iRow)) {
          if (mMapOfAllDigits[clusterCenterTime + iTime][clusterCenterRow + iRow][clusterCenterPad + iPad - signnum(iPad)] > mQThreshold) {
            updateTempCluster(mMapOfAllDigits[clusterCenterTime + iTime][clusterCenterRow + iRow][clusterCenterPad + iPad], clusterCenterPad + iPad, clusterCenterRow + iRow, clusterCenterTime + iTime);
          }
        } else if (std::abs(iTime) < std::abs(iRow) && std::abs(iPad) < std::abs(iRow)) {
          if (mMapOfAllDigits[clusterCenterTime + iTime][clusterCenterRow + iRow - signnum(iRow)][clusterCenterPad + iPad] > mQThreshold) {
            updateTempCluster(mMapOfAllDigits[clusterCenterTime + iTime][clusterCenterRow + iRow][clusterCenterPad + iPad], clusterCenterPad + iPad, clusterCenterRow + iRow, clusterCenterTime + iTime);
          }
        }
      }
    }
  }
  // At the end, out mTempCluster should contain all digits that were assigned to the cluster.
  // So before returning it, we update it one last time to calculate the correct means and sigmas.
  updateTempClusterFinal();
  return mTempCluster;
}
