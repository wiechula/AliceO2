// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DigitContainer.cxx
/// \brief Implementation of the Digit Container
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#include "TPCSimulation/DigitContainer.h"
#include <iostream>
#include <cassert>

using namespace o2::TPC;

void DigitContainer::addDigit(size_t eventID, size_t hitID, const CRU &cru, TimeBin timeBin, GlobalPadNumber globalPad, float charge)
{
  const int sector = cru.sector();
  mSector.setDigit(eventID, hitID, cru, timeBin, globalPad, charge);
}


void DigitContainer::fillOutputContainer(std::vector<Digit> *output, dataformats::MCTruthContainer<MCCompLabel> &mcTruth,
                                         std::vector<DigitMCMetaData> *debug, TimeBin eventTime, bool isContinuous, bool isFinal)
{
    mSector.fillOutputContainer(output, mcTruth, debug, mSectorID, eventTime, isContinuous);
    mSector.reset();
}



