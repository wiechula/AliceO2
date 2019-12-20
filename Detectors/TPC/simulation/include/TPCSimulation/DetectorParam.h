// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DetectorParam.h
/// \brief Implementation of the parameter class for the TPC detector configuration
/// \author Jens Wiechula, Jens.Wiechula@ikf.uni-frankfurt.de

#ifndef ALICEO2_TPC_DetectorParam_H_
#define ALICEO2_TPC_DetectorParam_H_

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"

namespace o2
{
namespace tpc
{

struct DetectorParam : public o2::conf::ConfigurableParamHelper<DetectorParam> {
  enum class SimulationType : uint8_t {
    Physics = 0, /// default pysics simulation
    Krypton = 1  /// Krypton simulation
  };

  SimulationType SimType{SimulationType::Physics}; ///< Simulation type

  O2ParamDef(DetectorParam, "TPCDetector");
};
} // namespace tpc
} // namespace o2

#endif // ALICEO2_TPC_DetectorParam_H_
