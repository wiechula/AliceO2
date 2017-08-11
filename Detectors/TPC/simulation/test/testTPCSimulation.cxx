// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file testTPCSimulation.cxx
/// \brief This task tests several small components of the TPC simulation
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#define BOOST_TEST_MODULE Test TPC Simulation
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "TPCSimulation/Point.h"
#include "TPCSimulation/DigitMC.h"
#include "TPCSimulation/DigitMCMetaData.h"

template <typename T>
using Point3D = ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<T>, ROOT::Math::DefaultCoordinateSystemTag>;

namespace o2 {
namespace TPC {

  /// \brief Trivial test of the initialization of a Point and its getters
  /// Precision: 1E-12 %
  BOOST_AUTO_TEST_CASE(Point_test)
  {
    Point testpoint(2.f, 3.f, 4.f, 5.f, 6, 7, 8);
    BOOST_CHECK_CLOSE(testpoint.GetX(),2.,1E-12);
    BOOST_CHECK_CLOSE(testpoint.GetY(),3.,1E-12);
    BOOST_CHECK_CLOSE(testpoint.GetZ(),4.,1E-12);
    BOOST_CHECK_CLOSE(testpoint.GetTime(),5.,1E-12);
    BOOST_CHECK_CLOSE(testpoint.GetEnergyLoss(),6,1E-12);
    BOOST_CHECK_CLOSE(testpoint.GetTrackID(),7.,1E-12);
    BOOST_CHECK_CLOSE(testpoint.GetDetectorID(),8.,1E-12);
  }

  /// \brief Trivial test of the initialization of a DigitMC and its getters
  /// Precision: 1E-12 %
  BOOST_AUTO_TEST_CASE(DigitMC_test)
  {
    const std::vector<long> testMC = {{1000001, 2000002, 3000033}};
    DigitMC testdigit(testMC, 1, 2.f, 3, 4, 5);
    BOOST_CHECK(testdigit.getCRU() == 1);
    BOOST_CHECK_CLOSE(testdigit.getCharge(), 2.f, 1E-12);
    BOOST_CHECK(testdigit.getRow() == 3);
    BOOST_CHECK(testdigit.getPad() == 4);
    BOOST_CHECK(testdigit.getTimeStamp() == 5);
    BOOST_CHECK(testdigit.getNumberOfMClabels() == testMC.size());
    BOOST_CHECK(testdigit.getMCEvent(0) == 1);
    BOOST_CHECK(testdigit.getMCEvent(1) == 2);
    BOOST_CHECK(testdigit.getMCEvent(2) == 3);
    BOOST_CHECK(testdigit.getMCTrack(0) == 1);
    BOOST_CHECK(testdigit.getMCTrack(1) == 2);
    BOOST_CHECK(testdigit.getMCTrack(2) == 33);
  }

  /// \brief Trivial test of the initialization of a DigitMCMetaData and its getters
  /// Precision: 1E-12 %
  BOOST_AUTO_TEST_CASE(DigitMCMetaData_test)
  {
    DigitMCMetaData testdigit(1.f, 2.f, 3.f, 4.f);
    BOOST_CHECK_CLOSE(testdigit.getRawADC(), 1.f, 1E-12);
    BOOST_CHECK_CLOSE(testdigit.getCommonMode(), 2.f,1E-12);
    BOOST_CHECK_CLOSE(testdigit.getPedestal(), 3.f,1E-12);
    BOOST_CHECK_CLOSE(testdigit.getNoise(), 4.f,1E-12);
  }
}
}
