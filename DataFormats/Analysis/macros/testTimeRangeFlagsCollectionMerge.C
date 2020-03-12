#include <iostream>

#include "DetectorsCommonDataFormats/DetID.h"
#include "DataFormatsAnalysis/FlagReasons.h"
#include "DataFormatsAnalysis/TimeRangeFlagsCollection.h"
#include "CCDB/CcdbApi.h"

using namespace o2::analysis;

using time_type = uint64_t;
using o2::detectors::DetID;
using o2::ccdb::CcdbApi;

void testTime(DetID detID, time_type time, const TimeRangeFlagsCollection<time_type>& coll);

void testTimeRangeFlagsCollectionMerge()
{
  // ===| set up reasons |======================================================
  auto& reasons = FlagReasons::instance();
  reasons.addReason("Bad");
  reasons.addReason("Bad for PID");
  reasons.addReason("Limited acceptance");
  reasons.print();
  std::cout << "\n";

  // ===| set up mask ranges |==================================================
  TimeRangeFlagsCollection<time_type, uint16_t> coll1;
  coll1.addTimeRangeFlags(DetID::ITS, 0, 50, 1);
  coll1.addTimeRangeFlags(DetID::ITS, 51, 100, 4); 
  coll1.addTimeRangeFlags(DetID::TPC, 0, 9, 6);
  TimeRangeFlagsCollection<time_type, uint16_t> coll2;
  coll2.addTimeRangeFlags(DetID::TPC, 0, 1, 6);
  coll2.addTimeRangeFlags(DetID::TPC, 10, 100, 4);

  printf("Collection one:\n");
  coll1.print();
  std::cout << "\n";
  printf("Collection two:\n");
  coll2.print();

  // ===| test if time has flags |==============================================
  std::cout << "==================| check times for flags |====================\n";
  std::cout << "                     " << DetID(DetID::ITS).getName() << "\n";
  testTime(DetID::ITS, 2, coll1);
  testTime(DetID::ITS, 50, coll1);
  std::cout << "\n";

  std::cout << "                     " << DetID(DetID::TPC).getName() << "\n";
  testTime(DetID::TPC, 2, coll2);
  testTime(DetID::TPC, 50, coll2);

  coll2.merge(coll1);
  std::cout << "\n";
  printf("NEW Collection two:\n");
  coll2.print();
  

}

void testTime(DetID detID, time_type time, const TimeRangeFlagsCollection<time_type>& coll)
{
  const TimeRangeFlags<time_type>* flags{nullptr};
  if ((flags = coll.findTimeRangeFlags(detID, time))) {
    std::cout << "Time " << time << " has the flags: " << flags->collectMaskReasonNames() << "\n";
  } else {
    std::cout << "Time " << time << " has NO flags\n";
  }
}

