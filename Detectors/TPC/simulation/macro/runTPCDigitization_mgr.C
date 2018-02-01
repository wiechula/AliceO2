/*
 * filterHits.C
 *
 *  Created on: Jan 30, 2018
 *      Author: swenzel
 */

// R__ADD_INCLUDE_PATH($O2_ROOT/../../ms_gsl/latest/include)
R__ADD_INCLUDE_PATH(/data/Work/software/alicesw/sw/ubuntu1404_x86-64/ms_gsl/latest/include)
R__ADD_INCLUDE_PATH($VC_ROOT/include)
R__ADD_LIBRARY_PATH($O2_ROOT/lib)
R__LOAD_LIBRARY(libSimulationDataFormat)

#include <Steer/HitProcessingManager.h>
#include <TPCSimulation/Digitizer.h>
#include <TPCSimulation/DigitizerTask.h>
#include <cassert>
#include <functional>
#include "ITSMFTSimulation/Hit.h"
#include "TPCBase/Sector.h"
#include "TPCSimulation/Point.h"
#include "TPCSimulation/SAMPAProcessing.h"

void getHits(const o2::steer::RunContext& context, std::vector<std::vector<o2::TPC::HitGroup>*>& hitvectors,
             std::vector<o2::TPC::TPCHitGroupID>& hitids,

             const char* branchname, float tmin /*NS*/, float tmax /*NS*/,
             std::function<float(float, float, float)>&& f)
{
  // f is some function taking event time + z of hit and returns final "digit" time

  // ingredients
  // *) tree

  auto br = context.getBranch(branchname);
  if (!br) {
    std::cerr << "No branch found\n";
    return;
  }

  auto& eventrecords = context.getEventRecords();

  auto nentries = br->GetEntries();
  hitvectors.resize(nentries, nullptr);

  // do the filtering
  for (int entry = 0; entry < nentries; ++entry) {
    if (tmin > f(eventrecords[entry].timeNS, 0, 0)) {
      continue;
    }
    if (tmax < f(eventrecords[entry].timeNS, 0, 250)) {
      break;
    }

    br->SetAddress(&hitvectors[entry]);
    br->GetEntry(entry);

    int groupid = -1;
    auto groups = hitvectors[entry];
    for (auto& singlegroup : *groups) {
      groupid++;
      auto zmax = singlegroup.mZAbsMax;
      auto zmin = singlegroup.mZAbsMin;
      // in case of secondaries, the time ordering may be reversed
      if (zmax < zmin)
        std::swap(zmax, zmin);
      assert(zmin <= zmax);
      // auto tof = singlegroup.
      float tmaxtrack = f(eventrecords[entry].timeNS, 0., zmin);
      float tmintrack = f(eventrecords[entry].timeNS, 0., zmax);
      std::cout << tmintrack << " & " << tmaxtrack << "\n";
      assert(tmaxtrack >= tmintrack);
      if (tmin > tmaxtrack || tmax < tmintrack) {
        std::cout << "DISCARDING " << groupid << " OF ENTRY " << entry << "\n";
        continue;
      }
      // need to record index of the group
      hitids.emplace_back(entry, groupid);
    }
  }
}

// TPC hit selection lambda
auto fTPC = [](float tNS, float tof, float z) {
  // returns time in NS
  return tNS + o2::TPC::SAMPAProcessing::getDriftTime(z) * 1000 + tof;
};

void runTPCDigitization(const o2::steer::RunContext& context)
{
  std::vector<std::vector<o2::TPC::HitGroup>*> hitvectorsleft;  // "TPCHitVector"
  std::vector<std::vector<o2::TPC::HitGroup>*> hitvectorsright; // "TPCHitVector"
  std::vector<o2::TPC::TPCHitGroupID> hitidsleft;               // "TPCHitIDs"
  std::vector<o2::TPC::TPCHitGroupID> hitidsright;

  o2::TPC::DigitizerTask task;
  task.Init2();

  // ===| open file and register branches |=====================================
  std::unique_ptr<TFile> file = std::unique_ptr<TFile>(TFile::Open("tpc_digi.root", "recreate"));
  TTree outtree("tpc_digits", "TPC digits");

  using digitType = std::vector<o2::TPC::Digit>;
  using mcType = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
  std::array<digitType*, 36> digitArrays;
  std::array<mcType*, 36> mcTruthArrays;
  std::array<TBranch*, 36> digitBranches;
  std::array<TBranch*, 36> mcTruthBranches;

  for (int sector = 0; sector < 36; ++sector) {
    digitArrays[sector] = new digitType;
    digitBranches[sector] = outtree.Branch(Form("TPCDigit_%i", sector), &digitArrays[sector]);

    // Register MC Truth container
    mcTruthArrays[sector] = new mcType;
    mcTruthBranches[sector] = outtree.Branch(Form("TPCDigitMCTruth_%i", sector), &mcTruthArrays[sector]);
  }

  const auto TPCDRIFT = 100000;
  for (int driftinterval = 0;driftinterval<1; ++driftinterval) {
    auto tmin = driftinterval * TPCDRIFT;
    auto tmax = (driftinterval + 1) * TPCDRIFT;

    std::cout << "Driftinterval " << driftinterval << " times: " << tmin << " - " << tmax << "\n";

    hitvectorsleft.clear();
    hitidsleft.clear();
    hitvectorsright.clear();
    hitidsright.clear();

    bool hasData = false;
    // loop over sectors
    for (int s = 0; s < 36; ++s) {
      task.setSector(s);
      task.setOutputData(digitArrays[s], mcTruthArrays[s]);
      std::stringstream sectornamestreamleft;
      sectornamestreamleft << "TPCHitsShiftedSector" << s;
      getHits(context, hitvectorsleft, hitidsleft, sectornamestreamleft.str().c_str(), tmin, tmax, fTPC);

      std::stringstream sectornamestreamright;
      sectornamestreamright << "TPCHitsShiftedSector" << int(o2::TPC::Sector::getLeft(s));
      getHits(context, hitvectorsright, hitidsright, sectornamestreamright.str().c_str(), tmin, tmax, fTPC);

      task.setData(&hitvectorsleft, &hitvectorsright, &hitidsleft, &hitidsright, &context);
      task.setTimeRange(tmin, tmax);
      task.Exec2("");

      digitBranches[s]->Fill();
      mcTruthBranches[s]->Fill();

      hasData |= hitidsleft.size() || hitidsright.size();
    }

    // condition to end:
    if (!hasData) {
      break;
    }
    task.FinishTask();
  }

  outtree.Fill();
  file->Write();
}

// void runITSDigitization(const o2::steer::RunContext& context)
//{
//  std::vector<o2::ITSMFT::Hit>* hitvector = nullptr;
//  o2::ITSMFT::DigitizerTask task(true);
//  task.Init();
//  auto br = context.getBranch("ITSHit");
//  assert(br);
//  br->SetAddress(&hitvector);
//  for (int entry = 0; entry < context.getNEntries(); ++entry) {
//    br->GetEntry(entry);
//    task.setData(hitvector, &context);
//    task.Exec("");
//  }
//  task.FinishTask();
//}

void runTPCDigitization_mgr()
{
  auto& hitrunmgr = o2::steer::HitProcessingManager::instance();
  hitrunmgr.addInputFile("o2sim.root");
  hitrunmgr.registerRunFunction(runTPCDigitization);

  // hitrunmgr.registerRunFunction(runITSDigitization);
  using Data_t = std::vector<o2::ITSMFT::Hit>;

  //  o2::ITS::DigitizerTask task;
  //  TGeoManager::Import("O2geometry.root");
  //  task.Init();
  //  hitrunmgr.registerDefaultRunFunction(o2::steer::defaultRunFunction(
  //	task, "ITSHit"));

  hitrunmgr.run();
}
