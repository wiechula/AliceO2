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

/// \file CalibdEdxSpec.cxx
/// \brief Workflow for time based dE/dx calibration.
/// \author Thiago Badaró <thiago.saramela@usp.br>

#include "TPCWorkflow/CalibdEdxSpec.h"

// o2 includes
#include "DataFormatsTPC/TrackTPC.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsCommonDataFormats/NameConf.h"
#include "DetectorsCalibration/Utils.h"
#include "Framework/Task.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/ConfigParamRegistry.h"
#include "TPCCalibration/CalibdEdx.h"

using namespace o2::framework;

namespace o2::tpc
{

class CalibdEdxDevice : public Task
{
 public:
  void init(framework::InitContext& ic) final
  {
    const int minEntriesSector = ic.options().get<int>("min-entries-sector");
    const int minEntries1D = ic.options().get<int>("min-entries-1d");
    const int minEntries2D = ic.options().get<int>("min-entries-2d");

    const int dEdxBins = ic.options().get<int>("dedxbins");
    const int zBins = ic.options().get<int>("zbins");
    const int angularBins = ic.options().get<int>("angularbins");
    const float mindEdx = ic.options().get<float>("min-dedx");
    const float maxdEdx = ic.options().get<float>("max-dedx");

    const bool dumpData = ic.options().get<bool>("file-dump");
    float field = ic.options().get<float>("field");

    const auto inputGRP = o2::base::NameConf::getGRPFileName();
    const auto grp = o2::parameters::GRPObject::loadFrom(inputGRP);
    if (grp != nullptr) {
      field = 5.00668f * grp->getL3Current() / 30000.;
      LOGP(info, "Using GRP file to set the magnetic field to {} kG", field);
    }

    mCalib = std::make_unique<CalibdEdx>(mindEdx, maxdEdx, dEdxBins, zBins, angularBins);
    mCalib->setApplyCuts(false);
    mCalib->setFitCuts({minEntriesSector, minEntries1D, minEntries2D});
    mCalib->setField(field);

    mDumpToFile = dumpData;
  }

  void run(ProcessingContext& pc) final
  {
    const auto tfcounter = o2::header::get<DataProcessingHeader*>(pc.inputs().get("tracks").header)->startTime;
    const auto tracks = pc.inputs().get<gsl::span<TrackTPC>>("tracks");

    LOGP(info, "Processing TF {} with {} tracks", tfcounter, tracks.size());

    mCalib->fill(tracks);
    // sendOutput(pc.outputs());
  }

  void endOfStream(EndOfStreamContext& eos) final
  {
    LOGP(info, "Finalizing calibration");
    mCalib->finalize();
    mCalib->print();
    // sendOutput(eos.outputs());

    if (mDumpToFile) {
      mCalib->getCalib().saveFile("calibdEdx.root");
    }
  }

 private:
  void sendOutput(DataAllocator& output)
  {
    // extract CCDB infos and calibration objects, convert it to TMemFile and send them to the output
    // TODO in principle, this routine is generic, can be moved to Utils.h
  }

  bool mDumpToFile{};
  std::unique_ptr<CalibdEdx> mCalib;
};

DataProcessorSpec getCalibdEdxSpec()
{
  std::vector<OutputSpec> outputs;
  outputs.emplace_back(ConcreteDataTypeMatcher{"TPC", "dEdxCalibData"});

  return DataProcessorSpec{
    "tpc-calib-dEdx",
    Inputs{
      InputSpec{"tracks", "TPC", "MIPS"},
    },
    outputs,
    adaptFromTask<CalibdEdxDevice>(),
    Options{
      {"min-entries-sector", VariantType::Int, 100, {"minimum entries per stack to fit each GEM stack individually, bellow this all sectors are integrated"}},
      {"min-entries-1d", VariantType::Int, 500, {"minimum entries per stack to fit 1D correction"}},
      {"min-entries-2d", VariantType::Int, 2500, {"minimum entries per stack to fit 2D correction"}},

      {"dedxbins", VariantType::Int, 100, {"number of dEdx bins"}},
      {"zbins", VariantType::Int, 20, {"number of Z bins"}},
      {"angularbins", VariantType::Int, 18, {"number of bins for angular data, like Tgl and Snp"}},
      {"min-dedx", VariantType::Float, 5.0f, {"minimum value for the dEdx histograms"}},
      {"max-dedx", VariantType::Float, 100.0f, {"maximum value for the dEdx histograms"}},

      {"field", VariantType::Float, 5.f, {"magnetic field"}},
      {"file-dump", VariantType::Bool, false, {"directly dump calibration to file"}}}};
}

} // namespace o2::tpc
