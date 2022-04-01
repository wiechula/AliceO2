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

/// @file IACProcessorSpec.cxx
/// @brief TPC Integrated Analogue Current processing
/// @author Jens Wiechula

#include <vector>

#include "Framework/Task.h"
#include "Framework/InputRecordWalker.h"
#include "Framework/WorkflowSpec.h"
#include "Framework/ConfigParamRegistry.h"
#include "DetectorsRaw/RDHUtils.h"
#include "DPLUtils/RawParser.h"
#include "CommonUtils/TreeStreamRedirector.h"

#include "TPCBase/RDHUtils.h"
#include "TPCCalibration/IACDecoder.h"
#include "TPCWorkflow/IACProcessorSpec.h"

using namespace o2::framework;

namespace o2::tpc
{

class IACProcessorDevice : public o2::framework::Task
{
 public:
  using FEEIDType = rdh_utils::FEEIDType;

  void init(o2::framework::InitContext& ic) final
  {
    if (ic.options().get<bool>("write-debug-tree")) {
      mDecoder.enableDebugTree();
    }
    mDecoder.setDebugLevel(ic.options().get<uint32_t>("debug-level"));
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    const auto creationTime = DataRefUtils::getHeader<DataProcessingHeader*>(pc.inputs().getFirstValid(true))->creation;

    uint32_t tfCounter = 0;
    std::vector<InputSpec> filter = {{"check", ConcreteDataTypeMatcher{o2::header::gDataOriginTPC, "RAWDATA"}, Lifetime::Timeframe}}; // TODO: Change to IDC when changed in DD
    for (auto const& ref : InputRecordWalker(pc.inputs(), filter)) {
      const auto* dh = DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
      // ---| extract hardware information to do the processing |---
      const auto feeId = (FEEIDType)dh->subSpecification;
      const auto link = rdh_utils::getLink(feeId);
      const uint32_t cruID = rdh_utils::getCRU(feeId);
      const auto endPoint = rdh_utils::getEndPoint(feeId);
      tfCounter = dh->tfCounter;

      // only select IDCs
      // ToDo: cleanup once IDCs will be propagated not as RAWDATA, but IDC.
      if (link != rdh_utils::IDCLinkID) {
        continue;
      }
      LOGP(info, "IDC Processing firstTForbit {:9}, tfCounter {:5}, run {:6}, feeId {:6} ({:3}/{}/{:2})", dh->firstTForbit, dh->tfCounter, dh->runNumber, feeId, cruID, endPoint, link);

      // ---| data loop |---
      const gsl::span<const char> raw = pc.inputs().get<gsl::span<char>>(ref);
      o2::framework::RawParser parser(raw.data(), raw.size());
      for (auto it = parser.begin(), end = parser.end(); it != end; ++it) {
        const auto size = it.size();
        if (size == sizeof(0)) {
          auto* rdhPtr = it.get_if<o2::header::RAWDataHeaderV6>();
          if (!rdhPtr) {
            throw std::runtime_error("could not get RDH from packet");
          }
          if (rdhPtr->packetCounter == 0) {
            mDecoder.setReferenceTime(creationTime); // TODO set proper time
          }
          continue;
        }
        assert(size == sizeof(idc::Container));
        auto data = it.data();
      }
    }

    if (mDebugStream) {
      writeDebugOutput(tfCounter);
    }
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
    LOGP(info, "endOfStream");
    if (mDebugStream) {
      mDebugStream->Close();
    }
  }

 private:
  iac::Decoder mDecoder;                                         ///< Decoder for IAC data
  std::unique_ptr<o2::utils::TreeStreamRedirector> mDebugStream; ///< debug output streamer

  //____________________________________________________________________________
  void writeDebugOutput(uint32_t tfCounter)
  {
  }
};

o2::framework::DataProcessorSpec getIACProcessorSpec()
{
  using device = o2::tpc::IACProcessorDevice;

  std::vector<OutputSpec> outputs;

  return DataProcessorSpec{
    "tpc-iac-processor",
    select("A:TPC/RAWDATA"),
    outputs,
    AlgorithmSpec{adaptFromTask<device>()},
    Options{
      {"write-debug-tree", VariantType::Bool, false, {"write a debug output tree."}},
      {"debug-level", VariantType::UInt32, 0u, {"amount of debug to show"}},
    } // end Options
  };  // end DataProcessorSpec
}
} // namespace o2::tpc
