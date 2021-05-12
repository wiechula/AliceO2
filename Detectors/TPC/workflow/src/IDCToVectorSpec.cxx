// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>
#include "fmt/format.h"

#include "TFile.h"
#include "Framework/Task.h"
#include "Framework/ControlService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/Logger.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/WorkflowSpec.h"
#include "Framework/InputRecordWalker.h"
#include "DPLUtils/RawParser.h"
#include "Headers/DataHeader.h"
#include "CommonUtils/TreeStreamRedirector.h"

#include "DataFormatsTPC/Defs.h"
#include "DataFormatsTPC/IDC.h"
#include "TPCBase/RDHUtils.h"
#include "TPCBase/Mapper.h"

using namespace o2::framework;
using o2::header::gDataOriginTPC;

namespace o2::tpc
{

class IDCToVectorDevice : public o2::framework::Task
{
 public:
  // TODO: Remove once in mapper and add Mapper:: in the places below
  static constexpr unsigned int PADSPERREGION[10]{1200, 1200, 1440, 1440, 1440, 1440, 1600, 1600, 1600, 1600};
  static constexpr unsigned int GLOBALPADOFFSET[10]{0, 1200, 2400, 3840, 5280, 6720, 8160, 9760, 11360, 12960};

  using FEEIDType = rdh_utils::FEEIDType;
  IDCToVectorDevice(const std::vector<uint32_t>& crus) : mCRUs(crus) {}

  void init(o2::framework::InitContext& ic) final
  {
    // set up ADC value filling
    if (ic.options().get<bool>("write-debug")) {
      mDebugStream = std::make_unique<o2::utils::TreeStreamRedirector>("idc_vector_debug.root", "recreate");
    }

    initIDC();
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    std::vector<InputSpec> filter = {{"check", ConcreteDataTypeMatcher{o2::header::gDataOriginTPC, "RAWDATA"}, Lifetime::Timeframe}}; // TODO: Change to IDC when changed in DD
    const auto& mapper = Mapper::instance();

    uint32_t heartbeatOrbit = 0;
    uint32_t heartbeatBC = 0;
    bool first = true;

    for (auto const& ref : InputRecordWalker(pc.inputs(), filter)) {
      const auto* dh = DataRefUtils::getHeader<o2::header::DataHeader*>(ref);

      // ---| extract hardware information to do the processing |---
      const auto feeId = (FEEIDType)dh->subSpecification;
      const auto link = rdh_utils::getLink(feeId);
      const uint32_t cruID = rdh_utils::getCRU(feeId);
      const auto endPoint = rdh_utils::getEndPoint(feeId);

      // only select IDCs
      // ToDo: cleanup once IDCs will be propagated not as RAWDATA, but IDC.
      if (link != rdh_utils::IDCLinkID) {
        continue;
      }
      LOGP(info, "IDC Processing firstTForbit {:9}, tfCounter {:5}, run {:6}, feeId {:6} ({:3}/{}/{:2})", dh->firstTForbit, dh->tfCounter, dh->runNumber, feeId, cruID, endPoint, link);

      if (std::find(mCRUs.begin(), mCRUs.end(), cruID) == mCRUs.end()) {
        LOGP(error, "IDC CRU {:3} not configured in CRUs, skipping", cruID);
        continue;
      }

      const CRU cru(cruID);
      const auto& partInfo = mapper.getPartitionInfo(cru.partition());
      const int fecLinkOffsetCRU = (partInfo.getNumberOfFECs() + 1) / 2;
      const int fecSectorOffset = partInfo.getSectorFECOffset();
      const GlobalPadNumber regionPadOffset = GLOBALPADOFFSET[cru.region()]; // ToDo: Add Mapper:: once in mapper
      const GlobalPadNumber numberPads = PADSPERREGION[cru.region()];        // ToDo: Add Mapper:: once in mapper
      int sampaOnFEC{}, channelOnSAMPA{};
      auto& idcVec = mIDCvectors[cruID];
      auto& infoVec = mIDCInfos[cruID];

      // ---| data loop |---
      const gsl::span<const char> raw = pc.inputs().get<gsl::span<char>>(ref);
      o2::framework::RawParser parser(raw.data(), raw.size());
      for (auto it = parser.begin(), end = parser.end(); it != end; ++it) {
        const auto size = it.size();
        assert(size == sizeof(idc::Container));
        auto data = it.data();
        auto& idcs = *((idc::Container*)(data));
        const uint32_t orbit = idcs.header.heartbeatOrbit;
        const uint32_t bc = idcs.header.heartbeatBC;

        auto infoIt = std::find(infoVec.begin(), infoVec.end(), orbit);
        if (!infoVec.size()) {
          infoVec.emplace_back(orbit, bc);
          infoIt = infoVec.end() - 1;
          //} else if (!infoVec.back().matches(orbit, bc)) {
        } else if (infoIt == infoVec.end()) {
          auto& lastInfo = infoVec.back();
          if ((orbit - lastInfo.heartbeatOrbit) != mNTFsIDC) {
            LOGP(error, "received packet with invalid jump in idc orbit ({} - {} == {} != {})", orbit, lastInfo.heartbeatOrbit, orbit - lastInfo.heartbeatOrbit, mNTFsIDC);
          }
          infoVec.emplace_back(orbit, bc);
          infoIt = infoVec.end() - 1;
        }

        // check if end poit was already processed
        //auto& lastInfo = infoVec.back();
        auto& lastInfo = *infoIt;
        if (lastInfo.wasEPseen(endPoint)) {
          LOGP(info, "Already received another data packet for CRU {}, ep {}, orbit {}, bc {}", cruID, endPoint, orbit, bc);
          continue;
        }

        lastInfo.setEPseen(endPoint);
        const size_t idcOffset = infoVec.size() - 1;
        LOGP(info, "processing IDCs for CRU {}, ep {}, orbit {}, bc {}, idcOffset {}, idcVec size {}, epSeen {:02b}", cruID, endPoint, orbit, bc, idcOffset, idcVec.size(), lastInfo.epSeen);

        for (uint32_t iLink = 0; iLink < idc::Links; ++iLink) {
          if (!idcs.hasLink(iLink)) {
            continue;
          }
          const int fecInSector = iLink + endPoint * fecLinkOffsetCRU + fecSectorOffset;

          for (uint32_t iChannel = 0; iChannel < idc::Channels; ++iChannel) {
            const auto val = idcs.getChannelValueFloat(iLink, iChannel);
            Mapper::getSampaAndChannelOnFEC(cruID, iChannel, sampaOnFEC, channelOnSAMPA);
            const GlobalPadNumber padInSector = mapper.globalPadNumber(fecInSector, sampaOnFEC, channelOnSAMPA);
            const GlobalPadNumber padInRegion = padInSector - regionPadOffset;
            const GlobalPadNumber vectorPosition = padInRegion + idcOffset * numberPads;
            LOGP(info, "filling channel {} idcVec[{} ({})] = {}", iChannel, vectorPosition, padInRegion, val);
            idcVec[vectorPosition] = val;
          }
        }
      }
    }

    if (mDebugStream) {
      writeDebugOutput();
    }
    snapshotIDCs(pc.outputs());
  }

  void endOfStream(o2::framework::EndOfStreamContext& ec) final
  {
    LOGP(info, "endOfStream");
    ec.services().get<ControlService>().readyToQuit(QuitRequest::Me);
    if (mDebugStream) {
      mDebugStream->Close();
    }
  }

 private:
  struct IDCInfo {
    IDCInfo() = default;
    IDCInfo(const IDCInfo&) = default;
    IDCInfo(uint32_t orbit, uint16_t bc) : heartbeatOrbit(orbit), heartbeatBC(bc) {}

    uint32_t heartbeatOrbit{0};
    uint16_t heartbeatBC{0};
    uint16_t epSeen{0};

    bool operator==(const uint32_t orbit) const { return (heartbeatOrbit == orbit); }
    bool operator==(const IDCInfo& inf) const { return (inf.heartbeatOrbit == heartbeatOrbit) && (inf.heartbeatBC == heartbeatBC) && (inf.epSeen == epSeen); }
    void setEPseen(uint32_t ep) { epSeen |= uint16_t(1 << ep); }
    bool wasEPseen(uint32_t ep) const { return epSeen & uint16_t(1 << ep); }
    bool matches(uint32_t orbit, int16_t bc) const { return ((heartbeatOrbit == orbit) && (heartbeatBC == bc)); }
    bool hasBothEPs() const { return epSeen == 3; }
  };

  uint32_t mNTFsIDC{12};                                         ///< Number of time frames over which IDCs are integrated
  static const uint32_t mMaxIDCPerTF{13};                        ///< maximum number of IDCs expected per TF
  uint32_t mIDCsReceived{0};                                     ///< number of IDC packages received for one TF
  std::vector<uint32_t> mCRUs;                                   ///< CRUs expected for this device
  std::unordered_map<uint32_t, std::vector<float>> mIDCvectors;  ///< decoded IDCs per cru for each pad in the region
  std::unordered_map<uint32_t, std::vector<IDCInfo>> mIDCInfos;  ///< follow if both end points for each CRU was seen
  std::unique_ptr<o2::utils::TreeStreamRedirector> mDebugStream; ///< debug output streamer

  //____________________________________________________________________________
  void snapshotIDCs(DataAllocator& output)
  {
    LOGP(info, "snapshotIDCs");
    // received data from all expected CRUs?
    //if (mIDCInfos.size() != mCRUs.size()) {
    //LOGP(fatal, "didn't receive data from all expected CRUs {} != {}", mIDCInfos.size(), mCRUs.size());
    //}

    // check integrety of data between CRUs
    size_t lastSize = 0;
    std::vector<IDCInfo> const* infVecComp = nullptr;
    std::vector<uint64_t> orbitBCInfo;

    for (const auto& [cru, infVec] : mIDCInfos) {

      for (const auto& inf : infVec) {
        if (!inf.hasBothEPs()) {
          LOGP(fatal, "IDC CRU {:3}: data missing at ({:8}, {:4}) for one or both end points {:02b}", cru, inf.heartbeatOrbit, inf.heartbeatBC, inf.epSeen);
        }
      }

      if (!infVecComp) {
        infVecComp = &infVec;
        lastSize = infVec.size();
        std::for_each(infVec.begin(), infVec.end(), [&orbitBCInfo](const auto& inf) { orbitBCInfo.emplace_back((uint64_t(inf.heartbeatOrbit) << 32) + uint64_t(inf.heartbeatBC)); });
        continue;
      }

      if (lastSize != infVec.size()) {
        LOGP(fatal, "IDC CRU {:3}: unequal number of IDC values {} != {}", cru, lastSize, infVec.size());
      }

      if (!std::equal(infVecComp->begin(), infVecComp->end(), infVec.begin())) {
        LOGP(fatal, "IDC CRU {:3}: mismatch in orbits");
      }
    }

    // send data
    for (const auto& [cru, idcVec] : mIDCvectors) {
      const header::DataHeader::SubSpecificationType subSpec{cru << 7};
      output.snapshot(Output{gDataOriginTPC, "IDCVECTOR", subSpec}, idcVec);
      output.snapshot(Output{gDataOriginTPC, "IDCORBITS", subSpec}, orbitBCInfo);
    }

    // clear output
    initIDC();
  }

  //____________________________________________________________________________
  void initIDC()
  {
    for (const auto cruID : mCRUs) {
      const CRU cru(cruID);
      const GlobalPadNumber numberPads = PADSPERREGION[cru.region()] * mMaxIDCPerTF; // ToDo: Add Mapper:: once in mapper
      auto& idcVec = mIDCvectors[cruID];
      idcVec.resize(numberPads);
      std::fill(idcVec.begin(), idcVec.end(), -1.f);

      auto& infosCRU = mIDCInfos[cruID];
      infosCRU.clear();
      //for (auto& info : infosCRU) {
      //info = IDCInfo{};
      //}
    }
  }

  //____________________________________________________________________________
  void writeDebugOutput()
  {
    mDebugStream->GetFile()->cd();
    auto& stream = (*mDebugStream) << "idcs";
    uint32_t seen = 0;
    for (auto cru : mCRUs) {
      if (mIDCInfos.find(cru) == mIDCInfos.end()) {
        continue;
      }
      auto& infos = mIDCInfos[cru];
      auto& idcVec = mIDCvectors[cru];

      for (int i = 0; i < infos.size(); ++i) {
        auto& info = infos[i];

        auto idcFirst = idcVec.begin() + i * PADSPERREGION[cru];
        auto idcLast = idcFirst + PADSPERREGION[cru];
        std::vector<float> idcs(idcFirst, idcLast);
        float mean = std::accumulate(idcs.begin(), idcs.end(), 0.f) / float(idcs.size());

        stream << "cru=" << cru
               << "epSeen=" << info.epSeen
               << "orbit=" << info.heartbeatOrbit
               << "bc=" << info.heartbeatBC
               << "idcs=" << idcs
               << "idc_mean=" << mean
               << "\n";
      }
    }
    //stream.getTree().Write();
  }
};

o2::framework::DataProcessorSpec getIDCToVectorSpec(const std::string inputSpec, std::vector<uint32_t> const& crus)
{
  using device = o2::tpc::IDCToVectorDevice;

  std::vector<OutputSpec> outputs;
  for (const uint32_t cru : crus) {
    const header::DataHeader::SubSpecificationType subSpec{cru << 7};
    outputs.emplace_back(gDataOriginTPC, "IDCVECTOR", subSpec, Lifetime::Timeframe);
    outputs.emplace_back(gDataOriginTPC, "IDCORBITS", subSpec, Lifetime::Timeframe);
  }

  return DataProcessorSpec{
    fmt::format("tpc-idc-to-vector"),
    select(inputSpec.data()),
    outputs,
    AlgorithmSpec{adaptFromTask<device>(crus)},
    Options{
      {"write-debug", VariantType::Bool, false, {"write a debug output tree."}},
    } // end Options
  };  // end DataProcessorSpec
}
} // namespace o2::tpc
