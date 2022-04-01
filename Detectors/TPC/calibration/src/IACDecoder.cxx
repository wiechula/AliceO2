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

#include <cassert>
#include <string_view>

#include "Framework/Logger.h"

#include "TPCCalibration/IACDecoder.h"

using namespace o2::tpc::iac;

//______________________________________________________________________________
bool Decoder::process(const char* data, size_t size)
{
  assert(size == sizeof(iac::packet));
  auto& iac = *(iac::packet*)data;
  const auto instance = iac.getInstance();
  if (instance >= Instances) {
    return true;
  }

  if (mDebugLevel & (uint32_t)DebugFlags::PacketInfo) {
    printPacketInfo(iac);
  }

  const auto packetInstance = iac.header.pktCount;
  const auto packetFE = iac.data.pktNumber;
  const bool isOK = iac.check();
  const auto dataWords = iac.getDataWords();
  const auto feIndex = iac.getFEIndex();
  auto& lastPacketInstance = mPktCountInstance[instance];
  auto& lastPacketFE = mPktCountFEs[feIndex];

  // check packet counters are increasing by one
  //
  if (lastPacketInstance && (packetInstance != (lastPacketInstance + 1))) {
    LOGP(error, "Packet for instance {} missing, last packet {}, this packet {}", instance, lastPacketInstance, packetInstance);
  }

  if (lastPacketFE && (packetFE != (lastPacketFE + 1))) {
    LOGP(error, "Packet for frontend {} missing, last packet {}, this packet {}", feIndex, lastPacketFE, packetFE);
  }
  lastPacketInstance = packetInstance;
  lastPacketFE = packetFE;

  if (isOK) {
    auto& dataStrings = mDataStrings[feIndex];
    dataStrings.insert(dataStrings.end(), dataWords.begin(), dataWords.end());
  } else {
    LOGP(error, "Problem in IAC data found, header check: {}, data check: {}", iac.header.check(), iac.data.check());
  }

  return isOK;
}

bool Decoder::decodeChannels(DecodedData& iacs, size_t& carry, int feid)
{
  const auto& data = mDataStrings[feid];
  const size_t dataSize = data.size();
  const size_t next = std::min(size_t(8 * 8), dataSize - carry);
  const size_t start = carry;
  while (carry < dataSize) {
    if (carry + 8 >= dataSize) {
      return false;
    }
    if (data[carry] >= '0' && data[carry] <= '7') {
      const uint32_t channel = data[carry] - '0';
      ++carry;
      uint32_t value = 0;
      for (int i = 0; i < 6; ++i) {
        const auto c = data[carry];
        uint32_t nibble = 0;
        if ((c >= '0') && (c <= '9')) {
          nibble = c - '0';
        } else if ((c >= 'A') && (c <= 'F')) {
          nibble = c - 'A' + 10;
        }
        value <<= 4;
        value |= (nibble & 0xF);
        ++carry;
      }
      int32_t valueSigned = value & 0x00FFFFFF;
      // negative value?
      if ((valueSigned >> 23) & 1) {
        valueSigned |= 0xff000000;
      }
      iacs.currents[channel] = valueSigned;

      if (data[carry] != '\n') {
        LOGP(warning, "Problem decoding data value for FE {}, channel {} at position {} / {}, dump: {}\n",
             feid, channel, carry, dataSize, std::string_view(&data[start], next));
      }
      ++carry;
    } else {
      return true;
    }
  }
  return true;
}

void Decoder::decode(int feid)
{
  auto& data = mDataStrings[feid];
  DecodedData decdata;
  DecodedData decAdditional;
  bool aligned{false};

  size_t carry = 0;
  size_t textSize = 0;
  const size_t dataSize = data.size();

  fmt::print("================ Processing feid {:2} with size {}  =================\n", feid, data.size());
  while (carry < dataSize) {
    if (!aligned) {
      while (data[carry] != '\n') {
        if (carry >= dataSize) {
          return;
        }
        ++carry;
      }
      ++carry;
      aligned = true;
    }
    const size_t next = std::min(size_t(20), dataSize - carry);
    // fmt::print("Checking position {} / {}, {}\n", carry, dataSize, std::string_view(&data[carry], next));

    if (data[carry] >= '0' && data[carry] <= '7') {
      if (!decodeChannels(decdata, carry, feid)) {
        break;
      }
    } else if (data[carry] == 'S') {
      // time stamp comes after channel data
      ++carry;
      std::string_view vd(&data[carry], 8);
      std::stringstream str(vd.data());
      str.flags(std::ios_base::hex);
      uint32_t timeStamp;
      str >> timeStamp;
      decdata.timeStamp = timeStamp;
      decAdditional.timeStamp = timeStamp;

      carry += 8;
      if (data[carry] != '\n' || data[carry + 1] != 's') {
        LOGP(warning, "Problem decoding time stamp for FE () at position {} / {}, dump: {}\n",
             feid, carry - 8, dataSize, std::string_view(&data[carry - 8], next));
      } else {
        if (mDebugStream && (mDebugLevel & (uint32_t)DebugFlags::StreamSingleFE)) {
          (*mDebugStream) << "d"
                          << "data=" << decdata
                          << "feid=" << feid
                          << "tscount=" << mTSCountFEs[feid]
                          << "\n";
        }
        ++mTSCountFEs[feid].first;

        // copy decoded data to output
        const auto refTime = timeStamp / TimeSampleDistance;
        const auto nSamples = mCurrentsTime.size();
        if ((mTSCountFEs[feid].second == 0) && (refTime > 0)) {
          LOGP(info, "Skipping initial data packet {} with time stamp {}", mTSCountFEs[feid].first, timeStamp);
        } else {
          ++mTSCountFEs[feid].second;
          if (nSamples < refTime + 1) {
            mCurrentsTime.resize(refTime + 1);
          }
          auto& currentData = mCurrentsTime[refTime];
          currentData.time = refTime; // TODO: set proper time in ms since epoch!!!
          currentData.timeStamps[feid] = timeStamp;
          auto& currents = currentData.currents;
          std::copy(decdata.currents.begin(), decdata.currents.end(), currents.begin() + feid * ChannelsPerFE);
        }
      }

      carry += 2;
      decdata.reset();
    } else if (const auto pos = mDecodeAdditional.find(data[carry]); (pos != std::string::npos) && mDebugStream) {
      // in case of debug stream output, decode additionally configured data streams
      const auto streamStart = carry;
      const char streamType = data[carry];
      const char endMarker = streamType + 32;
      ++carry;

      if (!decodeChannels(decAdditional, carry, feid)) {
        return;
      }

      if (data[carry] != endMarker) {
        const size_t next = std::min(size_t(20), dataSize - carry);
        LOGP(warning, "Problem decoding additional stream '{}' values for FE () at position {} / {}, dump: {}",
             streamType, feid, carry, dataSize, std::string_view(&data[streamStart], next));
      } else {
        const char treeName[2] = {streamType, '\0'};
        (*mDebugStream) << treeName
                        << "data=" << decAdditional
                        << "feid=" << feid
                        << "\n";
      }

      decAdditional.reset();
      ++carry;
    } else if (AllowedAdditionalStreams.find(data[carry] != std::string_view::npos)) {
      // skip stream if not configured or no debug stream
      const char streamType = data[carry];
      const char endMarker = streamType + 32;
      while (data[carry] != endMarker) {
        if (carry >= dataSize) {
          return;
        }
        ++carry;
      }
      ++carry;
    } else if (data[carry] >= 'a' && data[carry] <= 'z') {
      LOGP(info, "Skipping {}", data[carry]);
      ++carry;
      decdata.reset();
    } else {
      LOGP(error, "Can't interpret position {} / {}, {}, stopping decoding\n", carry, dataSize, std::string_view(&data[carry], next));
      break;
    }
  }
}

void Decoder::finalize()
{
  LOGP(info, "finalize");

  for (int feid = 0; feid < NumberFEs; ++feid) {
    decode(feid);
  }

  if (mDebugStream) {
    for (auto& currentData : mCurrentsTime) {
      (*mDebugStream) << "c"
                      << "values=" << currentData
                      << "\n";
    }
    mDebugStream->Close();
    mDebugStream.reset();
  }
}

void Decoder::printPacketInfo(const iac::packet& iac)
{
  const auto& header = iac.header;
  const auto& iacc = iac.data;

  LOGP(info, "{:>4} {:>4} {:>8} {:>8} -- {:>4} {:>4} {:>8} {:>8} {:>10} -- {:>4}\n", //
       "vers",                                                                       //
       "inst",                                                                       //
       "bc",                                                                         //
       "pktCnt",                                                                     //
       "feid",                                                                       //
       "size",                                                                       //
       "pktNum",                                                                     //
       "time",                                                                       //
       "crc32",                                                                      //
       "ok"                                                                          //
  );

  LOGP(info, "{:>4} {:>4} {:>8} {:>8} -- {:>4} {:>4} {:>8} {:>8} {:>#10x} -- {:>4b}\n", //
       header.version,                                                                  //
       header.instance,                                                                 //
       header.bunchCrossing,                                                            //
       header.pktCount,                                                                 //
       iacc.feid,                                                                       //
       iacc.pktSize,                                                                    //
       iacc.pktNumber,                                                                  //
       iacc.timeStamp,                                                                  //
       iacc.crc32,                                                                      //
       iacc.check()                                                                     //
  );
}
