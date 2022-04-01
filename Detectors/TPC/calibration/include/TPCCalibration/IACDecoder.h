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

/// \file IACDecoder.h
/// \brief Decoding of integrated analogue currents
/// @author Jens Wiechula, Jens.Wiechula@ikf.uni-frankfurt.de

#ifndef ALICEO2_IACDECODER_H
#define ALICEO2_IACDECODER_H

#include <algorithm>
#include <array>
#include <deque>
#include <unordered_map>
#include <vector>

#include "DataFormatsTPC/Defs.h"
#include "DataFormatsTPC/IAC.h"
#include "CommonUtils/TreeStreamRedirector.h"

using std::size_t;

namespace o2::tpc::iac
{

constexpr float ADCtoNanoAmp = 125000.f / 8388608.f;     ///< 125000 nA / std::pow(2,23) bits
constexpr uint32_t ChannelsPerFE = 8;                    ///< Channels per front-end card. One channel is one stack
constexpr size_t Instances = 2;                          ///< Number of instances to process
constexpr size_t NumberFEs = FEsPerInstance * Instances; ///< Total number of frontends to process
                                                         ///
struct DataPoint {
  long time{};                      ///< Reference time since epoch in ms
  std::vector<uint32_t> timeStamps; ///< time stamps for each FE in internal units
  std::vector<int32_t> currents;    ///< Current in signed ADC values, use IACDecoder::ADCtoNanoAmp, to convert to nA

  DataPoint()
  {
    currents.resize(GEMSTACKS);
    timeStamps.resize(NumberFEs);
  }

  void reset()
  {
    time = 0;
    std::fill(timeStamps.begin(), timeStamps.end(), 0u);
    std::fill(currents.begin(), currents.end(), 0);
  }

  ClassDefNV(DataPoint, 1);
};

class Decoder
{
 public:
  static constexpr std::string_view AllowedAdditionalStreams{"MRLIX"}; ///< Allowed additional data streams that can be decoded with debug stream enabled
  static constexpr uint32_t TimeSampleDistance = 16;                   ///< Number of samples between time data stamps

  Decoder()
  {
    mCurrentsTime.reserve(1000 * 60); ///< reserve ~60sec
  }

  enum class DebugFlags {
    PacketInfo = 0x01,      ///< Print packe information
    DumpFullStream = 0x10,  ///< Dump the data character streams
    StreamSingleFE = 0x100, ///< Stream debug output for each single FE
  };

  bool process(const char* data, size_t size);

  /// Finalize decoding of remaining data, close debug stream
  void finalize();

  void setDebugLevel(uint32_t level = (uint32_t)DebugFlags::PacketInfo) { mDebugLevel = level; }

  void setReferenceTime(long time) { mReferenceTime = time; }

  /// Set additional data to decode
  ///
  /// These data streams will only be decoded for debug purposes and written in the debug tree, if debug tree output is enabled
  /// \param additional additional data streams (case sensitive)
  /// Parameters are 'M': mean values over some integration time
  ///                'R': RMS values
  ///                'L': low voltage data
  ///                'I': min values
  ///                'X': max values
  void setDecodeAdditional(std::string_view additional)
  {
    for (const auto c : additional) {
      if (AllowedAdditionalStreams.find(c) != std::string_view::npos) {
        mDecodeAdditional = additional;
      }
    }
  }

  void enableDebugTree()
  {
    if (!mDebugStream) {
      mDebugStream = std::make_unique<o2::utils::TreeStreamRedirector>(mDebugOutputName.data(), "recreate");
    }
  }

 private:
  /// Decoded data of one FE
  struct DecodedData {
    uint32_t timeStamp{};
    std::array<int32_t, 8> currents{};

    void resetCurrent()
    {
      std::fill(currents.begin(), currents.end(), 0);
    }

    void reset()
    {
      timeStamp = 0;
      resetCurrent();
    }
  };

  long mReferenceTime{};                                            ///< reference time when the time stamp counter in the FE was reset
  std::array<std::deque<char>, NumberFEs> mDataStrings;             ///< ASCI data sent by FE
  std::array<uint32_t, Instances> mPktCountInstance;                ///< Packet counter for the instance
  std::array<uint32_t, NumberFEs> mPktCountFEs;                     ///< Packet counter for the single FEs
  std::array<std::pair<uint32_t, uint32_t>, NumberFEs> mTSCountFEs; ///< Counter how often the time stamp was seen for the single FEs, all / valid
  std::vector<DataPoint> mCurrentsTime;                             ///< Decoded currents per Stack for decoded times
  // std::vector<uint32_t> mRefTime;                                ///< reference time stamps
  std::unique_ptr<o2::utils::TreeStreamRedirector> mDebugStream; ///< Debug output streamer
  std::string mDecodeAdditional;                                 ///< Decode these additional data for debugging purposes
  std::string mDebugOutputName{"IAC_debug.root"};                ///< name of the debug output tree

  uint32_t mDebugLevel{0}; ///< Amount of debug information to print

  bool decodeChannels(DecodedData& iacs, size_t& carry, int feid);
  void decode(int feid);

  void printPacketInfo(const iac::packet& iac);
};

} // namespace o2::tpc::iac
#endif
