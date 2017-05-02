#ifndef ALICEO2_TPC_CALIBRAWBASE_H_
#define ALICEO2_TPC_CALIBRAWBASE_H_

/// \file   CalibRawBase.h
/// \author Jens Wiechula, Jens.Wiechula@ikf.uni-frankfurt.de

#include <vector>
#include <memory>

#include "TString.h"
#include "Rtypes.h"

#include "TPCReconstruction/GBTFrameContainer.h"


namespace o2
{
namespace TPC
{


/// \brief Base class for raw data calibrations
///
/// This class is the base class for raw data calibrations
/// It implements base raw reader functionality and calls 
/// an 'Update' function for each digit
///
/// origin: TPC
/// \author Jens Wiechula, Jens.Wiechula@ikf.uni-frankfurt.de

class CalibRawBase
{
  public:
    CalibRawBase() : mNevents(0), mTimeBinsPerCall(500) {;}

    /// Update function called once per digit
    ///
    /// \param sector
    virtual Int_t Update(const Int_t roc, const Int_t row, const Int_t pad,
                         const Int_t timeBin, const Float_t signal) = 0;

    /// add GBT frame container to process
    void addGBTFrameContainer(GBTFrameContainer *cont) { mGBTFrameContainers.push_back(std::unique_ptr<GBTFrameContainer>(cont)); }

    /// set number of time bins to process in one call to ProcessEvent
    void setTimeBinsPerCall(Int_t nTimeBins) { mTimeBinsPerCall = nTimeBins; }

    /// return the number of time bins processed in one call to ProcessEvent
    Int_t getTimeBinsPerCall() const { return mTimeBinsPerCall; }

    void setupContainers(TString fileInfo);
    
    /// Process one event with mTimeBinsPerCall length
    Bool_t ProcessEvent(); 
    
    /// Rewind the events
    void RewindEvents();

    /// Dump the relevant data to file
    virtual void dumpToFile(TString filename) {}
    
  private:
    size_t    mNevents;                //!< number of processed events
    Int_t     mTimeBinsPerCall;        //!< numver of time bins to process in ProcessEvent
    std::vector<std::unique_ptr<GBTFrameContainer>> mGBTFrameContainers; //! raw reader pointer

    virtual void ResetEvent() = 0;
    virtual void EndEvent() {++mNevents; }
};

//----------------------------------------------------------------
// Inline Functions
//----------------------------------------------------------------
inline Bool_t CalibRawBase::ProcessEvent()
{
  if (!mGBTFrameContainers.size()) return kFALSE;
  ResetEvent();
  
  // loop over raw readers, fill digits for 500 time bins and process
  // digits
  std::vector<DigitData> digits(80);
  for (auto& reader_ptr : mGBTFrameContainers) {
    auto reader = reader_ptr.get();
    for (int i=0; i<mTimeBinsPerCall; ++i) {
      if (reader->getData(digits)) {
        for (auto& digi : digits) {
          CRU cru(digi.getCRU());
          const int sector = cru.sector().getSector();
          // TODO: OROC case needs subtraction of number of pad rows in IROC
          const int row    = digi.getRow();// + mTPCmapper.getPadRegionInfo(cru.region()).getGlobalRowOffset();
          const int pad    = digi.getPad();
          if (row==255 || pad==255) continue;
          const int timeBin= i; //digi.getTimeStamp();
          const float signal = digi.getChargeFloat();
          //const FECInfo& fecInfo = mTPCmapper.getFECInfo(PadSecPos(sector, row, pad));
          //printf("Call update: %d, %d, %d, %d (%d), %.3f -- reg: %02d -- FEC: %02d, Chip: %02d, Chn: %02d\n", sector, row, pad, timeBin, i, signal, cru.region(), fecInfo.getIndex(), fecInfo.getSampaChip(), fecInfo.getSampaChannel());
          Update(sector, row, pad, timeBin, signal );
        }
      }
      else {
        return kFALSE;
      }
      digits.clear();
    }
  }

  EndEvent();
  return kTRUE;
}

} // namespace TPC

} // namespace o2
#endif
