void dumpDigits(TString fileInfo, TString outputFileName = "", TString pedestalFileName = "", Int_t nevents = 100, uint32_t verbosity = 0, uint32_t debugLevel = 0)
{
  using namespace o2::tpc;
  DigitDump dig;
  dig.setupContainers(fileInfo, verbosity, debugLevel);
  dig.setDigitFileName(outputFileName.Data());
  dig.setPedestalAndNoiseFile(pedestalFileName.Data());
  dig.addPadMask(0, 4, 42);
  dig.addPadMask(0, 28, 60);
  dig.addPadMask(0, 61, 69);

  //dig.setTimeBinRange(0, numberTimeBins);

  CalibRawBase::ProcessStatus status = CalibRawBase::ProcessStatus::Ok;
  for (Int_t i = 0; i < nevents; ++i) {
    status = dig.processEvent();
    cout << "Processing event " << i << " with status " << int(status) << '\n';
    if (status != CalibRawBase::ProcessStatus::Ok) {
      break;
    }
  }
}
