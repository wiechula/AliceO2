void runPulser(TString fileInfo, TString outputFileName = "", Int_t nevents = 100,
               Int_t adcMin = 0, Int_t adcMax = 1100,
               Int_t firstTimeBin = 0, Int_t lastTimeBin = 500,
               TString pedestalAndNoiseFile = "",
               uint32_t verbosity = 0, uint32_t debugLevel = 0)
{
  using namespace o2::tpc;
  // ===| set up calibration class |============================================
  CalibPulser calib;
  calib.setADCRange(adcMin, adcMax);
  calib.setupContainers(fileInfo, verbosity, debugLevel);
  calib.setTimeBinRange(firstTimeBin, lastTimeBin);
  calib.setDebugLevel();
  //calib.setDebugLevel(debugLevel);

  // ===| load pedestal if requested |==========================================
  if (!pedestalAndNoiseFile.IsNull()) {
    CalDet<float> dummy;
    CalDet<float>* pedestal = nullptr;
    CalDet<float>* noise = nullptr;
    TFile f(pedestalAndNoiseFile);
    if (!f.IsOpen() || f.IsZombie()) {
      std::cout << "Could not open noise and pedestal file " << pedestalAndNoiseFile << "\n";
      return;
    }

    f.GetObject("Pedestals", pedestal);
    f.GetObject("Noise", noise);
    if (!(pedestal && noise)) {
      std::cout << "Could not load pedestal and nosie from file " << pedestalAndNoiseFile << "\n";
      return;
    }
    calib.setPedestalAndNoise(pedestal, noise);
  }

  CalibRawBase::ProcessStatus status = CalibRawBase::ProcessStatus::Ok;

  for (Int_t i = 0; i < nevents; ++i) {
    status = calib.processEvent(i);
    cout << "Processing event " << i << " with status " << int(status) << '\n';
    //if ((status != CalibRawBase::ProcessStatus::Ok) && (status != CalibRawBase::ProcessStatus::IncompleteEvent)) {
    if ((status == CalibRawBase::ProcessStatus::LastEvent)) {
      break;
    }
  }
  calib.analyse();

  cout << "Number of processed events: " << calib.getNumberOfProcessedEvents() << '\n';
  cout << "Status: " << int(status) << '\n';
  if (outputFileName.IsNull()) {
    outputFileName = "Pulser.root";
  }

  calib.dumpToFile(outputFileName.Data());

  cout << "To display the Pulsers run: root.exe $calibMacroDir/drawPulser.C'(\"" << outputFileName << "\")'\n";
}
