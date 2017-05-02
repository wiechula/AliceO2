void runPedestal(TString fileInfo="GBTx0_Battery_Floating:0:0;GBTx1_Battery_Floating:1:0")
{
  using namespace o2::TPC;
  CalibPedestal ped;
  ped.setupContainers(fileInfo);

  //while (ped.ProcessEvent());
  for (Int_t i=0; i<5; ++i) {
    ped.ProcessEvent();
  }
  ped.analyse();

  ped.dumpToFile("testPedestal.root");

  const CalDet<float>& calPedestal = ped.getPedestal();
  const CalDet<float>& calNoise    = ped.getNoise();
  Painter::Draw(calPedestal);
  Painter::Draw(calNoise);
}
