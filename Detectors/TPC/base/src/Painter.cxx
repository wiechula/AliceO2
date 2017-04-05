#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"

#include "TPCBase/Mapper.h"
#include "TPCBase/CalDet.h"
#include "TPCBase/CalArray.h"
#include "TPCBase/Painter.h"

using namespace o2::TPC;

//template <class T>
void Painter::Draw(CalDet<T> calDet)
{
  using DetType = CalDet<T>;
  using CalType = CalArray<T>;

  auto name = calDet.getName().c_str();

  // ===| define histograms |===================================================
  auto hAside1D = new TH1F(Form("h_Aside_1D_%s", name), Form("%s (A-Side)", name),
                         300, -100, 100); //TODO: modify ranges

  auto hCside1D = new TH1F(Form("h_Cside_1D_%s", name), Form("%s (C-Side)", name),
                         300, -100, 100); //TODO: modify ranges

  auto hAside2D = new TH2F(Form("h_Aside_2D_%s", name), Form("%s (A-Side)", name),
                         300, -300, 300, 300, -300, 300);

  auto hCside2D = new TH2F(Form("h_Cside_2D_%s", name), Form("%s (C-Side)", name),
                         300, -300, 300, 300, -300, 300);


  const Mapper& mapper = Mapper::instance();

  for (auto& cal : calDet.getData()) {

    int calPadSubsetNumber = cal.getPadSubsetNumber();
    int row = -1;
    int pad = -1;
    int offset=0;
    Sector sector;
    switch (cal.getPadSubset()) {
      case PadSubset::ROC: {
        ROC roc(calPadSubsetNumber);
        offset = (roc.rocType() == RocType::IROC)?0:mapper.getPadsInIROC();
        sector = roc.getSector();
        break;
      }
      case PadSubset::Partition: {
        const int npartitions = mapper.getNumberOfPartitions();
        const int partition = calPadSubsetNumber%npartitions;
        const int rowOffset = mapper.getPartitionInfo(partition).getGlobalRowOffset();
        offset = mapper.globalPadNumber(PadPos(rowOffset, 0));
        sector = calPadSubsetNumber/npartitions;
        break;
      }
      case PadSubset::Region: {
        const int nregions = mapper.getNumberOfPadRegions();
        const int region = calPadSubsetNumber%mapper.getNumberOfPadRegions();
        const int rowOffset = mapper.getPadRegionInfo(region).getGlobalRowOffset();
        offset = mapper.globalPadNumber(PadPos(rowOffset, 0));
        sector = calPadSubsetNumber/nregions;
        break;
      }
    }

    printf("sector: %d, %d (%d)\n", calPadSubsetNumber, sector.getSector(), sector.side()==Side::C);
    auto hist2D = hAside2D;
    if (sector.side() == Side::C) {
      hist2D = hCside2D;
    }

    int padNumber = offset;
    for (const auto& val : cal.getData()) {
      const PadCentre& localCoord = mapper.padCentre(padNumber);
      const GlobalPosition3D pos = mapper.LocalToGlobal(LocalPosition3D(localCoord.getX(), localCoord.getY(), 0), sector);

      Int_t bin = hist2D->FindBin(pos.getX(), pos.getY());

      //hist2D->SetBinContent(bin, val);
      hist2D->Fill(pos.getX(), pos.getY(), val);
      if (val>0) {
        printf("%d, %.2f, %.2f, %.2f\n", sector.getSector(), pos.getX(), pos.getY(), val);
      }

      ++padNumber;
    }
  }

  // ===| Draw histograms |=====================================================
  auto c = new TCanvas(Form("c_%s", name));
  c->Divide(2,2);

  c->cd(1);
  hAside2D->Draw("colz");

  c->cd(2);
  hCside2D->Draw("colz");

  c->cd(3);
  hAside1D->Draw();

  c->cd(4);
  hCside1D->Draw();
}

//template <class T>
void Painter::Draw(CalArray<T>)
{
}
