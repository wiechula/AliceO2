// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file findKrBoxCluster.C
/// \brief This macro retrieves clusters from Krypton and X-Ray runs, input tpcdigits.root
/// \author Philip Hauer <philip.hauer@cern.ch>

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"

#include "TPCReconstruction/KrCluster.h"
#include "TPCReconstruction/KrBoxClusterFinder.h"
#include "TPCBase/Digit.h"

#include <array>
#include <iostream>
#include <tuple>
#include <vector>
#endif

void findKrBoxCluster()
{
  // Read the digits:
  TFile* file = new TFile("tpcdigits.root");
  TTree* tree = (TTree*)file->Get("o2sim");
  Long64_t nEntries = tree->GetEntries();
  std::cout << "The Tree has " << nEntries << " Entries." << std::endl;

  // Initialize File for later writing
  TFile* f = new TFile("boxClustersDefault.root", "RECREATE", "Clusters");
  // Initialize Tree which will later be put into the root file
  // Maybe it is better to save each sector in a seperate tree
  TTree* T = new TTree("T", "Clusters");
  // Tree will be filled with a vector of clusters
  std::vector<o2::tpc::KrCluster> vCluster{};
  T->Branch("cluster", &vCluster);

  // Probably an unnecessary complicated way of reading the evaluates
  // Since the X-Ray runs were performed during pre-commissioning, we could
  // only readout two sectors at once. In addition, the sectors changed for
  // different runs. Hence, this way:
  std::array<std::vector<o2::tpc::Digit>*, 2> DigitizedSignal; // this arry will hold the data of the tpc digits for 2 sectors
  std::array<std::vector<o2::tpc::Digit>*, 2> DigitizedSignal1;

  DigitizedSignal1[0] = nullptr;
  DigitizedSignal1[1] = nullptr;
  DigitizedSignal[0] = nullptr;
  DigitizedSignal[1] = nullptr;

  // this loop should find out which digits are non empty:
  string digit1 = "TPCDigit_1";
  string digit2 = "TPCDigit_1";
  for (int i = 0; i <= 35; i++) {
    string number = std::to_string(i);
    number = "TPCDigit_" + number;
    tree->SetBranchAddress(number.c_str(), &DigitizedSignal1[0]);
    tree->GetEntry(0);
    if (DigitizedSignal1[0]->size() != 0) {
      if (digit1 == "TPCDigit_1") {
        digit1 = number;
      } else {
        digit2 = number;
        break;
      }
    }
  }
  tree->SetBranchAddress(digit1.c_str(), &DigitizedSignal[0]);
  tree->SetBranchAddress(digit2.c_str(), &DigitizedSignal[1]);
  std::cout << digit1.c_str() << std::endl;
  std::cout << digit2.c_str() << std::endl;

  // Now everything can get processed
  // Loop over all events
  for (int iEvent = 0; iEvent < nEntries; ++iEvent) {
    std::cout << iEvent + 1 << "/" << nEntries << std::endl;
    tree->GetEntry(iEvent);
    // Each event consists of sectors (atm only two)
    for (auto sector : DigitizedSignal) {
      if (sector->size() != 0) {
        // Initialize boxCluster Object
        o2::tpc::KrBoxClusterFinder cluster(*sector);
        // Find Cluster centers
        std::vector<std::tuple<int, int, int>> localMaxima = cluster.findLocalMaxima();
        // Loop over cluster centers
        for (const std::tuple<int, int, int>& coords : localMaxima) {
          int padMax = std::get<0>(coords);
          int rowMax = std::get<1>(coords);
          int timeMax = std::get<2>(coords);
          // Build total cluster
          vCluster.push_back(cluster.buildCluster(padMax, rowMax, timeMax));
        }
        // Fill Tree
        T->Fill();
        vCluster.clear();
      }
    }
  }
  // Write Tree to file
  f->cd();
  T->Write();
  f->Close();
  return;
}
