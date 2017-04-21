
#include "Generators/AliceUnigenGenerator.h"

#include "FairMCEventHeader.h"

#include "FairPrimaryGenerator.h"
#include "FairLogger.h"
#include "FairRunSim.h"

#include "Generators/URun.h"
#include "Generators/UEvent.h"
#include "Generators/UParticle.h"

#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"

#include <iostream>
using namespace std;


// ------------------------------------------------------------------------
AliceUnigenGenerator::AliceUnigenGenerator()
  : FairGenerator(),
    fEvents(0),
    fInputFile(NULL),
    fFileName(""),
    fInTree(NULL),
    fEvent(NULL),
    fCM(kFALSE),
    fBetaCM(0.),
    fGammaCM(1.),
    fPhiMin(0.),
    fPhiMax(0.),
    fEventPlaneSet(kFALSE)
{
}
// ------------------------------------------------------------------------


// ------------------------------------------------------------------------
AliceUnigenGenerator::AliceUnigenGenerator(TString fileName)
  : FairGenerator(),
    fEvents(0),
    fInputFile(NULL),
    fFileName(fileName),
    fInTree(NULL),
    fEvent(NULL),
    fCM(kFALSE),
    fBetaCM(0.),
    fGammaCM(0.),
    fPhiMin(0.),
    fPhiMax(0.),
    fEventPlaneSet(kFALSE)
{
  LOG(INFO) << "AliceUnigenGenerator: Opening input file "
            << fFileName.Data() << FairLogger::endl;
  fInputFile = TFile::Open(fFileName,"read");
  if (NULL == fInputFile) {
    LOG(FATAL) << "AliceUnigenGenerator: Cannot open input file." << FairLogger::endl;
}

  // Get run description
  URun* run = (URun*) fInputFile->Get("run");
  if(NULL==run) {
    LOG(FATAL) << "AliceUnigenGenerator: No run description in input file." << FairLogger::endl;
  }
  run->Print();
  //
  // Check kinematics - This is not finalized since it is not needed: The Kinematics.root is always in the LAB of Alice (correct?)
  // 
if( fCM )
  {
     
// Boost cms into LHC lab frame (taken from AliGenMC.cxx)
   fBetaCM = TMath::TanH(- 0.5 * TMath::Log(Double_t(run->GetZProj())*Double_t(run->GetATarg()) /(Double_t(run->GetAProj())*Double_t(run->GetZTarg())  )));
   fGammaCM = 1./TMath::Sqrt(1. - fBetaCM*fBetaCM);
 }
  else {
//do nothing 
	  fBetaCM = 0.;
	  fGammaCM = 1.;
   }
   cout<< "AliceUnigenGenerator: boost parameters: "
              << "betaCM = "    << fBetaCM
              << ", gammaCM = " << fGammaCM << FairLogger::endl;
  
	  
  delete run;
  fInTree = (TTree*) fInputFile->Get("events");
  if(NULL == fInTree) {
    LOG(FATAL) << "AliceUnigenGenerator: No event tree in input file." << FairLogger::endl;
  }
  fEvent = new UEvent();
  fInTree->SetBranchAddress("event", &fEvent);
}
// ------------------------------------------------------------------------


// ------------------------------------------------------------------------
AliceUnigenGenerator::~AliceUnigenGenerator()
{
  CloseInput();
}
// ------------------------------------------------------------------------


// ------------------------------------------------------------------------
Bool_t AliceUnigenGenerator::ReadEvent(FairPrimaryGenerator* primGen)
{
  // Check for input file
  if(NULL==fInputFile || NULL==fInTree) {
    LOG(ERROR) << "AliceUnigenGenerator: Input file is not opened!" << FairLogger::endl;
    return kFALSE;
  }

  // If end of input file is reached : close it and abort run
  if(fEvents >= fInTree->GetEntries()) {
    LOG(INFO) << "AliceUnigenGenerator: End of input file reached" << FairLogger::endl;
    CloseInput();
    return kFALSE;
  }

  // Read entry
  fInTree->GetEntry(fEvents);

  LOG(INFO) << "AliceUnigenGenerator: Event " << fEvent->GetEventNr()
            << ", multiplicity " << fEvent->GetNpa() << FairLogger::endl;

  UParticle* particle;
  Double_t pz;
  Double_t pz1;
  Double_t phi = 0.;

  // ---> Generate rotation angle  D
  if ( fEventPlaneSet ) { phi = gRandom->Uniform(fPhiMin, fPhiMax); }

  // Set event id and impact parameter in MCEvent if not yet done
  FairMCEventHeader* event = primGen->GetEvent();
  if ( event && (! event->IsSet()) ) {
    event->SetEventID(fEvent->GetEventNr());
    event->SetB(fEvent->GetB());
    event->SetRotZ(phi);
    event->SetNPrim(fEvent->GetNpa());
    event->MarkSet(kTRUE);
  }

  // Loop over tracks in the current event
  for (Int_t itrack = 0; itrack < fEvent->GetNpa(); itrack++) {
    // Get particle
    particle = fEvent->GetParticle(itrack);
    // Boost
    pz = particle->Pz();
    if(fCM) {
      pz1 = fGammaCM*(pz - fBetaCM*particle->E());
    } else {
      pz1 = pz;
    }
    Double_t px = particle->Px();
    Double_t py = particle->Py();
	//  LOG(DEBUG2) << "Px before: "<< px << FairLogger::endl;
	//  LOG(DEBUG2) << "Py before: "<< py << FairLogger::endl;
    // Rotate momenta by event plane angle
    if ( fEventPlaneSet ) {
      Double_t pt = TMath::Sqrt(px*px + py*py);
      Double_t azim = TMath::ATan2(py,px);
      azim += phi;
      px = pt * TMath::Cos(azim);
      py = pt * TMath::Sin(azim);
	  //   LOG(DEBUG2) << "Px after: "<< px << FairLogger::endl;
	  //   LOG(DEBUG2) << "Py after: "<< py << FairLogger::endl;
    }

                    
	Int_t pdgIndex = particle->GetPdg();
	primGen->AddTrack(pdgIndex,
                      px,
                      py,
                      pz1,
                      particle->X(), particle->Y(), particle->Z());
  }
  
  fEvents += 1;

  return kTRUE;
}

// -----   Public method SetEventPlane   ----------------------------------
void AliceUnigenGenerator::SetEventPlane(Double_t phiMin, Double_t phiMax)
{
  fPhiMin = phiMin;
  fPhiMax = phiMax;
  fEventPlaneSet = kTRUE;
}

// ------------------------------------------------------------------------
void AliceUnigenGenerator::CloseInput()
{
  if(NULL != fInputFile) {
    LOG(INFO) << "AliceUnigenGenerator: Closing input file "
              << fFileName.Data() << FairLogger::endl;
    fInputFile->Close();
    fInputFile = NULL;
  }
}
// ------------------------------------------------------------------------


ClassImp(AliceUnigenGenerator);

