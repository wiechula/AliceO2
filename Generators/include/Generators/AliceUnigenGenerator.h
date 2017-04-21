
/*! \class AliceUnigenGenerator
    \brief A class to process the Alice MC 

	Class to process the output of the Unigenerator converter applied to ALICE kinematics.root
	Adopted from CbmUnigenGenerator by D. Kresan
	\author F. Sozzi
*/

#ifndef ALICE_UNIGEN_GENERATOR_H
#define ALICE_UNIGEN_GENERATOR_H

#include "TString.h"
#include "FairGenerator.h"

class TFile;
class TTree;
class UEvent;
class FairPrimaryGenerator;


class AliceUnigenGenerator : public FairGenerator {

private:
    Int_t         fEvents;           //! Current event number
    TFile        *fInputFile;        //! Input file
    TString       fFileName;         //! Input file Name
    TTree        *fInTree;           //! Input tree
    UEvent       *fEvent;            //! Input event
    Bool_t        fCM;               //! CM flag
    Double_t      fBetaCM;           //! CM velocity
    Double_t      fGammaCM;          //! CM gamma factor

    Double32_t fPhiMin, fPhiMax;          // Limits of event plane angle
    Bool_t     fEventPlaneSet;            // Flag whether event plane angle is used

    void CloseInput();

    AliceUnigenGenerator(const AliceUnigenGenerator&);
    AliceUnigenGenerator& operator=(const AliceUnigenGenerator&);

public:
	/// Constructor
	AliceUnigenGenerator();
    /// Constructor
	/// @param fileName Name of file in Unigenerator format
	AliceUnigenGenerator(TString fileName);
    /// Destructor
	virtual ~AliceUnigenGenerator();

  /** Public method SetEventPlane 
   **@param phiMin   Lower limit for event plane angle [rad]
   **@param phiMax   Upper limit for event plane angle [rad]
   **If set, an event plane angle will be generated with flat
   **distrtibution between phiMin and phiMax. 
   **/
  void SetEventPlane(Double_t phiMin, Double_t phiMax);

  /** Reads on event from the input file and pushes the tracks onto
   ** the stack. Abstract method in base class.
   ** @param pStack    pointer to the stack
   ** @return kTRUE if successful, kFALSE if not
   **/
  virtual Bool_t ReadEvent(FairPrimaryGenerator* primGen);

  
  //void SetBoostCMtoLab()
  ClassDef(AliceUnigenGenerator,1);
};


#endif
