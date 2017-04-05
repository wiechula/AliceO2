/// \file AliTPCUpgradeHwClusterer.cxx
/// \brief Hwclusterer for the TPC

#include "TPCSimulation/HwCluster.h"
#include "TPCSimulation/HwClusterer.h"
#include "TPCSimulation/HwClusterFinder.h"
#include "TPCSimulation/DigitMC.h"
#include "TPCSimulation/ClusterContainer.h"

#include "FairLogger.h"
#include "TMath.h"
#include "TError.h"   // for R__ASSERT()
#include "TClonesArray.h"
#include <vector>
#include <thread>
#include <mutex>

std::mutex g_display_mutex;

using namespace AliceO2::TPC;

HwClusterer::HwClusterer()
  : HwClusterer(Processing::Parallel, 0, 360, 0, false, 8, 8, 0)
{
}

//________________________________________________________________________
HwClusterer::HwClusterer(Processing processingType, int globalTime, int cru, float minQDiff,
    bool assignChargeUnique, int padsPerCF, int timebinsPerCF, int cfPerRow)
  : Clusterer()
  , mProcessingType(processingType)
  , mGlobalTime(globalTime)
  , mCRUs(cru)
  , mMinQDiff(minQDiff)
  , mAssignChargeUnique(assignChargeUnique)
  , mPadsPerCF(padsPerCF)
  , mTimebinsPerCF(timebinsPerCF)
  , mCfPerRow(cfPerRow)
{
}

//________________________________________________________________________
HwClusterer::~HwClusterer()
{
  LOG(DEBUG) << "Enter Destructor of HwClusterer" << FairLogger::endl;

  delete mClusterContainer;

  for (auto it : mClusterFinder) {
    for (auto itt : it) {
      for (auto ittt : itt) {
        delete ittt;
      }
    }
  }
}

//________________________________________________________________________
void HwClusterer::Init()
{
  LOG(DEBUG) << "Enter Initializer of HwClusterer" << FairLogger::endl;     
  // Test that init was not called before
  R__ASSERT(!mClusterContainer);


  /*
   * initalize all cluster finder
   */
  mCfPerRow = (int)ceil((Double_t)(mPadsMax+2+2)/(mPadsPerCF-2-2));
  mClusterFinder.resize(mCRUs);
  for (int iCRU = 0; iCRU < mCRUs; iCRU++){
    mClusterFinder[iCRU].resize(mRowsMax);
    for (int iRow = 0; iRow < mRowsMax; iRow++){
      mClusterFinder[iCRU][iRow].resize(mCfPerRow);
      for (int iCF = 0; iCF < mCfPerRow; iCF++){
        int padOffset = iCF*(mPadsPerCF-2-2)-2;
        mClusterFinder[iCRU][iRow][iCF] = new HwClusterFinder(iCRU,iRow,iCF,padOffset,mPadsPerCF,mTimebinsPerCF,mMinQDiff,mMinQMax,mRequirePositiveCharge);
        mClusterFinder[iCRU][iRow][iCF]->setAssignChargeUnique(mAssignChargeUnique);


        /*
         * Connect always two CFs to be able to comunicate found clusters. So
         * the "right" one can tell the one "on the left" which pads were
         * already used for a cluster.
         */
        if (iCF != 0) {
          mClusterFinder[iCRU][iRow][iCF]->setNextCF(mClusterFinder[iCRU][iRow][iCF-1]);
        }
      }
    }
  }


  /*
   * vector of HwCluster vectors, one vector for each CRU (possible thread)
   * to store the clusters found there
   */
  mClusterStorage.resize(mCRUs);


  /* 
   * vector of digit vectors, one vector for each CRU (possible thread) to
   * store there only those digits which are relevant for this particular 
   * CRU (thread)
   */
  mDigitContainer.resize(mCRUs);
  for (std::vector<std::vector<DigitMC*>>& dc : mDigitContainer ) dc.resize(mRowsMax);
  

  mClusterContainer = new ClusterContainer();
  mClusterContainer->InitArray("AliceO2::TPC::HwCluster");
}

//________________________________________________________________________
void HwClusterer::processDigits(
    const std::vector<std::vector<DigitMC*>>& digits,
    const std::vector<std::vector<HwClusterFinder*>>& clusterFinder, 
          std::vector<HwCluster>& cluster,
          int iCRU,
          int maxRows,
          int maxPads, 
          int maxTime)
{
//  std::thread::id this_id = std::this_thread::get_id();
//  g_display_mutex.lock();
//  std::cout << "thread " << this_id << " started.\n";
//  g_display_mutex.unlock();

  /*
   * prepare local storage
   */
  float iAllBins[maxTime][maxPads];
  for (int t = 0; t < maxTime; ++t) {
    for (int p = 0; p < maxPads; ++p) {
      iAllBins[t][p] = 0.0;
    }
  }
//  float iAllBins[maxTime][maxPads];
//  for (float (&timebin)[maxPads] : iAllBins) {
//    for (float &pad : timebin) {
//        pad = 0.0;
//    }
//  }

  for (int iRow = 0; iRow < maxRows; iRow++){

    /*
     * fill in digits
     */
    for (std::vector<DigitMC*>::const_iterator it = digits[iRow].begin(); it != digits[iRow].end(); ++it){
      const Int_t iTime         = (*it)->getTimeStamp();
      const Int_t iPad          = (*it)->getPad() + 2;  // offset to have 2 empty pads on the "left side"
      const Float_t charge      = (*it)->getCharge();
  
      //iAllBins[iTime][iPad] = charge;
      
      iAllBins[iTime][iPad] += charge;
    }
  
    /*
     * copy data to cluster finders
     */
    const Short_t iPadsPerCF = clusterFinder[iRow][0]->getNpads();
    const Short_t iTimebinsPerCF = clusterFinder[iRow][0]->getNtimebins();
    std::vector<std::vector<HwClusterFinder*>::const_reverse_iterator> cfWithCluster;
    for (int time = 0; time < maxTime; time++){
      for (int pad = 0; pad < maxPads; pad = pad + (iPadsPerCF -2 -2 )) {
        const Short_t cf = pad / (iPadsPerCF-2-2);
        clusterFinder[iRow][cf]->AddTimebin(&iAllBins[time][pad],time,(maxPads-pad)>=iPadsPerCF?iPadsPerCF:(maxPads-pad));
      }
      
      /*
       * search for clusters and store reference to CF if one was found
       */
      if (time % (iTimebinsPerCF-2 -2) == 0)  {

        /*  
         * ordering is important: from right to left, so that the CFs could inform each other if cluster was found
         */
        for (auto rit = clusterFinder[iRow].crbegin(); rit != clusterFinder[iRow].crend(); ++rit) {
          if ((*rit)->findCluster()) {
            cfWithCluster.push_back(rit);
          }
        }
      }
    }
  
    /*  
     * collect found cluster
     */
    for (std::vector<HwClusterFinder*>::const_reverse_iterator &cf_it : cfWithCluster) {
      std::vector<HwCluster>* cc = (*cf_it)->getClusterContainer();
      for (HwCluster& c : *cc){
        cluster.push_back(c);
      }
      (*cf_it)->clearClusterContainer();
    }

    /* 
     * remove digits again from storage
     */
    for (std::vector<DigitMC*>::const_iterator it = digits[iRow].begin(); it != digits[iRow].end(); ++it){
      const Int_t iTime       = (*it)->getTimeStamp();
      const Int_t iPad        = (*it)->getPad() + 2;  // offset to have 2 empty pads on the "left side"
  
      iAllBins[iTime][iPad] = 0.0;
    }
  }

//  g_display_mutex.lock();
//  std::cout << "thread " << this_id << " finished.\n";
//  g_display_mutex.unlock();
}

//________________________________________________________________________
ClusterContainer* HwClusterer::Process(TClonesArray *digits)
{
  R__ASSERT(mClusterContainer);
  mClusterContainer->Reset();


  /*
   * vector to store all threads for parallel processing
   * one thread per CRU (360 in total)
   */
  std::vector<std::thread> thread_vector;

  /*  
   * clear old storages
   */
  for (std::vector<HwCluster>& cs : mClusterStorage) cs.clear();
  for (std::vector<std::vector<DigitMC*>>& dc : mDigitContainer ) {
              for (std::vector<DigitMC*>& dcc : dc) dcc.clear();
  }

  int cru = 0;
  int iRow;
  int iPad;
  int iTimeBin;
  float charge;

  if (mProcessingType == Processing::Parallel)
    LOG(DEBUG) << std::thread::hardware_concurrency() << " concurrent threads are supported." << FairLogger::endl;

  /*  
   * Loop over digits
   */
  for (TIter digititer = TIter(digits).Begin(); digititer != TIter::End(); ++digititer) {
    DigitMC* digit = dynamic_cast<DigitMC*>(*digititer);

    /*
     * if CRU number of current digit changes, start processing (either
     * sequential or parallel) of all CRUs in between the last processed
     * one and the current one.
     */
    if (cru != digit->getCRU()) {
      for (int iCRU = cru; iCRU < digit->getCRU(); iCRU++) {
//        std::cout << "starting CRU " << iCRU << std::endl;
        if (mProcessingType == Processing::Parallel)
          thread_vector.push_back(
              std::thread(
                processDigits,                      // function name
                std::ref(mDigitContainer[iCRU]),    // digit container for individual CRUs
                std::ref(mClusterFinder[iCRU]),     // cluster finder for individual CRUs
                std::ref(mClusterStorage[iCRU]),    // container to store found clusters
                iCRU,                               // current CRU for deb. purposes
                mRowsMax,                           // max. numbers of rows per CRU
                mPadsMax+2+2,                       // max. numbers of pads in each row (+2 empty ones on each side)
                mTimeBinsMax                        // number of timebins to process
              )
            );
        else {
          processDigits(mDigitContainer[iCRU],mClusterFinder[iCRU],mClusterStorage[iCRU],iCRU,mRowsMax,mPadsMax+2+2,mTimeBinsMax);
        }
      }
    }



    /*
     * add current digit to storage
     */
    cru      = digit->getCRU();
    iRow     = digit->getRow();
    iPad     = digit->getPad();
    iTimeBin = digit->getTimeStamp();
    charge   = digit->getCharge();
//    if ((cru == 179)) {// && iRow == 5)){
//      printf("hw:  digi: %d, %d, %d, %d, %.2f\n", cru, iRow, iPad, iTimeBin, charge);
//    }
    mDigitContainer[digit->getCRU()][digit->getRow()].push_back(digit);
  }


  /*
   * process all remaining CRUs
   */
  for (int iCRU = cru; iCRU < mCRUs; iCRU++) {
//    std::cout << "starting CRU " << iCRU << " at the end" << std::endl;
    if (mProcessingType == Processing::Parallel)
      thread_vector.push_back(
          std::thread(
            processDigits,                      // function name
            std::ref(mDigitContainer[iCRU]),    // digit container for individual CRUs
            std::ref(mClusterFinder[iCRU]),     // cluster finder for individual CRUs
            std::ref(mClusterStorage[iCRU]),    // container to store found clusters
            iCRU,                               // current CRU for deb. purposes
            mRowsMax,                           // max. numbers of rows per CRU
            mPadsMax+2+2,                       // max. numbers of pads in each row (+2 empty ones on each side)
            mTimeBinsMax                        // number of timebins to process
          )
        );
    else {
      processDigits(mDigitContainer[iCRU],mClusterFinder[iCRU],mClusterStorage[iCRU],iCRU,mRowsMax,mPadsMax+2+2,mTimeBinsMax);
    }
  }


  /*
   * wait for threads to join
   */
  for (std::thread& t: thread_vector) {
    t.join();
  }

  /*
   * collect clusters from individual cluster finder
   */
  for (std::vector<HwCluster> cc : mClusterStorage) {
    if (cc.size() != 0) {  
      for (HwCluster& c : cc){
        mClusterContainer->AddCluster(c.getCRU(),c.getRow(),c.getQ(),c.getQmax(),
            c.getPadMean(),c.getTimeMean(),c.getPadSigma(),c.getTimeSigma());
////          if (c.getPadSigma() > 0.45 && c.getPadSigma() < 0.5) {
//          if ((c.getCRU() == 179)){// && c.getRow() == 5)){// && (int)c.getPadMean() == 103 && (int)c.getTimeMean() == 170) || 
////              (iCRU == 256 && iRow == 10 && (int)c.getPadMean() == 27 && (int)c.getTimeMean() == 181) ) { 
//           std::cout << "HwCluster - ";
////           c.Print(std::cout);
//           c.PrintDetails(std::cout);
//           std::cout << std::endl;
//          }
//          c.Print();
      }
    }
  }


  return mClusterContainer;
}
      
