#ifndef NuCorrector_h
#define NuCorrector_h 1

#include "EVENT/LCStrVec.h"
#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <math.h>

#include <set>
#include <vector>

class TFile;
class TH1F;
class TTree;

using namespace lcio ;
using namespace marlin ;


/**  Example processor for marlin.
 *
 *  If compiled with MARLIN_USE_AIDA
 *  it creates a histogram (cloud) of the MCParticle energies.
 *
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4>
 *  A histogram.
 *
 * @param CollectionName Name of the MCParticle collection
 *
 * @author F. Gaede, DESY
 * @version $Id: NuCorrector.h,v 1.4 2005-10-11 12:57:39 gaede Exp $
 */

class NuCorrector : public Processor {

 public:

  virtual Processor*  newProcessor() { return new NuCorrector ; }


  NuCorrector() ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   */
  virtual void processRunHeader( lcio::LCRunHeader *pLCRunHeader ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( EVENT::LCEvent *pLCEvent ) ;


  virtual void check( EVENT::LCEvent *pLCEvent ) ;

private:

  void Clear();

  void ExtractCollections(EVENT::LCEvent *pLCEvent);

  void FindEnu(EVENT::LCEvent *pLCEvent);

  /** Called after data processing for clean up.
   */
  virtual void end() ;


 protected:

  /** Input collection name.
   */
  std::string _colName{};

	int			m_nRun{} ;
	int			m_nEvt{} ;
	int			m_nRunSum{};
	int			m_nEvtSum{};

	int				m_nSLDecayTotal;                           ///<
    int				m_nSLDecayBHad;                           ///<
    int				m_nSLDecayCHad;                           ///<
	float				m_mcEnergyENu;                           ///<
    float				m_mcEnergyELep;                           ///<

	std::string		m_mcParticleCollection{};

//	TFile              *m_pTFile{};                             ///<
//    TTree              *m_pTTree{};                             ///<
    TH1F               *m_hmcEnergyENu{};                      ///<
	TH1F               *m_hmcEnergyELep{};                      ///<

	int			m_printing;
	std::string		m_rootFile{};                           ///<
	int			m_lookForQuarksWithMotherZ;

	TFile			*m_pTFile{};                             ///<
	TTree			*m_pTTree{};                             ///<

} ;

#endif
