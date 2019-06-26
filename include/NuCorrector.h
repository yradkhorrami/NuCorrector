#ifndef NuCorrector_h
#define NuCorrector_h 1

#include "EVENT/LCStrVec.h"
#include "marlin/Processor.h"
#include "IMPL/LCCollectionVec.h"
#include "lcio.h"
#include "TLorentzVector.h"
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

  void FormTLV(EVENT::LCEvent *pLCEvent);
  void CalculateNeutrinoEnergy(TLorentzVector mcHadron_tlv,TLorentzVector mcVisible_tlv,TLorentzVector mcNeutrino_tlv,TLorentzVector mcLepton_tlv);

  /** Called after data processing for clean up.
   */
  virtual void end() ;


 protected:

  /** Input collection name.
   */
	int					m_nRun{};
	int					m_nEvt{} ;
	int					m_nRunSum{};
	int					m_nEvtSum{};

	int					m_nSLDecayBHad;
	int					m_nSLDecayCHad;
	int					m_nSLDecayTotal;
	float					m_recENuPlus;
	float					m_recENuMinus;
	float					m_recENuClose;

	typedef std::vector<int>		IntVector;
	IntVector				m_BHadronIndex;
	IntVector				m_CHadronIndex;
	typedef std::vector<double>		DoubleVector;
	DoubleVector				m_mcEnergyENu;
	DoubleVector				m_mcEnergyELep;
	typedef std::vector<float>		FloatVector;
	FloatVector				m_recEnergyENuPlus;
	FloatVector				m_recEnergyENuMinus;
	FloatVector				m_recEnergyENuClose;
	float					m_EnergyCOM;

	std::string				m_mcParticleCollection{};
	std::string				m_SLDecaysCollection{};
	std::string				m_NuCorrCollection{};
	LCCollectionVec				*m_col_NuCorr{};

	TH1F					*m_hmcEnergyENu{};
	TH1F					*m_hmcEnergyELep{};

	int					m_printing;
	std::string				m_rootFile{};
	int					m_lookForQuarksWithMotherZ;

	TFile					*m_pTFile{};
	TTree					*m_pTTree{};

} ;

#endif
