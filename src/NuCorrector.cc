#include "NuCorrector.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Cluster.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#endif // MARLIN_USE_AIDA

#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

using namespace lcio ;
using namespace marlin ;


NuCorrector aNuCorrector ;

NuCorrector::NuCorrector() :

    Processor("NuCorrector"),
    m_mcEnergyENu(0.f),
    m_mcEnergyELep(0.f),
	m_nSLDecayTotal(0),
	m_nSLDecayBHad(0),
	m_nSLDecayCHad(0),
    m_pTTree(NULL),
    m_hmcEnergyENu(NULL),
	m_hmcEnergyELep(NULL),
    m_pTFile(NULL)

{

//      modify processor description
    _description = "NuCorrector finds semi-leptonic decays inside jets" ;

//      register steering parameters: name, description, class-variable, default value
	registerInputCollection( LCIO::MCPARTICLE,
							"MCParticleCollection" ,
							"Name of the MCParticle collection"  ,
							m_mcParticleCollection,
                            std::string("MCParticle")
                            );
    registerProcessorParameter("RootFile",
							"Name of the output root file",
							m_rootFile,
							std::string("NuCorrector.root")
							);

}

void NuCorrector::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    // usually a good idea to
    printParameters() ;

    m_nRun = 0 ;
    m_nEvt = 0 ;
    m_nRunSum = 0;
    m_nEvtSum = 0;
    this->Clear();

	m_pTFile = new TFile(m_rootFile.c_str(), "recreate");
	m_pTTree = new TTree("PfoAnalysisTree", "PfoAnalysisTree");
	m_pTTree->SetDirectory(m_pTFile);
	m_pTTree->Branch("run", &m_nRun, "run/I");
	m_pTTree->Branch("event", &m_nEvt, "event/I");
    m_pTTree->Branch("mcEnergyENu", &m_mcEnergyENu, "mcEnergyENu/F");
    m_pTTree->Branch("mcEnergyELep", &m_mcEnergyELep, "mcEnergyELep/F");
	m_pTTree->Branch("nBHadSLDecay", &m_nSLDecayBHad, "nBHadSLDecay/I");
	m_pTTree->Branch("nCHadSLDecay", &m_nSLDecayCHad, "nCHadSLDecay/I");
	m_pTTree->Branch("nSLDecayTotal", &m_nSLDecayTotal, "nSLDecayTotal/I");
	m_hmcEnergyENu = new TH1F("mcEnu", "energy of neutrino (from SLD)", 800, 0., 800.);
	m_hmcEnergyELep = new TH1F("mcElep", "energy of lepton (from SLD)", 800, 0., 800.);

	m_hmcEnergyENu->SetDirectory(m_pTFile);
	m_hmcEnergyELep->SetDirectory(m_pTFile);

}

void NuCorrector::processRunHeader( LCRunHeader *pLCRunHeader)
{
	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;
//	m_nRun++ ;
}

void NuCorrector::processEvent( LCEvent *pLCEvent )
{

	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	++m_nEvtSum;

	if ((m_nEvtSum % 100) == 0)
		std::cout << " processed events: " << m_nEvtSum << std::endl;

	this->Clear();
	this->ExtractCollections(pLCEvent);
	this->FindEnu(pLCEvent);
	m_pTTree->Fill();
}

void NuCorrector::check( LCEvent *pLCEvent )
{
//	nothing to check here - could be used to fill checkplots in reconstruction processor
}

void NuCorrector::end()
{

	m_pTFile->cd();
	m_pTTree->Write();
	m_hmcEnergyENu->Write();
	m_hmcEnergyELep->Write();

	m_pTFile->Close();
    delete m_pTFile;

//	std::cout << "NuCorrector::end()  " << name()
//	<< " processed " << _nEvt << " events in " << _nRun << " runs "
//	<< std::endl ;

}

void NuCorrector::Clear()
{
	m_mcEnergyENu = 0.f;
    m_mcEnergyELep = 0.f;
	m_nSLDecayTotal = 0;
	m_nSLDecayBHad = 0;
	m_nSLDecayCHad = 0;

}

void NuCorrector::ExtractCollections(EVENT::LCEvent *pLCEvent)
{
	try
	{
		const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(m_mcParticleCollection);
	}
	catch (...)
	{
		streamlog_out(WARNING) << "Could not extract mc particle collection " << m_mcParticleCollection << std::endl;
	}
}
void NuCorrector::FindEnu(EVENT::LCEvent *pLCEvent)
{
	try
	{
		const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(m_mcParticleCollection);

		for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
		{
			const EVENT::MCParticle *pMCParticle = dynamic_cast<EVENT::MCParticle*>(pLCCollection->getElementAt(i));

			if (NULL == pMCParticle)
				throw EVENT::Exception("Collection type mismatch");

			const int absPdgCode(std::abs(pMCParticle->getPDG()));
			if ((absPdgCode == 12) || (absPdgCode == 14) || (absPdgCode == 16))
			{
				if (!(pMCParticle->getParents().empty()))
				{
					for (unsigned int p = 0; p < (pMCParticle->getParents()).size(); ++p)
					{
						const int absParPdgCode(std::abs((pMCParticle->getParents())[p]->getPDG()));
						const EVENT::MCParticle *parentMCParticle = dynamic_cast<EVENT::MCParticle*>((pMCParticle->getParents())[p]);
						if ((floor(absParPdgCode/100)==5) || (floor(absParPdgCode/1000)==5))
						{
							for (unsigned int d = 0; d < (((pMCParticle->getParents())[p])->getDaughters()).size(); ++d)
							{
								const int absDauPdgCode(std::abs(((pMCParticle->getParents())[p]->getDaughters())[d]->getPDG()));
								const EVENT::MCParticle *daughterMCParticle = dynamic_cast<EVENT::MCParticle*>((parentMCParticle->getDaughters())[d]);
								if ((absDauPdgCode == 11) || (absDauPdgCode == 13) || (absDauPdgCode == 15))
								{
									++m_nSLDecayTotal;
									++m_nSLDecayBHad;
									m_mcEnergyENu += pMCParticle->getEnergy();
									m_mcEnergyELep += daughterMCParticle->getEnergy();
//									m_mcEnergyELep += ((pMCParticle->getParents())[p]->getDaughters())[d]->getEnergy();
									streamlog_out(DEBUG) << "One Semi-Leptonic decay of B-Hadron was found" << std::endl;
								}
							}
						}
						if ((floor(absParPdgCode/100)==4) || (floor(absParPdgCode/1000)==4))
						{
							for (unsigned int d = 0; d < (((pMCParticle->getParents())[p])->getDaughters()).size(); ++d)
							{
								const int absDauPdgCode(std::abs(((pMCParticle->getParents())[p]->getDaughters())[d]->getPDG()));
								const EVENT::MCParticle *daughterMCParticle = dynamic_cast<EVENT::MCParticle*>((parentMCParticle->getDaughters())[d]);
								if ((absDauPdgCode == 11) || (absDauPdgCode == 13) || (absDauPdgCode == 15))
								{
									++m_nSLDecayTotal;
									++m_nSLDecayCHad;
									m_mcEnergyENu = pMCParticle->getEnergy();
									m_mcEnergyELep += daughterMCParticle->getEnergy();
//									m_mcEnergyELep = ((pMCParticle->getParents())[p]->getDaughters())[d]->getEnergy();
									streamlog_out(DEBUG) << "One Semi-Leptonic decay of Charmed-Hadron was found" << std::endl;
								}
							}
						}
					}
				}
			}
		}
		streamlog_out(DEBUG) << "Number of Semi-Leptonic decay of B-Hadron: " << m_nSLDecayBHad << std::endl;
		streamlog_out(DEBUG) << "Number of Semi-Leptonic decay of C-Hadron: " << m_nSLDecayCHad << std::endl;
		streamlog_out(DEBUG) << "Total Number of Semi-Leptonic decays: " << m_nSLDecayTotal << std::endl;
	}
	catch (...)
     {
         streamlog_out(WARNING) << "Could not extract Semi-Leptonic decay " << std::endl;
     }
	 m_hmcEnergyENu->Fill(m_mcEnergyENu, 1.);
	 m_hmcEnergyELep->Fill(m_mcEnergyELep, 1.);
}
