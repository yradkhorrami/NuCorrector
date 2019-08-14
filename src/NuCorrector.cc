#include "NuCorrector.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"
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
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

using namespace lcio ;
using namespace marlin ;


NuCorrector aNuCorrector ;

NuCorrector::NuCorrector() :

	Processor("NuCorrector"),
	m_EnergyCOM(0.f),
	m_nSLDecayTotal(0),
	m_nSLDecayBHad(0),
	m_nSLDecayCHad(0),
	m_recENuPlus(0.f),
	m_recENuMinus(0.f),
	m_recENuClose(0.f),
	m_BHadronIndex{},
	m_CHadronIndex{},
	m_mcEnergyENu{},
	m_mcEnergyELep{},
	m_recEnergyENuPlus{},
	m_recEnergyENuMinus{},
	m_recEnergyENuClose{}

{

//      modify processor description
_description = "NuCorrector reconstructs the neytrinos produced in semi-leptonic decays inside jets" ;

//      register steering parameters: name, description, class-variable, default value
	registerInputCollection( 	LCIO::MCPARTICLE,
					"MCParticleCollection" ,
					"Name of the MCParticle collection"  ,
					m_mcParticleCollection,
					std::string("MCParticle")
					);

	registerProcessorParameter(	"Energy",
					"Center of Mass Energy",
					m_EnergyCOM,
					float(0.001f)
					);

	registerInputCollection(	LCIO::MCPARTICLE,
					"SemiLeptonicDecays",
					"Collection of semi-leptonic decays",
					m_SLDecaysCollection,
					std::string("SLDecay")
					);

	registerOutputCollection( 	LCIO::MCPARTICLE,
					"NeutrinoCorrection",
					"Collection of Corrected/Estimated neutrino energies from SemiLeptonic decays",
					m_NuCorrCollection,
					std::string("NuCorrect")
					);
}

void NuCorrector::init()
{

	streamlog_out(DEBUG) << "   init called  " << std::endl ;

//	usually a good idea to
	printParameters() ;

	m_nRun = 0 ;
	m_nEvt = 0 ;
	m_nRunSum = 0;
	m_nEvtSum = 0;
	this->Clear();

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

	m_col_NuCorr = new LCCollectionVec(LCIO::MCPARTICLE);

	this->Clear();
	this->ExtractCollections(pLCEvent);
	this->FormTLV(pLCEvent);
//	m_pTTree->Fill();

	m_col_NuCorr->parameters().setValue("nBSLD", (int)m_nSLDecayBHad);
	m_col_NuCorr->parameters().setValue("nCSLD", (int)m_nSLDecayCHad);
	m_col_NuCorr->parameters().setValue("nSLD", (int)m_nSLDecayTotal);
	m_col_NuCorr->parameters().setValue("recENuPlus", (float)m_recENuPlus);
	m_col_NuCorr->parameters().setValue("recENuMinus", (float)m_recENuMinus);
	m_col_NuCorr->parameters().setValue("recENuClose", (float)m_recENuClose);
	m_col_NuCorr->parameters().setValues("recEnergyENuPlus", (std::vector<float>)m_recEnergyENuPlus);
	m_col_NuCorr->parameters().setValues("recEnergyENuMinus", (std::vector<float>)m_recEnergyENuMinus);
	m_col_NuCorr->parameters().setValues("recEnergyENuClose", (std::vector<float>)m_recEnergyENuClose);
	pLCEvent->addCollection(m_col_NuCorr, m_NuCorrCollection);

}

void NuCorrector::check( LCEvent *pLCEvent )
{
//	nothing to check here - could be used to fill checkplots in reconstruction processor
}

void NuCorrector::end()
{

}

void NuCorrector::Clear()
{
	m_BHadronIndex.clear();
	m_CHadronIndex.clear();
	m_mcEnergyENu.clear();
	m_mcEnergyELep.clear();
	m_recEnergyENuPlus.clear();
	m_recEnergyENuMinus.clear();
	m_recEnergyENuClose.clear();
	m_nSLDecayTotal = 0;
	m_nSLDecayBHad = 0;
	m_nSLDecayCHad = 0;
	m_EnergyCOM = 0.f;
	m_recENuPlus = 0.f;
	m_recENuMinus = 0.f;
	m_recENuClose = 0.f;

}

void NuCorrector::ExtractCollections(EVENT::LCEvent *pLCEvent)
{
	LCCollection *SLDecaysCollection = NULL;
	try
	{
		const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(m_mcParticleCollection);
		SLDecaysCollection = pLCEvent->getCollection(m_SLDecaysCollection);
	}
	catch (...)
	{
		streamlog_out(WARNING) << "Could not extract mc particle collection " << m_mcParticleCollection << std::endl;
	}
}
void NuCorrector::FormTLV(EVENT::LCEvent *pLCEvent)
{
	try
	{
		const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(m_mcParticleCollection);
		const EVENT::LCCollection *SLDecaysCollection = pLCEvent->getCollection(m_SLDecaysCollection);
		m_EnergyCOM = pLCEvent->getParameters().getFloatVal("Energy");
		m_nSLDecayBHad = SLDecaysCollection->getParameters().getIntVal("nBSLD");
		m_nSLDecayCHad = SLDecaysCollection->getParameters().getIntVal("nCSLD");
		m_nSLDecayTotal = SLDecaysCollection->getParameters().getIntVal("nSLD");
		m_BHadronIndex = SLDecaysCollection->getParameters().getIntVals("BHadronIndex", m_BHadronIndex);
		m_CHadronIndex = SLDecaysCollection->getParameters().getIntVals("CHadronIndex", m_CHadronIndex);
		if (m_BHadronIndex.size()!=0)
		{
			for (unsigned int n_bsl = 0; n_bsl < m_BHadronIndex.size(); ++n_bsl)
			{
				const EVENT::MCParticle *pMCBHadron = dynamic_cast<EVENT::MCParticle*>(pLCCollection->getElementAt(m_BHadronIndex[n_bsl]));
				if (!((floor(std::abs(pMCBHadron->getPDG())/100)==5) || (floor(std::abs(pMCBHadron->getPDG())/1000)==5)))
					continue;

				TLorentzVector mcHadron_tlv = TLorentzVector(pMCBHadron->getMomentum(),pMCBHadron->getEnergy());
				TLorentzVector mcVisible_tlv(0,0,0,0);
				TLorentzVector mcNeutrino_tlv(0,0,0,0);
				TLorentzVector mcLepton_tlv(0,0,0,0);
				for (unsigned int d = 0; d < (pMCBHadron->getDaughters()).size(); ++d)
				{
					const EVENT::MCParticle *pBHadronDaughter =dynamic_cast<EVENT::MCParticle*>((pMCBHadron->getDaughters())[d]);
					const int absDauPdgCode(std::abs(pBHadronDaughter->getPDG()));
					if (absDauPdgCode == 12 || absDauPdgCode == 14 || absDauPdgCode == 16)
					{
						mcNeutrino_tlv = TLorentzVector(pBHadronDaughter->getMomentum(),pBHadronDaughter->getEnergy());
						for (unsigned int l = 0; l < (pMCBHadron->getDaughters()).size(); ++l)
						{
							if (std::abs(((pMCBHadron->getDaughters())[l])->getPDG()) == absDauPdgCode - 1)
							{
								mcLepton_tlv = TLorentzVector(((pMCBHadron->getDaughters())[l])->getMomentum(),((pMCBHadron->getDaughters())[l])->getEnergy());
							}
						}
					}
					if (absDauPdgCode != 12 && absDauPdgCode != 14 && absDauPdgCode != 16)
					{
						mcVisible_tlv += TLorentzVector(pBHadronDaughter->getMomentum(),pBHadronDaughter->getEnergy());
					}
				}
				this->CalculateNeutrinoEnergy(mcHadron_tlv,mcVisible_tlv,mcNeutrino_tlv,mcLepton_tlv);
			}
		}
		if (m_CHadronIndex.size()!=0)
		{
			for (unsigned int n_csl = 0; n_csl < m_CHadronIndex.size(); ++n_csl)
			{
				const EVENT::MCParticle *pMCCHadron = dynamic_cast<EVENT::MCParticle*>(pLCCollection->getElementAt(m_CHadronIndex[n_csl]));
				if (!((floor(std::abs(pMCCHadron->getPDG())/100)==4) || (floor(std::abs(pMCCHadron->getPDG())/1000)==4)))
					continue;

				TLorentzVector mcHadron_tlv = TLorentzVector(pMCCHadron->getMomentum(),pMCCHadron->getEnergy());
				TLorentzVector mcVisible_tlv(0,0,0,0);
				TLorentzVector mcNeutrino_tlv(0,0,0,0);
				TLorentzVector mcLepton_tlv(0,0,0,0);
				for (unsigned int d = 0; d < (pMCCHadron->getDaughters()).size(); ++d)
				{
					const EVENT::MCParticle *pCHadronDaughter =dynamic_cast<EVENT::MCParticle*>((pMCCHadron->getDaughters())[d]);
					const int absDauPdgCode(std::abs(pCHadronDaughter->getPDG()));
					if (absDauPdgCode == 12 || absDauPdgCode == 14 || absDauPdgCode == 16)
					{
						mcNeutrino_tlv = TLorentzVector(pCHadronDaughter->getMomentum(),pCHadronDaughter->getEnergy());
						for (unsigned int l = 0; l < (pMCCHadron->getDaughters()).size(); ++l)
						{
							if (std::abs(((pMCCHadron->getDaughters())[l])->getPDG()) == absDauPdgCode - 1)
							{
								mcLepton_tlv = TLorentzVector(((pMCCHadron->getDaughters())[l])->getMomentum(),((pMCCHadron->getDaughters())[l])->getEnergy());
							}
						}
					}
					if (absDauPdgCode != 12 && absDauPdgCode != 14 && absDauPdgCode != 16)
					{
						mcVisible_tlv += TLorentzVector(pCHadronDaughter->getMomentum(),pCHadronDaughter->getEnergy());
					}
				}
				this->CalculateNeutrinoEnergy(mcHadron_tlv,mcVisible_tlv,mcNeutrino_tlv,mcLepton_tlv);
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
}

void NuCorrector::CalculateNeutrinoEnergy(TLorentzVector mcHadron_tlv,TLorentzVector mcVisible_tlv,TLorentzVector mcNeutrino_tlv,TLorentzVector mcLepton_tlv)
{
	try
	{
		TVector3 direction((mcHadron_tlv.Vect()).Unit());
		TVector3 visible_p(mcVisible_tlv.Vect());
		TVector3 mcNeutrino_p(mcNeutrino_tlv.Vect());
		TVector3 mcLepton_p(mcLepton_tlv.Vect());

		TVector3 visible_p_parallel(visible_p.Dot(direction) * direction);
		TVector3 visible_p_normal(visible_p - visible_p_parallel);
		double visible_E(mcVisible_tlv.E());
		double Hadron_mass = mcHadron_tlv.M();
		double visible_mass = mcVisible_tlv.M();

		double _A = visible_p_parallel.Mag() * (pow(Hadron_mass,2) - pow(visible_mass,2) - 2 * visible_p_normal.Mag2());
		double _B = visible_E * std::sqrt(pow(pow(Hadron_mass,2) - pow(visible_mass,2),2) - 4 * pow(Hadron_mass,2) * visible_p_normal.Mag2());
		double _C = 2 * (pow(visible_E,2) - visible_p_parallel.Mag2());
		if (pow(pow(Hadron_mass,2) - pow(visible_mass,2),2) - 4 * pow(Hadron_mass,2) * visible_p_normal.Mag2()<0)
			_B = 0;

		TVector3 recNeutrino_p_parallel_plus(((_A + _B)/_C) * direction);
		TVector3 recNeutrino_p_parallel_minus(((_A - _B)/_C) * direction);
		TVector3 recNeutrino_p_normal(-1 * visible_p_normal);
		TVector3 recNeutrino_p_plus(recNeutrino_p_parallel_plus + recNeutrino_p_normal);
		TVector3 recNeutrino_p_minus(recNeutrino_p_parallel_minus + recNeutrino_p_normal);
		m_recEnergyENuPlus.push_back(recNeutrino_p_plus.Mag());
		m_recENuPlus += recNeutrino_p_plus.Mag();
		m_recEnergyENuMinus.push_back(recNeutrino_p_minus.Mag());
		m_recENuMinus += recNeutrino_p_minus.Mag();
		if (std::abs(recNeutrino_p_plus.Mag() - mcNeutrino_p.Mag()) < std::abs(recNeutrino_p_minus.Mag() - mcNeutrino_p.Mag()))
		{
			m_recEnergyENuClose.push_back(recNeutrino_p_plus.Mag());
			m_recENuClose += recNeutrino_p_plus.Mag();
		}
		else
		{
			m_recEnergyENuClose.push_back(recNeutrino_p_minus.Mag());
			m_recENuClose += recNeutrino_p_minus.Mag();
		}
		m_mcEnergyENu.push_back(mcNeutrino_p.Mag());
		m_mcEnergyELep.push_back(mcLepton_tlv.E());
	}

	catch (...)
	{
		streamlog_out(WARNING) << "Could not reconstruct the energy of neutrino in SemiLeptonic decay" << std::endl;
	}
}
