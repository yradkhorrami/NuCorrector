#include "NuCorrector.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#endif // MARLIN_USE_AIDA

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

using namespace lcio ;
using namespace marlin ;


NuCorrector aNuCorrector ;

NuCorrector::NuCorrector() :

        Processor("NuCorrector"),
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

//	std::cout << "NuCorrector::end()  " << name()
//	<< " processed " << _nEvt << " events in " << _nRun << " runs "
//	<< std::endl ;

}

void NuCorrector::Clear()
{

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

