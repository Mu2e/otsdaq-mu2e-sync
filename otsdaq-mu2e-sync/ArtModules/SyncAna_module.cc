// Author: M. Stortini
// Analyze time sync data

// Framework
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

// Offline
#include "Offline/RecoDataProducts/inc/MSDHit.hh"

// Online
#include "otsdaq/Macros/CoutMacros.h"
#include "otsdaq/Macros/ProcessorPluginMacros.h"
#include "otsdaq/MessageFacility/MessageFacility.h"
#include "otsdaq/NetworkUtilities/TCPSendClient.h"

// ROOT
#include <TBufferFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>

// TRACE
#include "trace.h"
#define TRACE_NAME "SyncAna"

using namespace mu2e;
namespace ots {
class SyncAna : public art::EDAnalyzer {
public:
	struct Config {
		// clang-format off
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<std::string> msdTag{Name("msdTag")   , Comment("Timing paddle hit collection tag")};
      fhicl::Atom<int>         diag  {Name("diagLevel"), Comment("Diagnostic level"), 0};
		// clang-format on
	};

	typedef art::EDAnalyzer::Table<Config> Parameters;

	explicit SyncAna(Parameters const& conf);

	void analyze(art::Event const& event) override;
	void beginRun(art::Run const&) override;
	void beginJob() override;
	void endJob() override;

	void bookHistograms(const int index);
	void fillHistograms(const int index);

private:
	// Histogram info
	enum { kMaxHists = 100 };
	struct Hist_t {
		TH1* nhits;
		TH1* delta_t;
		TH1* delta_t_min;
		TH1* delta_t_prev;
	};

	// Inputs
	std::string msdTag_;
	int diagLevel_, evtCounter_;
	art::ServiceHandle<art::TFileService> tfs_;

	// Data
	const MSDHitCollection* msdHits_;
	Hist_t* hist_[kMaxHists];  // histogram books

	double prev_hit_time_ = 0.;  // for comparing separation between hits (cosmics)
	int evt_since_last_hit_ = -1;
};
}  // namespace ots

ots::SyncAna::SyncAna(Parameters const& conf)
	: art::EDAnalyzer(conf)
	, msdTag_(conf().msdTag())
	, diagLevel_(conf().diag())
	, evtCounter_(0) {
	__COUT__ << "[SyncAna::" << __func__ << "] Constructor" << std::endl;
}

void ots::SyncAna::beginJob() {
	__COUT__ << "[SyncAna::" << __func__ << "] Beginning job" << std::endl;
	for (int index = 0; index < kMaxHists; ++index) hist_[index] = nullptr;

	bookHistograms(0);
}

void ots::SyncAna::bookHistograms(const int index) {
	if (index >= kMaxHists)
		throw cet::exception("BADCONFIG")
			<< "Histogram index " << index << " larger than maximum " << kMaxHists;
	if (hist_[index])
		throw cet::exception("BADCONFIG")
			<< "Histogram index " << index << " already booked!";
	hist_[index] = new Hist_t;

	// Create the histograms of interest
	art::TFileDirectory dir = tfs_->mkdir(std::format("hist_{}", index));
	Hist_t* Hist = hist_[index];
	Hist->nhits = dir.make<TH1D>("nhits", "N(MSD hits)", 50, 0, 50);
	Hist->delta_t = dir.make<TH1F>("delta_t", "#Deltat (ns)", 300, -50, 250);
	Hist->delta_t_min = dir.make<TH1F>("delta_t_min", "min #Deltat (ns)", 100, -50, 50);
	Hist->delta_t_prev = dir.make<TH1F>("delta_t_prev", "#Deltat (us)", 10000, 0., 5.e6);
}

void ots::SyncAna::fillHistograms(const int index) {
	if (!hist_[index])
		throw cet::exception("BADCONFIG")
			<< "Histogram index " << index << " not booked!";
	Hist_t* Hist = hist_[index];
	Hist->nhits->Fill(msdHits_->size());

	constexpr double time_evt = 1.e5;
	double deltaT = 0.0, time_prev(prev_hit_time_ - evt_since_last_hit_ * time_evt);
	for (size_t i = 0; i < msdHits_->size(); i++) {
		const double time = msdHits_->at(i).time();
		if (i > 0) {
			deltaT = time - time_prev;  // - 10000.;
			Hist->delta_t->Fill(deltaT);
		}
		else if (evt_since_last_hit_ >= 0) {
			Hist->delta_t_prev->Fill((time - time_prev) / 1000.);  // fill in us
		}
		time_prev = time;

		// look for the minimum delta T in the event
		double min_deltaT(-1.), time_min(0.);
		for (size_t j = i + 1; j < msdHits_->size(); j++) {
			const double time_j = msdHits_->at(j).time();
			const double dt = std::abs(time - time_j);
			if (min_deltaT < 0. || dt < min_deltaT) {
				min_deltaT = dt;
				time_min = time_j;
			}
		}
		if (min_deltaT >= 0. && min_deltaT < 1000.)
			Hist->delta_t_min->Fill(time - time_min);
	}
}

void ots::SyncAna::analyze(art::Event const& event) {
	++evtCounter_;
	TLOG(TLVL_DEBUG + 20) << "[SyncAna::" << __func__ << "] Analyzing event "
						  << evtCounter_ << ": " << event.id() << std::endl;

	art::Handle<mu2e::MSDHitCollection> msdH;
	if (!event.getByLabel(msdTag_, msdH) || !msdH.product()) {
		TLOG(TLVL_WARNING) << "[SyncAna::" << __func__
						   << "] No MSD hit collection found!\n";
		msdHits_ = nullptr;
	}
	else {
		msdHits_ = msdH.product();
		TLOG(TLVL_DEBUG + 21) << "[SyncAna::" << __func__ << "] Event " << event.id()
							  << " has " << msdHits_->size() << " MSD hits\n";
	}

	// fill the histograms of interest
	if (msdHits_) fillHistograms(0);

	if (!msdHits_->empty()) {  // new hit found
		prev_hit_time_ = msdHits_->back().time();
		evt_since_last_hit_ = 1;
	}
	else if (evt_since_last_hit_ >= 0)
		++evt_since_last_hit_;  // increment if no hit was found
}

void ots::SyncAna::endJob() {
	__COUT__ << "[SyncAna::" << __func__ << "] Ending job, saw " << evtCounter_
			 << " events\n";
}

void ots::SyncAna::beginRun(const art::Run&) {}

DEFINE_ART_MODULE(ots::SyncAna)
