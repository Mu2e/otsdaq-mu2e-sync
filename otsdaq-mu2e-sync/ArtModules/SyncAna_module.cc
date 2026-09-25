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

// ROOT
#include <TBufferFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TH1F.h>

namespace mu2e {
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

    void bookHistograms();
    void fillHistograms();

  private:
    //--------------------------------
    // Histogram info
    //--------------------------------
    struct ChannelHist_t { // per-channel histograms
      TH1* nhits;
      TH1* tot;
      TH1* delta_t; // time since last hit started
      TH1* delta_tend; // time since last hit ended
    };
    struct GlobalHist_t { // global histograms
      TH1* nhits;
      TH1* tot;
      TH1* delta_t; // between different channels
      TH2* t_v_t_evt;
      TH2* t_v_t;
    };

    //--------------------------------
    // Channel info
    //--------------------------------
    struct ChannelInfo_t {
      int hit_count = 0;
      int evt_since_last_hit = -1; // information about the last hit in the channel
      double prev_hit_time = 0.;
      double prev_hit_tot = 0.;
      double prev_time(bool end = false) {
        constexpr double time_evt = 1.e5; // 100 us FIXME
        if(evt_since_last_hit < 0) return 0.;
        double time = prev_hit_time - evt_since_last_hit*time_evt;
        if(end) time += prev_hit_tot;
        return time;
      }
    };

    //--------------------------------
    // Inputs
    //--------------------------------
    std::string msdTag_;
    int diagLevel_, evtCounter_;

    //--------------------------------
    // Data
    //--------------------------------
    const MSDHitCollection* msdHits_;
    ChannelHist_t channelHists_[2];
    GlobalHist_t globalHists_;
    ChannelInfo_t channelInfo_[2];

  };
} // namespace mu2e

mu2e::SyncAna::SyncAna(Parameters const& conf) :
  art::EDAnalyzer(conf), msdTag_(conf().msdTag()), diagLevel_(conf().diag()), evtCounter_(0) {
}

void mu2e::SyncAna::beginJob() {
  bookHistograms();
}

void mu2e::SyncAna::bookHistograms() {
  // Create the histograms of interest
  art::ServiceHandle<art::TFileService> tfs;
  globalHists_.nhits = tfs->make<TH1D>("nhits", "N(MSD hits)", 50, 0, 50);
  globalHists_.tot = tfs->make<TH1D>("tot", "Hit TOT", 400, 0, 2000);
  globalHists_.delta_t = tfs->make<TH1F>("delta_t", "#Deltat (ns)", 300, -50, 250);
  globalHists_.t_v_t_evt = tfs->make<TH2F>("t_v_t_evt", "Event hit time;Channel 0 time (ns); Channel 1 time (ns)", 500, 0, 100000., 500, 0, 100000.);
  globalHists_.t_v_t = tfs->make<TH2F>("t_v_t", "Global hit time;Channel 0 time (ns); Channel 1 time (ns)", 100, 0, 1.e5, 100, 0, 1.e5);
  for(int channel = 0; channel < 2; ++channel) {
    channelHists_[channel].nhits      = tfs->make<TH1D>(std::format("channel_{}_nhits", channel).c_str(), "N(MSD hits)", 50, 0, 50);
    channelHists_[channel].tot        = tfs->make<TH1D>(std::format("channel_{}_tot", channel).c_str(), "Hit TOT", 400, 0, 2000);
    channelHists_[channel].delta_t    = tfs->make<TH1F>(std::format("channel_{}_delta_t", channel).c_str(), "#Deltat (ns)", 1000, 0., 2.e6);
    channelHists_[channel].delta_tend = tfs->make<TH1F>(std::format("channel_{}_delta_tend", channel).c_str(), "#Deltat_{end} (ns)", 1000, 0., 2.e3);
  }
}

void mu2e::SyncAna::fillHistograms() {
  if(!msdHits_) return;
  globalHists_.nhits->Fill(msdHits_->size());


  double prev_time[2] = {channelInfo_[0].prev_time(), channelInfo_[1].prev_time()};
  double prev_tot [2] = {channelInfo_[0].prev_hit_tot, channelInfo_[1].prev_hit_tot};
  for(const auto& hit : *msdHits_) {
    const double time = hit.time();
    const int    ID   = hit.channelID();
    const double ptime = prev_time[ID];
    const double ptend = ptime + prev_tot[ID];
    if(ptime > 0.) {
      channelHists_[ID].delta_t->Fill(time - ptime);
      channelHists_[ID].delta_tend->Fill(time - ptend);
    }
    prev_time[ID] = time;
    prev_tot [ID] = hit.tot();
    globalHists_.tot->Fill(hit.tot());
    channelHists_[ID].tot->Fill(hit.tot());
  }
  for(int channel = 0; channel < 2; ++channel) {
    channelHists_[channel].nhits->Fill(channelInfo_[channel].hit_count);
  }

  // Fill 2D plots
  const size_t nhits = msdHits_->size();
  for(size_t i = 0; i < nhits; ++i) {
    for(size_t j = i+1; j < nhits; ++j) {
      const auto& hit_1 = msdHits_->at(i);
      const auto& hit_2 = msdHits_->at(j);
      if(hit_1.channelID() == hit_2.channelID()) continue;
      const double t_1 = (hit_1.channelID() == 0) ? hit_1.time() : hit_2.time();
      const double t_2 = (hit_1.channelID() == 0) ? hit_2.time() : hit_1.time();
      globalHists_.t_v_t_evt->Fill(t_1, t_2);
    }
  }

  // look for a coincidence
  for(size_t i = 0; i < msdHits_->size(); ++i) {
    const auto& hit_i = msdHits_->at(i);
    for(size_t j = i + 1; j < msdHits_->size(); ++j) {
      const auto& hit_j = msdHits_->at(j);
      if(hit_i.channelID() == hit_j.channelID()) continue;
      const double dt = hit_j.time() - hit_i.time();
      if(std::fabs(dt) < 1000.) globalHists_.delta_t->Fill(dt);
    }
  }
}

void mu2e::SyncAna::analyze(art::Event const& event) {
  ++evtCounter_;

  art::Handle<mu2e::MSDHitCollection> msdH;
  if(!event.getByLabel(msdTag_, msdH) || !msdH.product()) {
    msdHits_ = nullptr;
  } else {
    msdHits_ = msdH.product();
  }

  // Count hits by channel
  if(msdHits_) {
    for(const auto& hit : *msdHits_) {
      const int ID = hit.channelID();
      ++channelInfo_[ID].hit_count;
    }
  }

  // fill the histograms of interest
  if(msdHits_) fillHistograms();

  // Update per-channel last hit info at the end of the event
  for(int channel = 0; channel < 2; ++channel) {
    auto& info = channelInfo_[channel];
    if(info.hit_count > 0) { // set the most recent hit time to the last hit
      for(const auto& hit : *msdHits_) {
        const double time = hit.time();
        const double ID   = hit.channelID();
        if(ID != channel) continue;
        info.prev_hit_time = time;
        info.prev_hit_tot = hit.tot();
      }
      info.evt_since_last_hit = 1;
    } else { // increment events since last hit
      ++info.evt_since_last_hit;
    }
    // Reset the counter
    info.hit_count = 0;
  }
}

void mu2e::SyncAna::endJob() {
}

void mu2e::SyncAna::beginRun(const art::Run&) {}

DEFINE_ART_MODULE(mu2e::SyncAna)
