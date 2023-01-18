// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
// O2 includes

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/CaloClusters.h"

#include "../filterTables.h"

#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct photonFilter {
  static constexpr int nTrigs{4};
  enum trigs{kPhot, kEl, kPair, kNbar};

  Produces<aod::PhotonFilters> tags;

  Configurable<float> ePhot{"ePhot", 2., "Minimal photon energy (GeV)"};
  Configurable<float> ePhot{"eEl", 2., "Minimal electron energy (GeV)"};
  Configurable<float> mPair{"mPair", 0.5, "Minimal photon pair mass (GeV)"};

  HistogramRegistry events{"events", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {
    events.add("events", "Events analysed", HistType::kTH1F, {{10, 0., 10}});
  }

  void process(o2::aod::CaloClusters const& clusters,
               o2::aod::CPVClusters const& cpvs)
  {
    // Process all clusters per TF and fill corresponding collisions
    bool keepEvent[nTrigs]{false};

    // to trace switching to new BC
    o2::InteractionRecord prevIr;
    int prevColId = -2;
    o2::InteractionRecord ir;
    for (const auto& clu : clusters) {
      if (clu.caloType() != 0)
        continue;
      if (prevColId == -2) { // first cluster
        prevColId = clu.colId();
        ir.setFromLong(clu.bc().globalBC());
      }
      if (ir != prevIr) {                      // switch to next BC
        mHistManager.fill(HIST("events"), 0.); // Total BC with PHOS clusters scanned
        // Can not fill with variable, have to fill manually
        if (keepEvent[kPhot]) {
          mHistManager.fill(HIST("events"), 1.);
          if (keepEvent[kEl]) {
            mHistManager.fill(HIST("events"), 2.);
          }
          if (keepEvent[kPair]) {
            mHistManager.fill(HIST("events"), 3.);
          }
        }
        if (keepEvent[kEl]) {
          mHistManager.fill(HIST("events"), 4.);
        }
        if (keepEvent[kPair]) {
          mHistManager.fill(HIST("events"), 5.);
        }
        if (keepEvent[kNbar]) {
          mHistManager.fill(HIST("events"), 6.);
        }

        // fill
        tags(prevColId, keepEvent[kPhot], keepEvent[kEl], keepEvent[kPair], keepEvent[kNbar]);
        // reset everething
        for (int i = 0; i < nTrigs; i++) {
          keepEvent[i] = false;
        }
        prevIr = ir;
        prevColId = clu.colId();
      }

      // Scan current cluster
      //  photons
      keepEvent[kPhot] |= (clu.e() > ePhot);
      // charged clusters above threshold
      keepEvent[kEl] |= (clu.trackdist() < 2. && clu.e() > eEl); // 2: Distance to CPV cluster in sigmas
      // antineutrons
      keepEvent[kNbar] |= (clu.ncell() > 2 && clu.m02() > 0.2 && clu.e() > 0.7 && clu.trackdist() > 2. &&
                             (clu.e() < 2. && clu.m02() > 4.5 - clu.m20()) ||
                           (clu.e() > 2. && clu.m02() > 4. - clu.m20()));

      // inv mass
      auto clu2 = clu;
      ++clu2;
      for (; clu2 != clusters.end() && clu.bc().globalBC() == clu2.bc().globalBC(); clu2++) { // scan same event only
        double m = pow(clu.e() + clu2.e(), 2) - pow(clu.px() + clu2.px(), 2) -
                   pow(clu.py() + clu2.py(), 2) - pow(clu.pz() + clu2.pz(), 2);
        if (m > mPair * mPair) {
          keepEvent[kPair] |= true;
          break;
        }
      }
    }

    // last event was not accounted yet
    mHistManager.fill(HIST("events"), 0.); // Total BC scanned
    // Can not fill with variable, have to fill manually
    if (keepEvent[kPhot]) {
      mHistManager.fill(HIST("events"), 1.);
      if (keepEvent[kEl]) {
        mHistManager.fill(HIST("events"), 2.);
      }
      if (keepEvent[kPair]) {
        mHistManager.fill(HIST("events"), 3.);
      }
    }
    if (keepEvent[kEl]) {
      mHistManager.fill(HIST("events"), 4.);
    }
    if (keepEvent[kPair]) {
      mHistManager.fill(HIST("events"), 5.);
    }
    if (keepEvent[kNbar]) {
      mHistManager.fill(HIST("events"), 6.);
    }
    // fill
    tags(prevColId, keepEvent[kPhot], keepEvent[kEl], keepEvent[kPair], keepEvent[kNbar]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<photonFilter>(cfg)};
}
