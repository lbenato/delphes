/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class MDSClusterAssignment
  *
  *  Decides which LLPs give muon detector shower (MDS) clusters, from the
  *  per-LLP probabilities stored by MDSClusterProbability: grouping of
  *  overlapping LLPs, P = 1 - prod(1 - p_i) per group, tight veto, one random
  *  number per group. See MDSClusterAssignment.h and mds_cluster_efficiency.py.
  *
  *  \author Lisa Benato
  *
  */

#include "modules/MDSClusterAssignment.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesDtClusterFormula.h"
#include "classes/DelphesFactory.h"

#include "ExRootAnalysis/ExRootConfReader.h"

#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"
#include "TVector2.h"

#include <algorithm>
#include <memory>
#include <numeric>
#include <set>
#include <vector>

using namespace std;

namespace
{
// pseudorapidity as in CscClusterId / TreeWriter (+-999.9 along the beam axis)
Double_t SafeEta(const TLorentzVector &momentum)
{
  Double_t cosTheta = TMath::Abs(momentum.CosTheta());
  Double_t signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
  return (cosTheta == 1.0 ? signPz * 999.9 : momentum.Eta());
}

Int_t FindRoot(vector<Int_t> &parent, Int_t i)
{
  while(parent[i] != i)
  {
    parent[i] = parent[parent[i]];
    i = parent[i];
  }
  return i;
}
} // namespace

//------------------------------------------------------------------------------

MDSClusterAssignment::MDSClusterAssignment() :
  fMergeDeltaR(0.6)
{
  fDTTightVetoFormula = make_unique<DelphesDtClusterFormula>();
  fCSCTightVetoFormula = make_unique<DelphesDtClusterFormula>();
}

//------------------------------------------------------------------------------

MDSClusterAssignment::~MDSClusterAssignment()
{
}

//------------------------------------------------------------------------------

void MDSClusterAssignment::Init()
{
  fMergeDeltaR = GetDouble("MergeDeltaR", 0.6);

  // veto volumes of the tight selections, formulas in (decayR, decayZ): non-zero inside
  fDTTightVetoFormula->Compile(GetString("DTTightVetoFormula", "0.0"));
  fCSCTightVetoFormula->Compile(GetString("CSCTightVetoFormula", "0.0"));

  // all LLPs (with the probabilities set by the MDSClusterProbability instances)
  fInputArray = ImportArray(GetString("InputArray", "llpFilter/LLP"));

  // LLPs that can be grouped: those in any of these arrays (DT and CSC loose acceptances)
  ExRootConfParam param = GetParam("GroupInputArray");
  fGroupInputArrays.clear();
  for(Long_t i = 0; i < param.GetSize(); ++i)
    fGroupInputArrays.push_back(ImportArray(param[i].GetString()));

  fDTLooseArray = ExportArray(GetString("DTLooseOutputArray", "DtClusterLoose"));
  fCSCLooseArray = ExportArray(GetString("CSCLooseOutputArray", "CscClusterLoose"));
  fDTTightArray = ExportArray(GetString("DTTightOutputArray", "DtClusterTight"));
  fCSCTightArray = ExportArray(GetString("CSCTightOutputArray", "CscClusterTight"));
}

//------------------------------------------------------------------------------

void MDSClusterAssignment::Finish()
{
}

//------------------------------------------------------------------------------

void MDSClusterAssignment::Process()
{
  // LLPs that can be grouped, in the order of the input array
  set<const Candidate *> groupable;
  for(const TObjArray *array : fGroupInputArrays)
    for(Int_t i = 0; i < array->GetEntriesFast(); ++i)
      groupable.insert(static_cast<const Candidate *>(array->At(i)));

  vector<Candidate *> llps;
  for(Int_t i = 0; i < fInputArray->GetEntriesFast(); ++i)
  {
    Candidate *candidate = static_cast<Candidate *>(fInputArray->At(i));
    if(groupable.count(candidate)) llps.push_back(candidate);
  }
  const Int_t n = llps.size();
  if(n == 0) return;

  // group by deltaR < fMergeDeltaR, chained (union-find)
  vector<Double_t> eta(n), phi(n);
  for(Int_t i = 0; i < n; ++i)
  {
    eta[i] = SafeEta(llps[i]->Momentum);
    phi[i] = llps[i]->Momentum.Phi();
  }
  vector<Int_t> parent(n);
  iota(parent.begin(), parent.end(), 0);
  for(Int_t i = 0; i < n; ++i)
    for(Int_t j = i + 1; j < n; ++j)
    {
      Double_t dPhi = TVector2::Phi_mpi_pi(phi[i] - phi[j]);
      if(TMath::Sqrt((eta[i] - eta[j]) * (eta[i] - eta[j]) + dPhi * dPhi) < fMergeDeltaR)
        parent[FindRoot(parent, i)] = FindRoot(parent, j);
    }

  // groups ordered by their first member
  vector<vector<Int_t>> groups;
  vector<Int_t> groupOf(n, -1);
  for(Int_t i = 0; i < n; ++i)
  {
    Int_t root = FindRoot(parent, i);
    if(groupOf[root] < 0)
    {
      groupOf[root] = groups.size();
      groups.push_back(vector<Int_t>());
    }
    groups[groupOf[root]].push_back(i);
  }

  for(size_t g = 0; g < groups.size(); ++g)
  {
    const vector<Int_t> &members = groups[g];

    // P = 1 - prod(1 - p_i) per detector and category; tight veto
    Double_t qDTLoose = 1.0, qDTTight = 1.0, qCSCLoose = 1.0, qCSCTight = 1.0;
    Bool_t vetoDT = false, vetoCSC = false;
    for(Int_t i : members)
    {
      const Candidate *c = llps[i];
      qDTLoose *= 1.0 - c->ProbDTLoose;
      qDTTight *= 1.0 - c->ProbDTTight;
      qCSCLoose *= 1.0 - c->ProbCSCLoose;
      qCSCTight *= 1.0 - c->ProbCSCTight;
      const TLorentzVector &position = c->DecayPosition;
      Double_t decayR = TMath::Hypot(position.X(), position.Y());
      Double_t decayZ = TMath::Abs(position.Z());
      if(fDTTightVetoFormula->Eval(decayR, decayZ) != 0.0) vetoDT = true;
      if(fCSCTightVetoFormula->Eval(decayR, decayZ) != 0.0) vetoCSC = true;
    }
    Double_t pDTLoose = 1.0 - qDTLoose, pCSCLoose = 1.0 - qCSCLoose;
    Double_t pDTTight = vetoDT ? 0.0 : 1.0 - qDTTight;
    Double_t pCSCTight = vetoCSC ? 0.0 : 1.0 - qCSCTight;
    Double_t sum = pDTLoose + pCSCLoose;
    if(sum > 1.0)
    {
      pDTLoose /= sum;
      pCSCLoose /= sum;
    }
    sum = pDTTight + pCSCTight;
    if(sum > 1.0)
    {
      pDTTight /= sum;
      pCSCTight /= sum;
    }

    // one random number per group
    Double_t u = gRandom->Uniform();
    for(Int_t i : members)
    {
      llps[i]->ClusterGroup = g;
      llps[i]->ClusterGroupSize = members.size();
      llps[i]->ClusterRandom = u;
    }

    // representative: the member with the largest probability in that detector and category (first one if equal)
    auto representative = [&](Float_t Candidate::*prob) {
      Candidate *best = llps[members[0]];
      for(Int_t i : members)
        if(llps[i]->*prob > best->*prob) best = llps[i];
      return best;
    };
    if(u < pDTLoose)
      fDTLooseArray->Add(representative(&Candidate::ProbDTLoose));
    else if(u > 1.0 - pCSCLoose)
      fCSCLooseArray->Add(representative(&Candidate::ProbCSCLoose));
    if(u < pDTTight)
      fDTTightArray->Add(representative(&Candidate::ProbDTTight));
    else if(u > 1.0 - pCSCTight)
      fCSCTightArray->Add(representative(&Candidate::ProbCSCTight));
  }
}

//------------------------------------------------------------------------------
