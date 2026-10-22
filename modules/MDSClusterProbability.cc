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

/** \class MDSClusterProbability
  *
  *  Computes, for every LLP of the input array (an LLPFilter output defining
  *  the detector acceptance), the probability that it gives a muon detector
  *  shower (MDS) cluster, and stores it on the LLP candidate in one of four
  *  slots (Slot: DTLoose, DTTight, CSCLoose or CSCTight). It never draws a
  *  random number and never removes a candidate: the decision is taken by
  *  MDSClusterAssignment.
  *
  *  \author Lisa Benato
  *
  */

#include "modules/MDSClusterProbability.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesDtClusterFormula.h"
#include "classes/DelphesFactory.h"

#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TString.h"

#include <memory>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

MDSClusterProbability::MDSClusterProbability() :
  fSlot(kDTLoose), fUseID(false), fEtaCutMax(999.0)
{
  fFormula = make_unique<DelphesDtClusterFormula>();
  fIDFormula = make_unique<DelphesDtClusterFormula>();
  fEtaCutFormula = make_unique<DelphesDtClusterFormula>();
}

//------------------------------------------------------------------------------

MDSClusterProbability::~MDSClusterProbability()
{
}

//------------------------------------------------------------------------------

void MDSClusterProbability::Init()
{
  TString slot = GetString("Slot", "");
  if(slot == "DTLoose")
    fSlot = kDTLoose;
  else if(slot == "DTTight")
    fSlot = kDTTight;
  else if(slot == "CSCLoose")
    fSlot = kCSCLoose;
  else if(slot == "CSCTight")
    fSlot = kCSCTight;
  else
    throw runtime_error("MDSClusterProbability: Slot must be one of DTLoose, DTTight, CSCLoose, CSCTight");

  // cluster efficiency formula in (decayR, decayZ, Ehad, Eem, multiplicity)
  fFormula->Compile(GetString("EfficiencyFormula", "0.0"));

  // optional cut-based cluster ID (CSC): N-station efficiency in (decayR, decayZ, Ehad) and eta cuts
  TString idFormula = GetString("IDFormula", "");
  fUseID = (idFormula.Length() > 0);
  if(fUseID)
  {
    fIDFormula->Compile(idFormula);
    fEtaCutFormula->Compile(GetString("EtaCutFormula", "999.0"));
    fEtaCutMax = GetDouble("EtaCutMax", 999.0);
  }

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "llpFilter/LLP"));
  fItInputArray.reset(fInputArray->MakeIterator());

  // create output array (all input candidates, unchanged apart from the probability slot)
  fOutputArray = ExportArray(GetString("OutputArray", "LLP"));
}

//------------------------------------------------------------------------------

void MDSClusterProbability::Finish()
{
}

//------------------------------------------------------------------------------

void MDSClusterProbability::Process()
{
  Candidate *candidate;
  Double_t Ehad, Eem, decayR, decayZ, multiplicity, prob, nStation, eta, cosTheta, signPz;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    const TLorentzVector &candidateDecayPosition = candidate->DecayPosition;
    const TLorentzVector &momentum = candidate->Momentum;
    decayZ = abs(candidateDecayPosition.Z());
    decayR = sqrt(pow(candidateDecayPosition.X(), 2) + pow(candidateDecayPosition.Y(), 2));
    Ehad = candidate->Ehad;
    Eem = candidate->Eem;
    multiplicity = candidate->Multiplicity;

    prob = fFormula->Eval(decayR, decayZ, Ehad, Eem, multiplicity);

    if(fUseID)
    {
      // as in CscClusterId
      cosTheta = TMath::Abs(momentum.CosTheta());
      signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
      eta = (cosTheta == 1.0 ? signPz * 999.9 : momentum.Eta());
      nStation = fIDFormula->Eval(decayR, decayZ, Ehad);
      prob *= nStation * (abs(eta) < fEtaCutMax) + (1.0 - nStation) * (abs(eta) < fEtaCutFormula->Eval(decayR, decayZ));
    }

    switch(fSlot)
    {
      case kDTLoose: candidate->ProbDTLoose = prob; break;
      case kDTTight: candidate->ProbDTTight = prob; break;
      case kCSCLoose: candidate->ProbCSCLoose = prob; break;
      case kCSCTight: candidate->ProbCSCTight = prob; break;
    }

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
