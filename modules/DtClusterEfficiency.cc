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

/** \class DtClusterEfficiency
  *
  *  This module parameterizes the reconstruction efficiency of LLP decay
  *  clusters in the CMS barrel muon drift-tube (DT) chambers, analogous to
  *  CscClusterEfficiency for the endcap CSC system. In addition to the decay
  *  position and the hadronic/EM energy of the LLP, the efficiency formula
  *  can also depend on the decay-product multiplicity, i.e. the number of
  *  LLP daughters that were used to compute Ehad and Eem (set upstream by
  *  LLPFilter).
  *
  *  \author Lisa Benato
  *
  */

#include "modules/DtClusterEfficiency.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesDtClusterFormula.h"
#include "classes/DelphesFactory.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

DtClusterEfficiency::DtClusterEfficiency()
{
  fFormula = make_unique<DelphesDtClusterFormula>();
}

//------------------------------------------------------------------------------

DtClusterEfficiency::~DtClusterEfficiency()
{
}

//------------------------------------------------------------------------------

void DtClusterEfficiency::Init()
{
  // read DtClusterEfficiency formula
  fFormula->Compile(GetString("EfficiencyFormula", "1.0"));

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray.reset(fInputArray->MakeIterator());

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void DtClusterEfficiency::Finish()
{
}

//------------------------------------------------------------------------------

void DtClusterEfficiency::Process()
{
  Candidate *candidate;
  Double_t Ehad, Eem, decayR, decayZ, multiplicity;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    const TLorentzVector &candidateDecayPosition = candidate->DecayPosition;
    decayZ = abs(candidateDecayPosition.Z());
    decayR = sqrt(pow(candidateDecayPosition.X(), 2) + pow(candidateDecayPosition.Y(), 2));
    Ehad = candidate->Ehad;
    Eem = candidate->Eem;
    multiplicity = candidate->Multiplicity;
    // apply an efficency formula
    if(gRandom->Uniform() > fFormula->Eval(decayR, decayZ, Ehad, Eem, multiplicity)) continue;

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
