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

#ifndef MDSClusterProbability_h
#define MDSClusterProbability_h

/** \class MDSClusterProbability
 *
 *  Computes, for every LLP of the input array (an LLPFilter output defining
 *  the detector acceptance), the probability that it gives a muon detector
 *  shower (MDS) cluster, and stores it on the LLP candidate in one of four
 *  slots (Slot: DTLoose, DTTight, CSCLoose or CSCTight). It never draws a
 *  random number and never removes a candidate: the decision is taken by
 *  MDSClusterAssignment, which combines the probabilities of overlapping LLPs
 *  and of both detectors with one random number per group.
 *
 *  The probability is EfficiencyFormula(decayR, decayZ, Ehad, Eem,
 *  multiplicity). If IDFormula is set (cut-based CSC cluster ID, as in
 *  CscClusterId), it is multiplied by
 *    P_ID = n * (|eta| < EtaCutMax) + (1 - n) * (|eta| < EtaCutFormula),
 *  with n = IDFormula(decayR, decayZ, Ehad) the efficiency of the N-station
 *  requirement and eta the LLP pseudorapidity.
 *
 *  \author Lisa Benato
 *
 */
#include "classes/DelphesModule.h"

#include <memory>

class TIterator;
class TObjArray;
class DelphesDtClusterFormula;

class MDSClusterProbability: public DelphesModule
{
public:
  MDSClusterProbability();
  ~MDSClusterProbability();

  void Init();
  void Process();
  void Finish();

private:
  enum ESlot { kDTLoose, kDTTight, kCSCLoose, kCSCTight };
  ESlot fSlot;

  std::unique_ptr<DelphesDtClusterFormula> fFormula; //!
  std::unique_ptr<DelphesDtClusterFormula> fIDFormula; //!
  std::unique_ptr<DelphesDtClusterFormula> fEtaCutFormula; //!
  bool fUseID;
  Double_t fEtaCutMax;

  std::unique_ptr<TIterator> fItInputArray; //!

  const TObjArray *fInputArray = nullptr; //!

  TObjArray *fOutputArray = nullptr; //!

  ClassDef(MDSClusterProbability, 1)
};

#endif
