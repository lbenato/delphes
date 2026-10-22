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

#ifndef MDSClusterAssignment_h
#define MDSClusterAssignment_h

/** \class MDSClusterAssignment
 *
 *  Decides which LLPs give muon detector shower (MDS) clusters, from the
 *  per-LLP probabilities stored by MDSClusterProbability (ProbDTLoose,
 *  ProbDTTight, ProbCSCLoose, ProbCSCTight):
 *
 *  - LLPs in any of the GroupInputArray arrays (the DT and CSC loose
 *    acceptances) are grouped when deltaR(LLP_i, LLP_j) < MergeDeltaR,
 *    chained; close LLPs give at most one cluster per detector;
 *  - a group gives a cluster with P = 1 - prod_i (1 - p_i), per detector and
 *    category (loose, tight);
 *  - tight: P = 0 if a member decays in the DT (CSC) veto volume,
 *    DTTightVetoFormula (CSCTightVetoFormula) in (decayR, decayZ) non-zero;
 *  - if P_DT + P_CSC > 1 in a category, both are scaled to sum to 1;
 *  - one random number u per group: a DT cluster if u < P_DT, a CSC cluster
 *    if u > 1 - P_CSC, the same u for the loose and the tight category.
 *  The cluster is represented by the member with the largest probability in
 *  that detector and category. Every grouped LLP gets ClusterGroup,
 *  ClusterGroupSize and ClusterRandom (= u).
 *
 *  Reference implementation: mds_cluster_efficiency.py.
 *
 *  \author Lisa Benato
 *
 */
#include "classes/DelphesModule.h"

#include <memory>
#include <vector>

class TIterator;
class TObjArray;
class DelphesDtClusterFormula;

class MDSClusterAssignment: public DelphesModule
{
public:
  MDSClusterAssignment();
  ~MDSClusterAssignment();

  void Init();
  void Process();
  void Finish();

private:
  Double_t fMergeDeltaR;

  std::unique_ptr<DelphesDtClusterFormula> fDTTightVetoFormula; //!
  std::unique_ptr<DelphesDtClusterFormula> fCSCTightVetoFormula; //!

  const TObjArray *fInputArray = nullptr; //!
  std::vector<const TObjArray *> fGroupInputArrays; //!

  TObjArray *fDTLooseArray = nullptr; //!
  TObjArray *fCSCLooseArray = nullptr; //!
  TObjArray *fDTTightArray = nullptr; //!
  TObjArray *fCSCTightArray = nullptr; //!

  ClassDef(MDSClusterAssignment, 1)
};

#endif
