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

#ifndef DtClusterEfficiency_h
#define DtClusterEfficiency_h

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
#include "classes/DelphesModule.h"

#include <memory>

class TIterator;
class TObjArray;
class DelphesDtClusterFormula;

class DtClusterEfficiency: public DelphesModule
{
public:
  DtClusterEfficiency();
  ~DtClusterEfficiency();

  void Init();
  void Process();
  void Finish();

private:
  std::unique_ptr<DelphesDtClusterFormula> fFormula; //!

  std::unique_ptr<TIterator> fItInputArray; //!

  const TObjArray *fInputArray = nullptr; //!

  TObjArray *fOutputArray = nullptr; //!

  ClassDef(DtClusterEfficiency, 1)
};

#endif
