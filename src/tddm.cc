/*
 *  Copyright (C) 2023 Joachim Vandekerckhove <joachim@uci.edu>
 *                     Adriana F. Chavez De la Pena <achavezd@uci.edu>
 *
 *  When using this module, please cite:
 *      Wabersich, D., & Vandekerckhove, J. (2014). Extending JAGS: A
 *      tutorial on adding custom distributions to JAGS (with a
 *      diffusion model example). Behavior Research Methods, 46, 15-28.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 *
 */
#include <module/Module.h>
#include <distributions/DTDDM.h>

using std::vector;

namespace jags {
namespace tddm {

class TDDMModule : public Module {
  public:
    TDDMModule();
    ~TDDMModule();
};

TDDMModule::TDDMModule() : Module("tddm")
{
  DTDDM *tddmdist;
  tddmdist = new DTDDM();
  //load distributions
  insert(tddmdist);
}

TDDMModule::~TDDMModule()
{
  vector<Function*> const &fvec = functions();
  for (unsigned int i = 0; i < fvec.size(); ++i) {
    delete fvec[i];
  }
  vector<Distribution*> const &dvec = distributions();
  for (unsigned int i = 0; i < dvec.size(); ++i) {
    delete dvec[i];
  }
}

} // namespace tddm
} // namespace jags

jags::tddm::TDDMModule _tddm_module;
