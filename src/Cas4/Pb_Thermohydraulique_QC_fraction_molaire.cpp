/****************************************************************************
* Copyright (c) 2015, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : Pb_Thermohydraulique_QC_fraction_molaire.cpp
// Directory : $PEMFC_ROOT/src/Cas4
//
/////////////////////////////////////////////////////////////////////////////

#include <Pb_Thermohydraulique_QC_fraction_molaire.h>

Implemente_instanciable( Pb_Thermohydraulique_QC_fraction_molaire, "Pb_Thermohydraulique_QC_fraction_molaire", Pb_Thermohydraulique_QC ) ;

Sortie& Pb_Thermohydraulique_QC_fraction_molaire::printOn( Sortie& os ) const
{
  Pb_Thermohydraulique_QC::printOn( os );
  return os;
}

Entree& Pb_Thermohydraulique_QC_fraction_molaire::readOn( Entree& is )
{
  Pb_Thermohydraulique_QC::readOn( is );
  return is;
}

int Pb_Thermohydraulique_QC_fraction_molaire::nombre_d_equations() const
{
  return 3;
}
const Equation_base& Pb_Thermohydraulique_QC_fraction_molaire::equation(int i) const
{
  assert ((i>=0) && (i<=2));
  if (i == 0)
    return eq_hydraulique;
  else
    {
      if (i==1)
        return eq_thermique;
      else
        return eq_fraction_molaire_;
    }
}
Equation_base& Pb_Thermohydraulique_QC_fraction_molaire::equation(int i)
{
  assert ((i>=0) && (i<=2));
  if (i == 0)
    return eq_hydraulique;
  else
    {
      if (i==1)
        return eq_thermique;
      else
        return eq_fraction_molaire_;
    }
}


int Pb_Thermohydraulique_QC_fraction_molaire::verifier()
{
  return  Pb_Thermohydraulique_QC::verifier();
}
