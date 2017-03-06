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
// File      : Terme_Puissance_Thermique_implicite_VDF_Elem.cpp
// Directory : pemfc pemfc/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Terme_Puissance_Thermique_implicite_VDF_Elem.h>
#include <Champ_couplage_ionique_electrique.h>
#include <Matrice_Morse.h>

Implemente_instanciable( Terme_Puissance_Thermique_implicite_VDF_Elem, "Puissance_Thermique_implicite_VDF_P0_VDF", Terme_Puissance_Thermique_VDF_Elem ) ;

Sortie& Terme_Puissance_Thermique_implicite_VDF_Elem::printOn( Sortie& os ) const
{
  Terme_Puissance_Thermique_VDF_Elem::printOn( os );
  return os;
}

Entree& Terme_Puissance_Thermique_implicite_VDF_Elem::readOn( Entree& is )
{
  Terme_Puissance_Thermique_VDF_Elem::readOn( is );
  return is;
}
void Terme_Puissance_Thermique_implicite_VDF_Elem::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& matrice) const
{
  const DoubleVect& vol = ref_cast(Zone_VF,equation().zone_dis().valeur()).volumes();
  const  Champ_couplage_ionique_electrique& P= ref_cast(Champ_couplage_ionique_electrique,la_puissance.valeur());
  const DoubleTab& der=P.get_derivee();
  int taille=der.dimension(0);	
  for (int i=0;i<taille;i++)
	{
	  matrice.coef(i,i)+=der(i)*vol(i);
	  assert(der(i)>=0);
	}
}
