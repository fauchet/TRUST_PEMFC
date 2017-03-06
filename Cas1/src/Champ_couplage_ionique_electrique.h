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
// File      : Champ_couplage_ionique_electrique.h
// Directory : pemfc pemfc/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Champ_couplage_ionique_electrique_included
#define Champ_couplage_ionique_electrique_included

#include <Champ_Fonc_P0_base.h>
#include <Ref_Probleme_base.h>


/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Champ_couplage_ionique_electrique
//
// <Description of class Champ_couplage_ionique_electrique>
//
/////////////////////////////////////////////////////////////////////////////

class Champ_couplage_ionique_electrique : public Champ_Fonc_P0_base
{

  Declare_instanciable( Champ_couplage_ionique_electrique ) ;

public :
  virtual void mettre_a_jour(double);
  virtual void associer_zone_dis_base(const Zone_dis_base&) ;
  virtual const Zone_dis_base& zone_dis_base() const ;
  const DoubleTab& get_derivee() const { return derivee_; } ;
protected :
  REF(Probleme_base) pb_ionique_,pb_electrique_;
  int sens_;
  REF(Zone_dis_base) zdb_;
  DoubleTab derivee_;
};

inline void Champ_couplage_ionique_electrique::associer_zone_dis_base(const Zone_dis_base& la_zone_dis_base)
{
  zdb_=la_zone_dis_base;
}

inline const Zone_dis_base& Champ_couplage_ionique_electrique::zone_dis_base() const
{
  return zdb_.valeur();
}

#endif /* Champ_couplage_ionique_electrique_included */
