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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Loi_Etat_Melange_GP_Fraction_Molaire.h
// Directory:   $TRUST_ROOT/src/ThHyd/Quasi_Compressible
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Loi_Etat_Melange_GP_Fraction_Molaire_included
#define Loi_Etat_Melange_GP_Fraction_Molaire_included

#include <Loi_Etat_Melange_GP.h>
#include <Champ_Inc_base.h>
#include <List.h>
#include <Ref_Champ_Inc_base.h>
#include <DoubleTab.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//     classe Loi_Etat_Melange_GP_Fraction_Molaire
//     Cette classe represente la loi d'etat pour un melange de gaz parfaits.
//
// .SECTION voir aussi
//     Fluide_Quasi_Compressible Loi_Etat_base Loi_Etat_GP
//////////////////////////////////////////////////////////////////////////////


class Loi_Etat_Melange_GP_Fraction_Molaire : public Loi_Etat_Melange_GP
{
  Declare_instanciable_sans_constructeur(Loi_Etat_Melange_GP_Fraction_Molaire);

public :

  Loi_Etat_Melange_GP_Fraction_Molaire();
  void associer_fluide(const Fluide_Quasi_Compressible&);
  virtual  void associer_inconnue(const Champ_Inc_base& inconnue);
  virtual  void associer_cg(const Champ_Inc_base& inconnue);
  void calculer_Cp();
  virtual void calculer_tab_Cp(DoubleTab& cp) const;
  void calculer_lambda();
  void calculer_alpha();
  void calculer_mu();
  void calculer_mu_sur_Sc();
  void calculer_mu0();

  void calculer_masse_volumique();
  double calculer_masse_volumique(double,double) const;
  virtual double calculer_masse_volumique_case(double P,double T,double r, int som) const;

  virtual void initialiser_inco_ch();
  void calculer_masse_molaire();
  void calculer_masse_molaire(DoubleTab& M) const;
  virtual void rabot(int futur=0);
  inline const DoubleTab& masse_molaire() const;
  inline DoubleTab& masse_molaire();
  const ArrOfDouble& get_masses_molaires() const
  {
    return Mi_ ;
  };
protected :
  int rabot_,verif_fraction_,nb_famille_grains_;
  double Sc_;
  //  DoubleTab Masse_mol_mel;
  ArrOfDouble Mi_;
  REF(Champ_Inc_base) ref_X_,ref_cg_;
};


inline const DoubleTab& Loi_Etat_Melange_GP_Fraction_Molaire::masse_molaire() const
{
  return Masse_mol_mel;
}

inline DoubleTab& Loi_Etat_Melange_GP_Fraction_Molaire::masse_molaire()
{
  return Masse_mol_mel;
}

#endif
