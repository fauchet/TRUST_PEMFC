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
// File      : Champ_couplage_ionique_electrique.cpp
// Directory : pemfc pemfc/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Champ_couplage_ionique_electrique.h>
#include <Param.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Interprete.h>
#include <Domaine.h>
#include <Champ_Generique_base.h>
#include <Zone_VDF.h>

Implemente_instanciable( Champ_couplage_ionique_electrique, "Champ_couplage_ionique_electrique", Champ_Fonc_P0_base ) ;
// XD Champ_couplage_ionique_electrique field_base Champ_couplage_ionique_electrique 1 not_set

Sortie& Champ_couplage_ionique_electrique::printOn( Sortie& os ) const
{
  Champ_Fonc_P0_base::printOn( os );
  return os;
}

Entree& Champ_couplage_ionique_electrique::readOn( Entree& is )
{
  // Champ_Fonc_P0_base::readOn( is );
  Param param(que_suis_je());
  Nom nom_pb_ionique,nom_pb_electrique,nom_pb;
  param.ajouter("nom_pb_ionique",&nom_pb_ionique,Param::REQUIRED);  // XD_ADD_P chaine not_set
  param.ajouter("nom_pb_electrique",&nom_pb_electrique,Param::REQUIRED);   // XD_ADD_P chaine not_set
  param.ajouter("nom_mon_pb",&nom_pb,Param::REQUIRED);   // XD_ADD_P chaine not_set

  param.lire_avec_accolades(is);
  pb_ionique_=ref_cast(Probleme_base,interprete().objet(nom_pb_ionique));
  pb_electrique_=ref_cast(Probleme_base,interprete().objet(nom_pb_electrique));
  const Probleme_base& mon_pb=ref_cast(Probleme_base,interprete().objet(nom_pb));
  
  valeurs_.resize(0);
  mon_pb.domaine().zone(0).creer_tableau_elements(valeurs_);
  associer_zone_dis_base(mon_pb.domaine_dis().zone_dis(0));
  //mettre_a_jour(0.);
  fixer_nb_comp(1);
  corriger_unite_nom_compo();

  
  sens_=1;
  if (nom_pb!=nom_pb_electrique)
    sens_=-1;
  return is;
}


double eval_f(const double& psi,const double& phi,const double& t)
{
  double nFsurRT=2*96500/(8.314*353.15);
  double io=1e-5;
  double alpha=0.5;
  double Erev=1.18;
  double eta=(psi-phi)-Erev;
  
  double res=0;
  double x1=alpha*nFsurRT*eta;
  if (x1>50)
    x1=50;
  res+=exp(x1);
  double x2=-(1-alpha)*nFsurRT*eta;
  if (x2>50)
    x2=50;
  res-=exp(x2);
  //  double t2=2.;
  res*=2.5e7*io; // *0.1*(1.+9.*((t/t2)*(t<t2)+(t>=t2)));
  //Cerr<<t<<" jjjjjjj "<<(1.+9.*((t/100.)*(t<100.)+(t>=100.)))<<finl; //abort();
  return res;
  
}





double eval_der_f(const double& psi,const double& phi,const double& t)
{
  double nFsurRT=2*96500/(8.314*353.15);
  double io=1e-5;
  double alpha=0.5;
  double Erev=1.18;
  double eta=(psi-phi)-Erev;
  
  double res=0;
  double x1=alpha*nFsurRT*eta;
  if (x1>50)
    x1=50;
  res+=exp(x1)*alpha*nFsurRT;
  double x2=-(1-alpha)*nFsurRT*eta;
  if (x2>50)
    x2=50;
  res+=exp(x2)*(1-alpha)*nFsurRT;
  //  double t2=2.;
  res*=2.5e7*io; // *0.1*(1.+9.*((t/t2)*(t<t2)+(t>=t2)));
  //Cerr<<t<<" jjjjjjj "<<(1.+9.*((t/100.)*(t<100.)+(t>=100.)))<<finl; //abort();
  return res;
  
}
void Champ_couplage_ionique_electrique::mettre_a_jour( double un_temps )
{
 
  changer_temps(un_temps);
  const Champ_base& psi=pb_electrique_.valeur().equation(0).inconnue();
  const Champ_base& phi=pb_ionique_.valeur().equation(0).inconnue();

  DoubleTab phi_inter,psi_inter;
  //  int is_vef=(phi.que_suis_je().find("P0")==-1);
  if (sens_==1)
    {
//#define old 
      // je suis su cote electrique
#ifdef old
      valeurs()=-5e37;
      affecter_(phi);
      phi_inter=valeurs();
#else
      Champ sto;
      phi_inter=pb_electrique_.valeur().get_champ_post("phi_inter").get_champ(sto).valeurs();
#endif  
      if (psi.que_suis_je().find("P0")==-1)
	 psi_inter=pb_electrique_.valeur().get_champ_post("psi_elem").get_champ(sto).valeurs();
      else
	psi_inter.ref(psi.valeurs());
      
    }
  else
    {
      // je suis su cote ionique
#ifdef old
      valeurs()=-5e37;
      affecter_(psi);
      psi_inter=valeurs();
#else
      Champ sto;
      psi_inter=pb_ionique_.valeur().get_champ_post("psi_inter").get_champ(sto).valeurs();
#endif
      if (phi.que_suis_je().find("P0")==-1)
	phi_inter=pb_ionique_.valeur().get_champ_post("phi_elem").get_champ(sto).valeurs();
      else
	phi_inter.ref(phi.valeurs());
      
    }
  DoubleTab& val=valeurs();
  if (!derivee_.get_md_vector().non_nul())
    derivee_=val;
  derivee_=0;
  for (int i=0;i<val.dimension(0);i++)
    {
      if ((psi_inter(i)>-1e30)&&(phi_inter(i)>-1e30))
	{
	  val(i)=eval_f(psi_inter(i),phi_inter(i),un_temps);
	 // derivee_(i)=eval_der_f(psi_inter(i),phi_inter(i),un_temps);
	}
      else
	val(i)=0.;
    }
  if (sens_==1)
    val*=-1;
  else 
    assert(sens_==-1);

  Cerr<<"iiiiiiiiiiiii "<<mp_min_vect(val)<<" "<<mp_max_vect(val)<<finl;
}
