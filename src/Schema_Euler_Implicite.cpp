/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
// File:        Schema_Euler_Implicite.cpp
// Directory:   $TRUST_ROOT/src/Kernel/Schemas_Temps
// Version:     /main/57
//
//////////////////////////////////////////////////////////////////////////////

#include <Schema_Euler_Implicite.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Probleme_Couple.h>
#include <Milieu_base.h>
#include <Param.h>
#include <Debog.h>
#include <LecFicDiffuse.h>
#include <communications.h>
#include <string>
Implemente_instanciable(Schema_Euler_Implicite,"Schema_euler_implicite|Scheme_euler_implicit",Schema_Implicite_base);

//     printOn()
/////

Sortie& Schema_Euler_Implicite::printOn(Sortie& s) const
{
  return  Schema_Implicite_base::printOn(s);
}

//// readOn
//

Entree& Schema_Euler_Implicite::readOn(Entree& s)
{
  nb_ite_sans_accel_min_=-1;
  nb_ite_max=200;
  residu_old_=0;
  facsec_max_=DMAXFLOAT;
  Schema_Implicite_base::readOn(s);
  if(!le_solveur.non_nul())
    {
      Cerr << "A solver must be selected." << finl;
      Cerr << "Syntax : " << finl
           << "Solveur solver_name [ solver parameters ] " << finl;
      exit();
    }
  if (diffusion_implicite())
    {
      Cerr << "diffusion_implicite option cannot be used with an implicit time scheme." << finl;
      exit();
    }
  return s;
}

void calcul_fac_sec(double& residu_,double& residu_old,double& facsec_,const double& facsec_max_,int& nb_ite_sans_accel_,int nb_ite_sans_accel_min)
{
  double rap=1.2;
  int nb_ite_sans_accel_max=10;
  nb_ite_sans_accel_max=20000;

  if (residu_old==0)
    {
      residu_old=residu_;
      nb_ite_sans_accel_=0;
    }
  else if (((residu_old>rap*residu_)&&(nb_ite_sans_accel_>nb_ite_sans_accel_min))||(nb_ite_sans_accel_>nb_ite_sans_accel_max))
    {
      facsec_*=sqrt(rap);
      residu_old=residu_;
      facsec_=min(facsec_,facsec_max_);
      nb_ite_sans_accel_=0;
    }
  nb_ite_sans_accel_++;
  /*
    else if (residu_>residu_old*rap)
    {
    facsec_/=sqrt(rap);
    residu_old=residu_;
    facsec_=min(facsec_,facsec_max_);
    }
  */
  // Process::Journal()<<" residu" <<residu_<<" "<<residu_old<<" "<<facsec_<<finl;
  // if (residu_old!=0) facsec_=facsec_*residu_old/eps;
  //
}

void Schema_Euler_Implicite::set_param(Param& param)
{
  param.ajouter("max_iter_implicite",&nb_ite_max);
  param.ajouter("nb_ite_sans_accel_min",&nb_ite_sans_accel_min_);
  param.ajouter("facsec_max", &facsec_max_);
  Schema_Implicite_base::set_param(param);
}

bool Schema_Euler_Implicite::initTimeStep(double dt)
{
  Schema_Temps_base::initTimeStep(dt);
  residu_=0;
  return true;
}
////////////////////////////////
//                            //
// Caracteristiques du schema //
//                            //
////////////////////////////////

// Description:
//    Renvoie le nombre de valeurs temporelles a conserver.
int Schema_Euler_Implicite::nb_valeurs_temporelles() const
{
  return 3 ;
}

// Description:
//    Renvoie le nombre de valeurs temporelles futures.
//    Ici : n+1, donc 1.
int Schema_Euler_Implicite::nb_valeurs_futures() const
{
  return 1 ;
}

// Description:
//    Renvoie le le temps a la i-eme valeur future.
//    Ici : t(n+1)
double Schema_Euler_Implicite::temps_futur(int i) const
{
  assert(i==1);
  return temps_courant()+pas_de_temps();
}

// Description:
//    Renvoie le le temps le temps que doivent rendre les champs a
//    l'appel de valeurs()
//    Ici : t(n+1)
double Schema_Euler_Implicite::temps_defaut() const
{
  return temps_courant()+pas_de_temps();
}

/////////////////////////////////////////
//                                     //
// Fin des caracteristiques du schema  //
//                                     //
/////////////////////////////////////////

void Schema_Euler_Implicite::Initialiser_Champs(Probleme_base& pb)
{
  int i;
  for( i=0; i<pb.nombre_d_equations(); i++)
    {
      DoubleTab& passe = pb.equation(i).inconnue().passe();
      DoubleTab& present = pb.equation(i).inconnue().valeurs();
      if ( (pb.equation(i).que_suis_je()!="Euler") || temps_courant_==0.)
        passe = present;
    }
}

int Schema_Euler_Implicite::mettre_a_jour()
{
  calcul_fac_sec(residu_,residu_old_,facsec_,facsec_max_,nb_ite_sans_accel_,nb_ite_sans_accel_min_);
  return Schema_Temps_base::mettre_a_jour();
}

int Schema_Euler_Implicite::Iterer_Pb(Probleme_base& pb,int compteur)
{
  int ii;
  bool convergence_pb = true;
  bool convergence_eqn = false;
  int drap = (pb.milieu().que_suis_je()=="Fluide_Quasi_Compressible");
  Cout << "=======================================================================================" << finl;
  Cout << "Schema_Euler_Implicite: Implicit iteration " << compteur << " on the " << pb.que_suis_je() << " problem " << pb.le_nom() << " :" << finl;
  Cout << "=======================================================================================" << finl;
  for(int i=0; i<pb.nombre_d_equations(); i++)
    {
      //renverse l'ordre des equations si Faiblement compressible
      if ((drap)&&(i<2))
        ii = 1-i;
      else
        ii = i;
      Equation_base& eqn= pb.equation(ii);
      DoubleTab& present = eqn.inconnue().valeurs();
      DoubleTab& futur = eqn.inconnue().futur();
      double temps=temps_courant_+dt_;

      // imposer_cond_lim   sert pour la pression et pour les echanges entre pbs
      eqn.zone_Cl_dis()->imposer_cond_lim(eqn.inconnue(),temps_courant()+pas_de_temps());
      Cout<<"Solving " << eqn.que_suis_je() << " equation :" << finl;
      const DoubleTab& inut=futur;
      convergence_eqn=le_solveur.valeur().iterer_eqn(eqn, inut, present, dt_, compteur);
      convergence_pb = convergence_pb&&convergence_eqn;
      futur=present;
      eqn.zone_Cl_dis()->imposer_cond_lim(eqn.inconnue(),temps_courant()+pas_de_temps());
      present=futur;

      // La ligne suivante realise:
      // MAJ NS (donc MAJ inc)
      // MAJ modele de turbulence donc k-eps
      //   eqn.inconnue().mettre_a_jour(temps);
      //   eqn.inconnue().reculer();
      eqn.inconnue().valeur().Champ_base::changer_temps(temps);
      Cout << finl;
    }
  return (convergence_pb==true);
}

void Schema_Euler_Implicite::test_stationnaire(Probleme_base& pb)
{
  for(int i=0; i<pb.nombre_d_equations(); i++)
    {
      // traitement de la convergence en temps
      DoubleTab& passe = pb.equation(i).inconnue().passe();
      DoubleTab& present = pb.equation(i).inconnue().valeurs();
      DoubleTab& futur = pb.equation(i).inconnue().futur();
      futur = present;
      pb.equation(i).zone_Cl_dis()->imposer_cond_lim(pb.equation(i).inconnue(),temps_courant()+pas_de_temps());
      present -= passe;
      present/=dt_;
      update_critere_statio(present, pb.equation(i));
      present = passe;
    }
}

void recommencer_pb(Probleme_base& pb)
{
  for(int i=0; i<pb.nombre_d_equations(); i++)
    {
      DoubleTab& passe = pb.equation(i).inconnue().passe();
      DoubleTab& present = pb.equation(i).inconnue().valeurs();
      present=passe;
    }
}


bool Schema_Euler_Implicite::iterateTimeStep(bool& converged)
{
  Probleme_base& prob=pb_base();

  int convergence_pb = 0;
  Initialiser_Champs(prob);

  int ok=0;
  int compteur;
  while (!ok)
    {
      compteur=0;
      Cout<<" "<<finl;
      //Cout<<"Schema_Euler_Implicite : solving of problem "<<prob.que_suis_je()<<finl;
      while ((!convergence_pb)&&(compteur<nb_ite_max))
        {
          compteur++;
          prob.updateGivenFields();
          convergence_pb = Iterer_Pb(prob,compteur);
          Cout<<" "<<finl;;
        }
      if ((!convergence_pb)&&(compteur==nb_ite_max))
        {
          Cout<<"!!! Schema_Euler_Implicite has not converged at t="<< temps_courant_ << " with dt =" << dt_<< " !!!" << finl;
          double rap=sqrt(2.);
          facsec_/=rap;
          dt_/=rap;
          Cout << "!!! It restarts with dt_ =" << dt_ << " !!!" << finl;
          ok=0;
          prob.abortTimeStep();
          prob.initTimeStep(dt_);
          // PL pour corriger plantage danbs U_in_var_impl
          const int nb_eqn=prob.nombre_d_equations();
          for(int ii=0; ii<nb_eqn; ii++)
            {
              Equation_base& eqn=prob.equation(ii);
              eqn.inconnue().mettre_a_jour(temps_courant_+dt_);
            }
          // Fin
          recommencer_pb(prob);
        }
      else
        {
          Cout<<"The "<<prob.que_suis_je()<<" problem " << prob.le_nom() << " has converged after "<<compteur<<" implicit iterations."<<finl;
          ok=1;
        }

    }
  //modification : prise en compte de la possibilite de modification du dt dans Itere_Pb
  //cas du solveur pour compressible : division du pas de temps par 2
  // si la convergence n'est pas atteinte en un nombre de pas de temps donne
  //  double temps = dt_ + temps_courant_;
  test_stationnaire(prob);

  converged=true;

  return true;
}

int Schema_Euler_Implicite::faire_un_pas_de_temps_pb_couple(Probleme_Couple& pbc)
{
  // Modif B.M. : Si on fait la sauvegarde entre derivee en temps inco et mettre a jour,
  //  un calcul avec reprise n'est pas equivalent au calcul ininterrompu
  //  (front-tracking notamment). Donc je mets la sauvegarde au debut du pas de temps.
  //if (lsauv())
  //  for (int i=0;i<pbc.nb_problemes();i++)
  //    ref_cast(Probleme_base,pbc.probleme(i)).sauver();

  int convergence_pbc = 0;
  int convergence_pb;
  int i;

  for(i=0; i<pbc.nb_problemes(); i++)
    Initialiser_Champs(ref_cast(Probleme_base,pbc.probleme(i)));


  int ok=0;
  int compteur;

  while (!ok)
    {
      compteur=0;

      while ((!convergence_pbc )&&(compteur<nb_ite_max))
        {
          compteur++;
          convergence_pbc=1;


          for(i=0; i<pbc.nb_problemes()*1; i++)
            {
              Probleme_base& pb = ref_cast(Probleme_base,pbc.probleme(i));
              const int nb_eqn=pb.nombre_d_equations();
              for(int ii=0; ii<nb_eqn; ii++)
                {
                  pb.equation(ii).zone_Cl_dis()->calculer_coeffs_echange(temps_courant_+pas_de_temps());
                }

            }

          for(i=0; i<pbc.nb_problemes(); i++)
            {
              pbc.probleme(i).updateGivenFields();
              convergence_pb = Iterer_Pb(ref_cast(Probleme_base,pbc.probleme(i)),compteur);
              convergence_pbc = convergence_pbc * convergence_pb ;
            }

        }
      if ((!convergence_pbc)&&(compteur==nb_ite_max))
        {
          Cout<<"The process is restarted. " <<finl;
          facsec_/=2;
          dt_/=2;
          ok=0;
          for(int i2=0; i2<pbc.nb_problemes(); i2++)
            {
              ref_cast(Probleme_base,pbc.probleme(i2)).abortTimeStep();
              ref_cast(Probleme_base,pbc.probleme(i2)).initTimeStep(dt_);
              recommencer_pb(ref_cast(Probleme_base,pbc.probleme(i2)));
            }
        }
      else
        {
          Cout<<"The "<<pbc.que_suis_je()<<" problem " << pbc.le_nom() << " has converged after "<<compteur<<" implicit iterations."<<finl;
          ok=1;
        }
    }

  double residu_1=0;
  for(i=0; i<pbc.nb_problemes(); i++)
    {
      residu_1=residu_;
      test_stationnaire(ref_cast(Probleme_base,pbc.probleme(i)));
      residu_=max(residu_,residu_1);
    }

  return 1;
}

int Schema_Euler_Implicite::faire_un_pas_de_temps_eqn_base(Equation_base& eqn)
{
  DoubleTab& passe = eqn.inconnue().passe();
  DoubleTab& present = eqn.inconnue().valeurs();
  DoubleTab& futur = eqn.inconnue().futur();
  int compteur = 0;
  bool convergence_eqn = false;
  passe = present;
  // futur=present;
  // sert pour la pression et les couplages
  eqn.zone_Cl_dis().imposer_cond_lim(eqn.inconnue(),temps_courant()+pas_de_temps());
  //present=futur;

  compteur=0;
  Cout << finl;
  while ((!convergence_eqn)&&(compteur<nb_ite_max))
    {
      compteur++;
      Cout<<"==================================================================================" << finl;
      Cout<<"Schema_Euler_Implicite: Implicit iteration " << compteur << " on the "<<eqn.que_suis_je() << " equation of the problem "<< eqn.probleme().le_nom()<< " :" <<finl;
      Cout<<"==================================================================================" << finl;
      const DoubleTab& inut=futur;
      convergence_eqn=le_solveur.valeur().iterer_eqn(eqn, inut, present, dt_, compteur);
      futur=present;
      eqn.zone_Cl_dis().imposer_cond_lim(eqn.inconnue(),temps_courant()+pas_de_temps());
      present=futur;
    }
  present -= passe;
  present/=dt_;

  update_critere_statio(present, eqn);


  present = passe;

  return 1;
}

int Schema_Euler_Implicite::reprendre(Entree&)
{
  Nom nom_fichier(nom_du_cas());
  nom_fichier+=".dt_ev";
  LecFicDiffuse test;
  if (!test.ouvrir(nom_fichier))
    {
      Cerr << "*****************************************" << finl;
      Cerr << "***************** WARNING *************" << finl;
      Cerr << "File " << nom_fichier << " does not exist." << finl;
      Cerr << "In order to restart a calculation carried out with an implicit time scheme" << finl;
      Cerr << "it is preferable to re-read the .dt_ev file to pick up some informations" << finl;
      Cerr << "and in particular the facsec of the previous calculation." << finl;
      Cerr << "*****************************************" << finl;
      return 1;
    }

  // Reprise du facsec et du residu dans le fichier .dt_ev s'il existe
  double facsec_lu=1.;
  double residu_lu=0;
  double facsec_lu_old;

  // Test ouverture du fichier
  if (je_suis_maitre())
    {
      EFichier fichier;
      fichier.ouvrir(nom_fichier);
      Nom chaine;
      double temps=0;
      double dt;
      std::string ligne;
      // Si en tete on lit
      fichier >> chaine;
      if (chaine=="#")
        {
          // On lit la ligne complete
          std::getline(fichier.get_ifstream(), ligne);
        } // Sinon on reouvre
      else
        {
          fichier.ouvrir(nom_fichier);
        }
      // Recherche du pas de temps precedant tinit_
      fichier >> temps;
      double residu_old_old=0;
      while (!fichier.eof() && temps<tinit_)
        {

          facsec_lu_old = facsec_lu;
          fichier >> dt;
          fichier >> facsec_lu;
          fichier >> residu_lu;
          if (residu_old_==0)
            residu_old_=residu_lu;
          if (residu_old_old==0)
            residu_old_old=residu_lu;
          // On lit le reste de la ligne
          std::getline(fichier.get_ifstream(), ligne);
          fichier >> temps;

          if (facsec_lu_old != facsec_lu)
            {
              residu_old_ = residu_old_old;

              nb_ite_sans_accel_ = 1;
            }
          residu_old_old = residu_lu;
          nb_ite_sans_accel_++;
        }
    }
// GF
  calcul_fac_sec(residu_lu,residu_old_,facsec_lu,facsec_max_,nb_ite_sans_accel_,nb_ite_sans_accel_min_);
  envoyer_broadcast(facsec_lu, 0 /* pe source */);
  envoyer_broadcast(residu_lu, 0 /* pe source */);
  envoyer_broadcast(nb_ite_sans_accel_, 0 /* pe source */);
  envoyer_broadcast(residu_old_, 0 /* pe source */);

  // On prend le facsec lu uniquement s'il est entre les
  // bornes specifiees dans le jeu de donnees
  // En effet, on peut faire un calcul explicite puis une reprise en implicite
  residu_ = residu_lu;
  if (facsec_lu>=facsec_)
    {
      if (facsec_lu<=facsec_max_)
        facsec_ = facsec_lu;
      else
        facsec_ = facsec_max_;
      Cerr << "Facsec after reading in " << nom_fichier << " : " << facsec_ << finl;
      Cerr << "Residu after reading in " << nom_fichier << " : " << residu_ << finl;
    }
  else
    {
      if (Process::je_suis_maitre())
        {
          Cerr << "The readen facsec in " << nom_fichier << " : " << facsec_lu << finl;
          Cerr << "is not used since it is lower than the facsec from the data set : " << facsec_ << finl;
        }
    }
  return 1;
}
