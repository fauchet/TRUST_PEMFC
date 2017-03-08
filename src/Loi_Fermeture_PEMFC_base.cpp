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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Loi_Fermeture_PEMFC_base.cpp
// Directory : $PEMFC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Loi_Fermeture_PEMFC_base.h>
#include <Param.h>
#include <Discretisation_base.h>
#include <Champ_P1NC.h>
#include <Milieu_base.h>
#include <Probleme_base.h>
#include <Zone_VF.h>
#include <Debog.h>
#include <Check_espace_virtuel.h>


Implemente_base( Loi_Fermeture_PEMFC_base, "Loi_Fermeture_PEMFC_base", Loi_Fermeture_base ) ;

Sortie& Loi_Fermeture_PEMFC_base::printOn( Sortie& os ) const
{
  Loi_Fermeture_base::printOn( os );
  return os;
}

Entree& Loi_Fermeture_PEMFC_base::readOn( Entree& is )
{
  is_cas4_=0;
  Loi_Fermeture_base::readOn( is );
  return is;
}




void Loi_Fermeture_PEMFC_base::discretiser(const Discretisation_base& dis)

{
  Loi_Fermeture_base::discretiser(dis);
}




#include <Linear_algebra_tools_impl.h>


double  eval_matrice_diffusion(const double& cO2,const double& cvap,const double& cN2, DoubleTab& Diff,int elem)
{
  // Valeur
  double  T=353.15;
  double  R=8.314;
  double  Rp=1e-5;
  double  epsilon=0.7;
  double  MO2=32e-3;
  double  MN2=28e-3;
  double  Mvap=18e-3;
  double  mu=1.2e-5;
  double  tau=3;
  double  K=6e-12;
  double  Krg=1;
  double pi=3.1415926;
  //
  double  cg = cO2 + cN2 + cvap;
  double  Xvap =cvap/cg ;
  double  XO2=cO2/cg;
  double  XN2=cN2/cg;
  double  Pg=cg*R*T;
  double  DO2N2 = 6.43e-5*pow( T ,1.823)/Pg;
  double  DO2vap = 4.26e-6*pow( T ,2.334)/Pg;
  double  DN2vap = 4.45e-6*pow( T ,2.334)/Pg;

  double DN2O2=DO2N2;
  double DvapO2=DO2vap;
  double DvapN2=DN2vap;


  double  DKO2=2./3.*Rp*sqrt(8*R*T/(pi*MO2));
  double  DKvap=2./3.*Rp*sqrt(8*R*T/(pi*Mvap));
  double  DKN2=2./3.*Rp*sqrt(8*R*T/(pi*MN2));
  double  invDAO2vap=1./DO2vap+1./DKO2;
  double  invDAO2N2=1./DO2N2+1./DKO2;
  double  invDAvapO2=1./DvapO2+1./DKvap;
  double  invDAvapN2=1./DvapN2+1./DKvap;
  double  invDAN2O2=1./DN2O2+1./DKN2;
  double  invDAN2vap=1./DN2vap+1./DKN2;
  double  AK =0.75/Rp*sqrt(pi*R*T*0.5);
  double  ss=(XO2*sqrt(MO2)+Xvap*sqrt(Mvap)+XN2*sqrt(MN2));
  double  Ac=epsilon*mu/(cg*tau*tau*K*Krg*ss);
  double  AA=1./(1./Ac+1./AK);


  // Matrice
  double a11=-XN2*invDAO2N2-Xvap*invDAO2vap;
  double a12=XO2*invDAvapO2;
  double a13=XO2*invDAN2O2;

  double a21=Xvap*invDAO2vap;
  double a22=-XO2*invDAvapO2-XN2*invDAvapN2;
  double a23=Xvap*invDAN2vap;

  double a31= -AA*sqrt(MO2);

  double a32= -AA*sqrt(Mvap);
  double a33= -AA*sqrt(MN2);



  Matrice33 M(a11,a12,a13,
              a21,a22,a23,
              a31,a32,a33);
  Matrice33 invM;
  Matrice33::inverse(M,invM);
  // Cerr<< invM(0,0)<<finl;
#if 1
  Matrice33 coef(1-XO2,-XO2,-XO2,
                 -Xvap,1-Xvap,-Xvap,
                 R*T,R*T,R*T);
#else
  Matrice33 coef(1,0,0,
                 0,1-Xvap,0,
                 R*T,R*T,R*T);
#endif
  /* cg=1;
  Diff(elem,0*3+0)=cg*invM(0,0)+R*T*invM(0,2);
  Diff(elem,0*3+1)=cg*invM(0,1)+R*T*invM(0,2);
  Diff(elem,0*3+2)=R*T*invM(0,2);
  Diff(elem,1*3+0)=cg*invM(1,0)+R*T*invM(1,2);
  Diff(elem,1*3+1)=cg*invM(1,1)+R*T*invM(1,2);
  Diff(elem,1*3+2)=R*T*invM(1,2);
  Diff(elem,2*3+0)=cg*invM(2,0)+R*T*invM(2,2);
  Diff(elem,2*3+1)=cg*invM(2,1)+R*T*invM(2,2);
  Diff(elem,2*3+2)=R*T*invM(2,2);
  */
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      {
        double p=0;
        for (int k=0; k<3; k++)
          p+=invM(i,k)*coef(k,j);
        p*=epsilon/(tau*tau);
        Diff(elem,i*3+j)=-p;
      }

  // double epsilon=0.7;
  //double tau=3.;

  return 1.;
}

double  eval_matrice_diffusion_cas4(const double& XO2,const double& Xvap,const double& XN2b, DoubleTab& Diff,int elem)
{
  // Valeur
  double  T=353.15;
  double  R=8.314;
  //
  //double  Pg=1.5e5;
  double  DO2N2 = 6.43e-5*pow( T ,1.823);
  double  DO2vap = 4.26e-6*pow( T ,2.334);
  double  DN2vap = 4.45e-6*pow( T ,2.334);

  double DN2O2=DO2N2;
  double DvapO2=DO2vap;
  double DvapN2=DN2vap;

  double XN2=1.-Xvap-XO2;
//  XN2=XN2b;

  double  invDAO2vap=1./DO2vap;
  double  invDAO2N2=1./DO2N2;
  double  invDAvapO2=1./DvapO2;
  double  invDAvapN2=1./DvapN2;
  double  invDAN2O2=1./DN2O2;
  double  invDAN2vap=1./DN2vap;

  // Matrice
  double a11=-XN2*invDAO2N2-Xvap*invDAO2vap;
  double a12=XO2*invDAvapO2;
  double a13=XO2*invDAN2O2;

  double a21=Xvap*invDAO2vap;
  double a22=-XO2*invDAvapO2-XN2*invDAvapN2;
  double a23=Xvap*invDAN2vap;

  Matrice33 invM;
#if 1
  double a31= 1;

  double a32= 1;
  double a33=1;



  Matrice33 M(a11,a12,a13,
              a21,a22,a23,
              a31,a32,a33);
  Matrice33::inverse(M,invM);
  Matrice33 coef(1./(R*T),0,0,
                 0,1./(R*T),0,
                 0,0,0);

  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      {
        double p=0;
        for (int k=0; k<3; k++)
          p+=invM(i,k)*coef(k,j);
        Diff(elem,i*3+j)=-p;
      }
#else


  double a31=XN2*invDAO2N2;
  double a32=XN2*invDAvapN2;
  double a33=-XO2*invDAO2N2-Xvap*invDAvapN2;

  Matrice33 M(a11,a12,a13,
              a21,a22,a23,
              a31,a32,a33);
  for (int j=0; j<3; j++)
    {
      double t=0;
      for (int i=0; i<3; i++)
        t+=M(i,j);
      if (!est_egal(t,0)) throw;
    }
  for (int i=2; i<3; i++)
    for (int j=0; j<3; j++)
      {
        M(i,j)+=1.e3;
      }
  Matrice33::inverse(M,invM);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      {
        double p=0;
        p+=invM(i,j);
        Diff(elem,i*3+j)=-p;
      }
#endif

  for (int i=0; i<3*0; i++)
    for (int j=0; j<3; j++)
      {
        Cout <<"Diff "<<i<<" "<<j<<" "<< Diff(elem,i*3+j) <<" "<<M (i,j)<<finl;
      }
  //Process::exit();
  assert(Diff(elem,0)>=0);

  return 1.;
}


void  bidouille_nu(DoubleTab& nu,const DoubleTab&   inconnue_org,const Zone_VF& zone_VF,int is_cas4)
{

  int nb_comp=inconnue_org.dimension(1);
  assert(nb_comp==3);
  ArrOfDouble c(3);


  const IntTab& elem_faces=zone_VF.elem_faces();
  int nbf=elem_faces.dimension(1);

  double invnbf=1./nbf;
  int nb_elem=zone_VF.zone().nb_elem();
  assert(nb_elem==nu.dimension(0));
  assert(zone_VF.nb_faces()==inconnue_org.dimension(0));

  assert(elem_faces.dimension(0)==nb_elem);
  for (int elem=0; elem<nb_elem; elem++)
    {
      c=0;
      for (int nc=0; nc<nb_comp; nc++)
        for (int f=0; f<nbf; f++)
          {
            int face=elem_faces(elem,f);
            c(nc)+=inconnue_org(face,nc);
          }
      c*=invnbf;
      if (is_cas4)
        eval_matrice_diffusion_cas4(c[0],c[1],c[2],  nu, elem);
      else
        eval_matrice_diffusion(c[0],c[1],c[2],  nu, elem);
    }
  nu.echange_espace_virtuel();
  double maxc=local_max_vect(nu);
  if (local_min_vect(nu)<0)
    {
      Cerr<<"oooooo !!!!!!!!!!!! "<< local_min_vect(nu) <<finl;
      //abort();
    }
  Cerr<<"iii "<<nu(0,0)<<finl;
  Cerr<<"max nu "<<maxc<<finl;

}

// copy of discretiation_tools
void Discretisation_tools_cells_to_faces(const Champ_base& He,  Champ_base& Hf)
{
  DoubleTab& tabHf=Hf.valeurs();
  const DoubleTab& tabHe=He.valeurs();
  Debog::verifier("element_face entreee",tabHe);
  assert_espace_virtuel_vect(tabHe);
  const Zone_dis_base& zone_dis_base=He.zone_dis_base();

  const Zone_VF& zone_vf= ref_cast(Zone_VF,zone_dis_base);
  // en realite on fait P1B vers face
  //assert(tabHe.dimension_tot(0)==zone_dis_base.nb_elem_tot());
  assert(tabHf.dimension_tot(0)==zone_vf.nb_faces_tot());

  const IntTab& elem_faces=zone_vf.elem_faces();
  const DoubleVect& volumes=zone_vf.volumes();
  const DoubleVect& volumes_entrelaces=zone_vf.volumes_entrelaces();



  tabHf=0;
  int nb_face_elem=elem_faces.dimension(1);
  int nb_elem_tot=zone_dis_base.nb_elem_tot();
  int nb_dim=tabHe.nb_dim();


  double coeffb=nb_face_elem;
  double coeffi=coeffb;
  if (zone_vf.que_suis_je()=="Zone_VDF")
    {
      coeffb=1;
      coeffi=2;
    }

  if (nb_dim==1)
    {
      for (int ele=0; ele<nb_elem_tot; ele++)
        {
          for (int s=0; s<nb_face_elem; s++)
            {

              tabHf(elem_faces(ele,s))+=tabHe(ele)*volumes(ele);
              //	      Cerr<<elem_faces(ele,s)<<" "<<tabHe(ele)*volumes(ele)<<" "<<tabHf(elem_faces(ele,s))<<finl;
            }
        }
      for (int f=0; f<zone_vf.premiere_face_int(); f++)
        tabHf(f)/=volumes_entrelaces(f)*coeffb;
      for (int f=zone_vf.premiere_face_int(); f<zone_vf.nb_faces(); f++)
        tabHf(f)/=volumes_entrelaces(f)*coeffi;

    }
  else
    {

      if (tabHf.nb_dim()==1)
        {
          assert(coeffi==2);
          assert(coeffb==1);
          // VDF
          //abort();
          for (int ele=0; ele<nb_elem_tot; ele++)
            {
              for (int s=0; s<nb_face_elem; s++)
                {
                  int face=elem_faces(ele,s);
                  //for (int comp=0;comp<nb_comp;comp++)
                  int comp=zone_vf.orientation()[face];
                  tabHf(face)+=tabHe(ele,comp)*volumes(ele);
                }
            }
          for (int f=0; f<zone_vf.premiere_face_int(); f++)
            tabHf(f)/=volumes_entrelaces(f)*coeffb;
          for (int f=zone_vf.premiere_face_int(); f<zone_vf.nb_faces(); f++)
            tabHf(f)/=volumes_entrelaces(f)*coeffi;
        }
      else
        {
          //abort();
          int nb_comp=tabHf.dimension(1);
          for (int ele=0; ele<nb_elem_tot; ele++)
            {
              for (int s=0; s<nb_face_elem; s++)
                {
                  int face=elem_faces(ele,s);
                  for (int comp=0; comp<nb_comp; comp++)
                    tabHf(face,comp)+=tabHe(ele,comp)*volumes(ele);
                }
            }

          for (int f=0; f<zone_vf.premiere_face_int(); f++)
            for (int comp=0; comp<nb_comp; comp++)
              tabHf(f,comp)/=volumes_entrelaces(f)*coeffb;
          for (int f=zone_vf.premiere_face_int(); f<zone_vf.nb_faces(); f++)
            for (int comp=0; comp<nb_comp; comp++)
              tabHf(f,comp)/=volumes_entrelaces(f)*coeffi;

        }
    }
  tabHf.echange_espace_virtuel();
//  Cerr<<min_array(tabHe)<<" elem  "<<max_array(tabHe)<<finl;
//  Cerr<<min_array(tabHf)<<" face  "<<max_array(tabHf)<<finl;

  Debog::verifier("element_face sortie",tabHf);
}



void Loi_Fermeture_PEMFC_base::calculer_Ni_ud(const double& temps)
{
  const  DoubleTab& inco= equation().inconnue().valeurs();
  int nb_comp=inco.dimension(1);
  // calcul du gradient
  DoubleTab grad(0, nb_comp, Objet_U::dimension);
  equation().zone_dis().zone().creer_tableau_elements(grad);

  DoubleTab& Ni=Ni_.valeurs();
  DoubleTab& ude=ud_.valeurs();


  const Champ_base& ch_rho_face=mon_probleme().get_champ("masse_volumique");
  const DoubleTab& rho_face=ch_rho_face.valeurs();
  assert(rho_face.get_md_vector()==Um_.valeurs().get_md_vector());

  const Champ_base& vitesse=mon_probleme().get_champ("vitesse");

  ArrOfDouble mi(3);
  mi[0]=32e-3;
  mi[1]=18e-3;
  mi[2]=28e-3;
  const DoubleTab& nu=diffu_.valeurs();
  Champ_P1NC::calcul_gradient(inco,grad,ref_cast(Zone_Cl_VEF,equation().zone_Cl_dis().valeur()));
  for (int elem=0; elem<grad.dimension(0); elem++)
    {

      for (int j=0; j<dimension; j++)
        {
          ude(elem,j)=0;

          for (int i=0; i<nb_comp; i++)

            {
              double p=0;
              for (int k=0; k<nb_comp; k++)
                p+=nu(elem,i*nb_comp+k)*grad(elem,k,j);

              Ni(elem,i*dimension+j)=-p;

              ude(elem,j)-=p*mi[i] ;
            }
          // ude(elem,j)/=rho_elem(elem);
        }
      for (int j=0; j<dimension; j++)
        {
          double pt=0;
          for (int i=0; i<nb_comp; i++)
            pt+=Ni(elem,i*dimension+j);
          if (!est_egal(pt,0))
            {
              Cerr<<"Ni oooo "<< pt<<" "<<Ni(elem,0*dimension+j)<<" "<<Ni(elem,1*dimension+j)<<" "<<Ni(elem,2*dimension+j)<< finl;

              //   exit();
            }

        }
    }
  ude.echange_espace_virtuel();
  // on devrait appeler Discretsisation_tools
  Discretisation_tools_cells_to_faces(ud_,Um_);
  Um_.valeurs()*=-1;
  tab_divide_any_shape(Um_.valeurs(),rho_face);
  Um_.valeurs()+=vitesse.valeurs();
  if (1)
    Um_.valeurs()=vitesse.valeurs();

  Ni_.changer_temps(temps);
  ud_.changer_temps(temps);
  Um_.changer_temps(temps);
}
void Loi_Fermeture_PEMFC_base::mettre_a_jour(double temps)
{
  const Zone_VF& zone_VF = ref_cast(Zone_VF,equation().zone_dis().valeur());
  const DoubleTab& inconnue_org=equation().inconnue().valeurs();
  //  int nb_comp=inconnue_org.dimension(1);
  // diffusivity_read_from_datafile_.valeur().dimensionner(1,nb_comp*nb_comp);
  DoubleTab& nu=diffu_.valeurs();
  diffu_->changer_temps(temps);

  bidouille_nu(nu,inconnue_org,zone_VF,is_cas4_);
  if (is_cas4_)
    calculer_Ni_ud(temps);
}



