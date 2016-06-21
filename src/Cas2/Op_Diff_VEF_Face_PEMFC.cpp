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
// File      : Op_Diff_VEF_Face_PEMFC.cpp
// Directory : $PEMFC_ROOT/src/Cas2
//
/////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_VEF_Face_PEMFC.h>
#include <Param.h>
#include <EChaine.h>

Implemente_instanciable( Op_Diff_VEF_Face_PEMFC, "Op_Diff_VEFPEMFC_const_P1NC", Op_Diff_VEF_Face_Matricial ) ;

Sortie& Op_Diff_VEF_Face_PEMFC::printOn( Sortie& os ) const
{
  Op_Diff_VEF_Face_Matricial::printOn( os );
  return os;
}

Entree& Op_Diff_VEF_Face_PEMFC::readOn( Entree& is )
{
  //Op_Diff_VEF_Face_Matricial::readOn( is );
  //diffusivity_read_from_datafile_.fixer_nb_comp(9);
  //diffusivity_read_from_datafile_.typer("Champ_uniforme");
  EChaine es(" Champ_uniforme 9 0 0 0 0 0 0 0 0 0 ");
  es >> diffusivity_read_from_datafile_;

  Op_Diff_VEF_Face_Matricial::readOn( is );
  //Param param(que_suis_je());

  //param.lire_avec_accolades_depuis(is);

  return is;
}
void Op_Diff_VEF_Face_PEMFC::set_param(Param& param)
{
  // abort();
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
        Diff(elem,i*3+j)=p;
      }


  return 1.;
}


void  bidouille_nu(DoubleTab& nu,const DoubleTab&   inconnue_org,const Zone_VEF zone_VEF)
{

  int nb_comp=inconnue_org.dimension(1);
  assert(nb_comp==3);
  ArrOfDouble c(3);


  const IntTab& elem_faces=zone_VEF.elem_faces();
  int nbf=elem_faces.dimension(1);

  double invnbf=1./nbf;
  int nb_elem=zone_VEF.zone().nb_elem();
  assert(nb_elem==nu.dimension(0));
  assert(zone_VEF.nb_faces()==inconnue_org.dimension(0));

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

      eval_matrice_diffusion(c[0],c[1],c[2],  nu, elem);
    }
  nu.echange_espace_virtuel();
  double epsilon=0.7;
  double tau=3.;
  nu*=epsilon/tau/tau;
  nu*=-1;
  double maxc=local_max_vect(nu);
  if (local_min_vect(nu)<0)
    {
      Cerr<<"oooooo !!!!!!!!!!!! "<< local_min_vect(nu) <<finl;
      //abort();
    }
  Cerr<<"iii "<<nu(0,0)<<finl;
  Cerr<<"max nu "<<maxc<<finl;

}




void Op_Diff_VEF_Face_PEMFC::mettre_a_jour(double)
{
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
  const DoubleTab& inconnue_org=equation().inconnue().valeurs();
  //  int nb_comp=inconnue_org.dimension(1);
  // diffusivity_read_from_datafile_.valeur().dimensionner(1,nb_comp*nb_comp);
  remplir_nu(nu_);

  bidouille_nu(nu_,inconnue_org,zone_VEF);
}






