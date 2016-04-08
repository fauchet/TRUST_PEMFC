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
// File:        Champ_implementation_P0.cpp
// Directory:   $TRUST_ROOT/src/Kernel/Champs_dis
// Version:     /main/29
//
//////////////////////////////////////////////////////////////////////////////

#include <Champ_implementation_P0.h>
#include <Domaine.h>
#include <Check_espace_virtuel.h>
#include <Zone_VF.h>
#include <Champ_base.h>

DoubleVect& Champ_implementation_P0::valeur_a_elem(const DoubleVect& position, DoubleVect& result, int poly) const
{
  const Champ_base& ch_base       = le_champ();
  int            nb_components = ch_base.nb_comp();
  const DoubleTab&  values        = ch_base.valeurs();

  assert(result.size() == nb_components);
  assert(poly>=0);
  assert(poly<values.dimension_tot(0));

  if (nb_components == 1)
    {
      if (values.nb_dim() == 1)
        {
          result(0) = values(poly);
        }
      else if (values.nb_dim() == 2)
        {
          assert(values.dimension(1) == 1);
          result(0) = values(poly,0);
        }
      else
        {
          Cerr << "Error TRUST in Champ_implementation_P0::valeur_a_elem()" << finl;
          Cerr << "The dimensions of the field are inconsistent" << finl;
          Process::exit();
        }
    }
  else
    {
      assert(values.nb_dim() == 2);
      assert(values.dimension(1) == nb_components);
      for (int i=0; i<nb_components; i++)
        {
          result(i) = values(poly,i);
        }
    }

  return result;
}

double Champ_implementation_P0::valeur_a_elem_compo(const DoubleVect& position, int poly, int ncomp) const
{
  const Champ_base& ch_base       = le_champ();
  int            nb_components = ch_base.nb_comp();
  const DoubleTab&  values        = ch_base.valeurs();

  assert(ncomp>=0);
  assert(ncomp<nb_components);
  assert(poly>=0);
  assert(poly<values.dimension_tot(0));

  double result = 0;

  if (nb_components == 1)
    {
      if (values.nb_dim() == 1)
        {
          result = values(poly);
        }
      else if (values.nb_dim() == 2)
        {
          assert(values.dimension(1) == 1);
          result = values(poly,0);
        }
      else
        {
          Cerr << "Error TRUST in Champ_implementation_P0::valeur_a_elem_compo()" << finl;
          Cerr << "The dimensions of the field are inconsistent" << finl;
          Process::exit();
        }
    }
  else
    {
      assert(values.nb_dim() == 2);
      assert(values.dimension(1) == nb_components);
      result = values(poly,ncomp);
    }

  return result;
}

DoubleTab& Champ_implementation_P0::valeur_aux_elems(const DoubleTab& positions, const IntVect& polys, DoubleTab& result) const
{
  const Champ_base& ch_base       = le_champ();
  int            nb_components = ch_base.nb_comp();
  const DoubleTab&  values        = ch_base.valeurs();

  int nb_polys = polys.size();

  assert((result.dimension_tot(0) == nb_polys)||(result.dimension(0) == nb_polys));

  if (nb_polys==0)
    return result;
  if (nb_components == 1)
    {
      assert(values.nb_dim() > 0);
      assert(result.nb_dim() > 0);
      assert(values.nb_dim() <= 2);
      assert(result.nb_dim() <= 2);

      if (values.nb_dim() == 1)
        {
          if (result.nb_dim() == 1)
            {
              for (int i=0; i<nb_polys; i++)
                {
                  int cell = polys(i);
                  assert(cell<values.dimension_tot(0));
                  // result(i) = (cell == -1) ? 0 : values(cell);
                  if (cell!=-1)
                    result(i) =values(cell);
                }
            }
          else if (result.nb_dim() == 2)
            {
              assert(result.dimension(1) == 1);
              for (int i=0; i<nb_polys; i++)
                {
                  int cell = polys(i);
                  assert(cell<values.dimension_tot(0));
                  //result(i,0) = (cell == -1) ? 0 : values(cell);
                  if (cell!=-1)
                    result(i,0) =values(cell);
                }
            }
          else
            {
              Cerr << "Error TRUST in Champ_implementation_P0::valeur_aux_elems()" << finl;
              Cerr << "The DoubleTab result has more than 2 entries" << finl;
              Process::exit();
            }
        }
      else if (values.nb_dim() == 2)
        {
          assert(values.dimension(1) == 1);
          if (result.nb_dim() == 1)
            {
              for (int i=0; i<nb_polys; i++)
                {
                  int cell = polys(i);
                  assert(cell<values.dimension_tot(0));
                  if (cell!=-1)
                    result(i) = values(cell,0);
                  //  result(i) = (cell == -1) ? 0 : values(cell,0);
                }
            }
          else if (result.nb_dim() == 2)
            {
              assert(result.dimension(1) == 1);
              for (int i=0; i<nb_polys; i++)
                {
                  int cell = polys(i);
                  assert(cell<values.dimension_tot(0));
                  if (cell!=-1)
                    result(i,0) = values(cell,0);
                  //result(i,0) = (cell == -1) ? 0 : values(cell,0);
                }
            }
          else
            {
              Cerr << "Error TRUST in Champ_implementation_P0::valeur_aux_elems()" << finl;
              Cerr << "The DoubleTab result has more than 2 entries" << finl;
              Process::exit();
            }
        }
      else
        {
          Cerr << "Error TRUST in Champ_implementation_P0::valeur_aux_elem()" << finl;
          Cerr << "The dimensions of the field are inconsistent" << finl;
          Process::exit();
        }
    }
  else
    {
      assert(values.nb_dim() == 2);
      assert(result.nb_dim() == 2);
      assert(values.dimension(1) == nb_components);
      assert(result.dimension(1) == nb_components);

      for (int i=0; i<nb_polys; i++)
        {
          int cell = polys(i);
          assert(cell<values.dimension_tot(0));
          for (int j=0; j<nb_components; j++)
            {
              if (cell!=-1)
                result(i,j) =  values(cell,j);
              //result(i,j) = (cell == -1) ? 0 : values(cell,j);
            }
        }
    }

  return result;
}

DoubleVect& Champ_implementation_P0::valeur_aux_elems_compo(const DoubleTab& positions, const IntVect& polys, DoubleVect& result, int ncomp) const
{
  const Champ_base& ch_base       = le_champ();
  int            nb_components = ch_base.nb_comp();
  const DoubleTab&  values        = ch_base.valeurs();

  int nb_polys = polys.size();

  assert(ncomp>=0);
  assert(ncomp<nb_components);
  assert(result.size() == nb_polys);

  if (nb_components == 1)
    {
      if (values.nb_dim() == 1)
        {
          for (int i=0; i<nb_polys; i++)
            {
              int cell = polys(i);
              assert(cell<values.dimension_tot(0));
              if (cell!=-1)
                result(i) =  values(cell);
              // result(i) = (cell == -1) ? 0 : values(cell);
            }
        }
      else if (values.nb_dim() == 2)
        {
          assert(values.dimension(1) == 1);
          for (int i=0; i<nb_polys; i++)
            {
              int cell = polys(i);
              assert(cell<values.dimension_tot(0));
              if (cell!=-1)
                result(i) =  values(cell,0);
              //  result(i) = (cell == -1) ? 0 : values(cell,0);
            }
        }
      else
        {
          Cerr << "Error TRUST in Champ_implementation_P0::valeur_aux_elems_compo()" << finl;
          Cerr << "The dimensions of the field are inconsistent" << finl;
          Process::exit();
        }
    }
  else
    {
      assert(values.nb_dim() == 2);
      assert(values.dimension(1) == nb_components);
      for (int i=0; i<nb_polys; i++)
        {
          int cell = polys(i);
          assert(cell<values.dimension_tot(0));
          if (cell!=-1)
            result(i) = values(cell,ncomp);
          // result(i) = (cell == -1) ? 0 : values(cell,ncomp);
        }
    }

  return result;
}

DoubleTab& Champ_implementation_P0::remplir_coord_noeuds(DoubleTab& positions) const
{
  const Zone& zone = get_zone_geom();
  positions.resize(zone.nb_elem(), Objet_U::dimension);
  zone.calculer_centres_gravite(positions);
  return positions;
}

int Champ_implementation_P0::remplir_coord_noeuds_et_polys(DoubleTab& positions, IntVect& polys) const
{
  const Zone& zone = get_zone_geom();
  int nb_elem = zone.nb_elem();
  polys.resize(nb_elem);
  remplir_coord_noeuds(positions);
  for(int i=0; i<nb_elem; i++)
    {
      polys(i)=i;
    }
  return 1;
}

DoubleTab& Champ_implementation_P0::valeur_aux_sommets_impl(DoubleTab& result) const
{
  const Champ_base& ch_base       = le_champ();
  int            nb_components = ch_base.nb_comp();
  const DoubleTab&  values        = ch_base.valeurs();

  const Zone& zone              = get_zone_geom();
  int      nb_cells          = zone.nb_elem_tot();
  int      nb_nodes          = zone.nb_som();
  int      nb_nodes_per_cell = zone.nb_som_elem();

  ArrOfInt count(nb_nodes);

  assert((result.dimension_tot(0) == nb_nodes)||(result.dimension(0) == nb_nodes));

  result = 0;
  count  = 0;

  if (nb_components == 1)
    {
      assert(values.nb_dim() > 0);
      assert(result.nb_dim() > 0);
      assert(values.nb_dim() <= 2);
      assert(result.nb_dim() <= 2);

      if (values.nb_dim() == 1)
        {
          if (result.nb_dim() == 1)
            {
              for (int i=0; i<nb_cells; i++)
                {
                  for (int j=0; j<nb_nodes_per_cell; j++)
                    {
                      int node = zone.sommet_elem(i,j);
                      if (node < nb_nodes)
                        {
                          count[node]++;
                          result(node) += values(i);
                        }
                    }
                }
              for (int i=0; i<nb_nodes; i++)
                {
                  result(i) /= count[i];
                }
            }
          else if (result.nb_dim() == 2)
            {
              assert(result.dimension(1) == 1);
              for (int i=0; i<nb_cells; i++)
                {
                  for (int j=0; j<nb_nodes_per_cell; j++)
                    {
                      int node = zone.sommet_elem(i,j);
                      if (node < nb_nodes)
                        {
                          count[node]++;
                          result(node,0) += values(i);
                        }
                    }
                }
              for (int i=0; i<nb_nodes; i++)
                {
                  result(i,0) /= count[i];
                }
            }
          else
            {
              Cerr << "Error TRUST in Champ_implementation_P0::valeur_aux_sommets_impl()" << finl;
              Cerr << "The DoubleTab result has more than 2 entries" << finl;
              Process::exit();
            }
        }
      else if (values.nb_dim() == 2)
        {
          assert(values.dimension(1) == 1);
          if (result.nb_dim() == 1)
            {
              for (int i=0; i<nb_cells; i++)
                {
                  for (int j=0; j<nb_nodes_per_cell; j++)
                    {
                      int node = zone.sommet_elem(i,j);
                      if (node < nb_nodes)
                        {
                          count[node]++;
                          result(node) += values(i,0);
                        }
                    }
                }
              for (int i=0; i<nb_nodes; i++)
                {
                  result(i) /= count[i];
                }
            }
          else if (result.nb_dim() == 2)
            {
              assert(result.dimension(1) == 1);
              for (int i=0; i<nb_cells; i++)
                {
                  for (int j=0; j<nb_nodes_per_cell; j++)
                    {
                      int node = zone.sommet_elem(i,j);
                      if (node < nb_nodes)
                        {
                          count[node]++;
                          result(node) += values(i,0);
                        }
                    }
                }
              for (int i=0; i<nb_nodes; i++)
                {
                  result(i) /= count[i];
                }
            }
          else
            {
              Cerr << "Error TRUST in Champ_implementation_P0::valeur_aux_sommets_impl()" << finl;
              Cerr << "The DoubleTab result has more than 2 entries" << finl;
              Process::exit();
            }
        }
      else
        {
          Cerr << "Error TRUST in Champ_implementation_P0::valeur_aux_sommets_impl()" << finl;
          Cerr << "The dimensions of the field are inconsistent" << finl;
          Process::exit();
        }
    }
  else
    {
      assert(result.nb_dim() == 2);
      assert(values.nb_dim() == 2);
      assert(result.dimension(1) == nb_components);
      assert(values.dimension(1) == nb_components);

      for (int i=0; i<nb_cells; i++)
        {
          for (int j=0; j<nb_nodes_per_cell; j++)
            {
              int node = zone.sommet_elem(i,j);
              if (node < nb_nodes)
                {
                  count[node]++;
                  for (int k=0; k<nb_components; k++)
                    {
                      result(node,k) += values(i,k);
                    }
                }
            }
        }
      for (int i=0; i<nb_nodes; i++)
        {
          for (int j=0; j<nb_components; j++)
            {
              result(i,j) /= count[i];
            }
        }
    }

  return result;
}

DoubleVect& Champ_implementation_P0::valeur_aux_sommets_compo_impl(DoubleVect& result, int ncomp) const
{
  const Champ_base& ch_base       = le_champ();
  int            nb_components = ch_base.nb_comp();
  const DoubleTab&  values        = ch_base.valeurs();

  const Zone& zone              = get_zone_geom();
  int      nb_cells          = zone.nb_elem_tot();
  int      nb_nodes          = zone.nb_som();
  int      nb_nodes_per_cell = zone.nb_som_elem();

  ArrOfInt count(nb_nodes);
  result=0;
  assert(ncomp>=0);
  assert(ncomp<nb_components);
  assert(result.size() == nb_nodes);

  if (nb_components == 1)
    {
      if (values.nb_dim() == 1)
        {
          for (int i=0; i<nb_cells; i++)
            {
              for (int j=0; j<nb_nodes_per_cell; j++)
                {
                  int node = zone.sommet_elem(i,j);
                  if (node < nb_nodes)
                    {
                      count[node]++;
                      result(node) += values(i);
                    }
                }
            }
          for (int i=0; i<nb_nodes; i++)
            {
              result(i) /= count[i];
            }
        }
      else if (values.nb_dim() == 2)
        {
          assert(values.dimension(1) == 1);
          for (int i=0; i<nb_cells; i++)
            {
              for (int j=0; j<nb_nodes_per_cell; j++)
                {
                  int node = zone.sommet_elem(i,j);
                  if (node < nb_nodes)
                    {
                      count[node]++;
                      result(node) += values(i,0);
                    }
                }
            }
          for (int i=0; i<nb_nodes; i++)
            {
              result(i) /= count[i];
            }
        }
      else
        {
          Cerr << "Error TRUST in Champ_implementation_P0::valeur_aux_sommets_compo_impl()" << finl;
          Cerr << "The dimensions of the field are inconsistent" << finl;
          Process::exit();
        }
    }
  else
    {
      assert(values.nb_dim() == 2);
      assert(values.dimension(1) == nb_components);
      for (int i=0; i<nb_cells; i++)
        {
          for (int j=0; j<nb_nodes_per_cell; j++)
            {
              int node = zone.sommet_elem(i,j);
              if (node < nb_nodes)
                {
                  count[node]++;
                  result(node) += values(i,ncomp);
                }
            }
        }
      for (int i=0; i<nb_nodes; i++)
        {
          result(i) /= count[i];
        }
    }

  return result;
}


int Champ_implementation_P0::imprime_P0(Sortie& os, int ncomp) const
{
  const Champ_base& cha=le_champ();
  int nb_compo_=cha.nb_comp();
  const Zone& la_zone=get_zone_geom();
  int nb_elem_tot = la_zone.nb_elem_tot();
  const DoubleTab& val = cha.valeurs();
  int elem;
  // On recalcule les centres de gravite :
  DoubleTab pos_centre;
  la_zone.calculer_centres_gravite(pos_centre);
  os << nb_elem_tot << finl;
  for (elem=0; elem<nb_elem_tot; elem++)
    {
      if (Objet_U::dimension==3)
        os << pos_centre(elem,0) << " " << pos_centre(elem,1) << " " << pos_centre(elem,2) << " " ;
      if (Objet_U::dimension==2)
        os << pos_centre(elem,0) << " " << pos_centre(elem,1) << " " ;
      if (nb_compo_ == 1)
        os << val(elem) << finl;
      else
        os << val(elem,ncomp) << finl;
    }
  os << finl;
  Cout << "Champ_P0_implementation::imprime_P0 END >>>>>>>>>> " << finl;
  return 1;
}

int Champ_implementation_P0::affecter_(const Champ_base& ch)
{
  if (le_champ().a_une_zone_dis_base() && ch.a_une_zone_dis_base() && le_champ().zone_dis_base()==ch.zone_dis_base())
    {
      // Meme support donc on utilise une methode plus rapide pour affecter_
      // que la methode generale dans Champ_Fonc_base::affecter_ ou Champ_Inc_base::affecter_
      DoubleTab noeuds;
      IntVect polys;
      remplir_coord_noeuds_et_polys(noeuds, polys);
      ch.valeur_aux_elems(noeuds, polys, le_champ().valeurs());
      return 1;
    }
  else
    {
      return 0;
    }
}
