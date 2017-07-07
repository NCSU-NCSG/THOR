/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

/*
 *  this user object writes the mesh to file
 */

#include "DumpMesh.h"

// MOOSE includes
#include "MooseMesh.h"

// System headers
#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>

template<>
InputParameters validParams<DumpMesh>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addParam<std::string >("filename", "mesh_out.dat", "Filename of mesh output file.");
  return params;
}

DumpMesh::DumpMesh(const InputParameters & params) :
    ElementUserObject(params)
{
}

void
DumpMesh::initialize()
{
}

void
DumpMesh::threadJoin(const UserObject & uo)
{
  /*
  if (_correct_XS)
  {
    const DumpMesh & pps = dynamic_cast<const DumpMesh &>(uo);

    for (std::map<SubdomainID, Real>::iterator it=_mesh_volumes.begin();
         it != _mesh_volumes.end(); it++)
    {
      std::map<SubdomainID, Real>::const_iterator itt = pps._mesh_volumes.find(it->first);
      if (itt != pps._mesh_volumes.end())
        it->second += itt->second;
    }
  }
  */
}

void
DumpMesh::execute()
{
}

void
DumpMesh::finalize()
{
  // Write the mesh here.
  if (processor_id() == 0)
  {
    libMesh::MeshBase* si_mesh = & _mesh.getMesh();
    std::ofstream fid;
    fid.open(getParam<std::string>("filename").c_str(),std::ios::out);

    // Print the header
    fid << si_mesh->n_nodes() << std::endl;
    fid << si_mesh->n_elem() << std::endl;
    fid << std::endl << std::endl;

    // nodes
    MeshBase::const_node_iterator nodes_it    = si_mesh->active_nodes_begin();
    MeshBase::const_node_iterator nodes_end   = si_mesh->active_nodes_end();

    unsigned int l = 0;
    std::map<dof_id_type, unsigned int> id2node;
    for (; nodes_it != nodes_end; ++nodes_it)
    {
      l++;
      Node * node = *nodes_it;
      dof_id_type node_id = node->id();
      id2node[node_id] = l;
      fid << l << " " << std::scientific << std::setw(14) << std::setprecision(6) << (*node)(0)
               << " " << std::scientific << std::setw(14) << std::setprecision(6) << (*node)(1)
               << " " << std::scientific << std::setw(14) << std::setprecision(6) << (*node)(2) << std::endl;
    }

    // elements
    MeshBase::const_element_iterator el_it = si_mesh->active_elements_begin();
    MeshBase::const_element_iterator el_end = si_mesh->active_elements_end();

    unsigned int e = 0;
    std::set<subdomain_id_type> all_subdomain_ids;
    for (; el_it != el_end; ++el_it)
    {
      e++;
      Elem * elem = *el_it;
      fid << e << " " <<  elem->subdomain_id() << " " << 0 << std::endl;
      all_subdomain_ids.insert(elem->subdomain_id());
    }

    std::cout << "--These subdomains are in the mesh--" << std::endl;
    std::set<subdomain_id_type>::iterator it;
    for (it = all_subdomain_ids.begin(); it != all_subdomain_ids.end(); ++it)
      std::cout << "Subdomain ID " << *it << std::endl;

    std::vector<SubdomainName> subdomain_names;
    subdomain_names.push_back("M-1-TRI");
    subdomain_names.push_back("M-2-TRI");
    subdomain_names.push_back("M-3-TRI");
    subdomain_names.push_back("M-4-TRI");
    subdomain_names.push_back("M-5-TRI");
    subdomain_names.push_back("M-6-TRI");
    subdomain_names.push_back("M-7-TRI");
    for (unsigned int j = 0; j < subdomain_names.size(); j++)
        std::cout << "Subdomain Name " << subdomain_names[j] << " => (ID) " << _mesh.getSubdomainID(subdomain_names[j]) << std::endl;

    el_it = si_mesh->active_elements_begin();
    e = 0;
    for (; el_it != el_end; ++el_it)
    {
      e++;
      Elem * elem = *el_it;
      fid << e << " ";
      for (unsigned int n = 0; n < elem->n_nodes(); n++)
      {
        Node * node = elem->get_node(n);
        dof_id_type node_id = node->id();
        fid << std::setw(10)<< id2node[node_id] << " ";
      }
      fid << std::endl;
    }

    fid.close();
  }
}
