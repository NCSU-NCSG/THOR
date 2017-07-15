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

#ifndef DUMPMESH_H
#define DUMPMESH_H

#include "ElementUserObject.h"

//Forward Declarations
class DumpMesh;

template<>
InputParameters validParams<DumpMesh>();

class DumpMesh : public ElementUserObject
{
public:
  DumpMesh(const InputParameters & params);

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & uo);
  virtual void finalize();

};
#endif /* DumpMesh_H */
