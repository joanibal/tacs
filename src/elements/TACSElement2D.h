/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_ELEMENT_2D_H
#define TACS_ELEMENT_2D_H

#include "TACSElement.h"
#include "TACSElementModel.h"
#include "TACSElementBasis.h"

class TACSElement2D : public TACSElement {
 public:
  static const int MAX_VARS_PER_NODE = 8;

  TACSElement2D( TACSElementModel *_model, TACSElementBasis *_basis );
  ~TACSElement2D();

  // Get the layout properties of the element
  int getVarsPerNode();
  int getNumNodes();
  ElementLayout getLayoutType();

  /**
    Retrieve the global design variable numbers associated with this element
  */
  int getDesignVarNums( int elemIndex, int dvLen, int dvNums[] );

  /**
    Set the element design variables from the design vector
  */
  void setDesignVars( int elemIndex, int dvLen, const TacsScalar dvs[] );

  /**
    Get the element design variables values
  */
  void getDesignVars( int elemIndex, int dvLen, TacsScalar dvs[] );

  /**
    Get the lower and upper bounds for the design variable values
  */
  void getDesignVarRange( int elemIndex, int dvLen,
                          TacsScalar lb[], TacsScalar ub[] );

  /**
    Add the residual to the provided vector
  */
  void addResidual( int elemIndex, double time, const TacsScalar *Xpts,
                    const TacsScalar *vars, const TacsScalar *dvars,
                    const TacsScalar *ddvars, TacsScalar *res );

  /**
    Add the residual and Jacobians to the arrays
  */
  void addJacobian( int elemIndex, double time,
                    double alpha, double beta, double gamma,
                    const TacsScalar *Xpts, const TacsScalar *vars,
                    const TacsScalar *dvars, const TacsScalar *ddvars,
                    TacsScalar *res, TacsScalar *mat );

  /**
    Add the derivative of the product of the adjoint variables w.r.t.
    the material design variables
  */
  void addAdjResProduct( int elemIndex, double time, double scale,
                        const TacsScalar psi[], const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[],
                        int dvLen, TacsScalar dvSens[] );

  /**
    Add the derivative of the product of the adjoint variables and the
    residuals with respect to the node locations
  */
  void addAdjResXptProduct( int elemIndex, double time, double scale,
                            const TacsScalar psi[], const TacsScalar Xpts[],
                            const TacsScalar vars[], const TacsScalar dvars[],
                            const TacsScalar ddvars[], TacsScalar fXptSens[] );

  /**
    Compute the output data for visualization
  */
  void getOutputData( int elemIndex, ElementType etype, int write_flag,
                      const TacsScalar Xpts[], const TacsScalar vars[],
                      const TacsScalar dvars[], const TacsScalar ddvars[],
                      int ld_data, TacsScalar *data );

 private:
  TACSElementModel *model;
  TACSElementBasis *basis;
};

#endif // TACS_ELEMENT_2D_H
