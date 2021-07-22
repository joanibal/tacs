
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

#ifndef TACS_KS_HEAT_FLUX_H
#define TACS_KS_HEAT_FLUX_H

#include "TACSFunction.h"
#include <map>
/*
  The following class implements the methods from TACSFunction.h
  necessary to calculate the KS function of temperature over the
  domain of some finite element model.

  Each class should only ever be passed to a single instance of
  TACS. If the KS function needs to be calculated for separate
  instances, this should be handled by separate instances of
  KSTemperature.

  The arguments to the KSTemperature class are:

  ksWeight:  the ks weight used in the calculation

  optional argumzents:

  alpha: scaling factor for the integration
*/
class TACSKSHeatFlux : public TACSFunction {
 public:
  static const int MAX_SURFACE_INDEX = 31;
  enum KSHeatFluxType { DISCRETE, CONTINUOUS,
                           PNORM_DISCRETE, PNORM_CONTINUOUS };
  enum KSRefType { MAX, MIN};


  TACSKSHeatFlux( TACSAssembler *_assembler, int *_elem_index,
                int *_surface_index, int _num_elems, double ksWeight,
                     double alpha=1.0  );
  ~TACSKSHeatFlux();

  /**
    Retrieve the name of the function
  */
  const char* getObjectName();


  // Set parameters for the KS function
  // ----------------------------------
  void setKSHeatFluxType( enum KSHeatFluxType type );
  void setKSRefType( enum KSRefType type );
  double getParameter();
  void setParameter( double _ksWeight );


  // Set the value of the heatflux offset for numerical stability
  // -----------------------------------------------------------
  void setRefHeatFluxOffset( TacsScalar _refHeatFlux ){
    refHeatFlux = _refHeatFlux;
  }

  /**
    Get the reference Heatflux value (either the max or min depneding on refType)
  */
  TacsScalar getReferenceHeatFlux();

  /**
     Initialize the function for the given type of evaluation
  */
  void initEvaluation( EvaluationType ftype );

  /**
     Perform an element-wise integration over this element.
  */
  void elementWiseEval( EvaluationType ftype,
                        int elemIndex, TACSElement *element,
                        double time, TacsScalar scale,
                        const TacsScalar Xpts[], const TacsScalar vars[],
                        const TacsScalar dvars[], const TacsScalar ddvars[] );

  /**
     Finalize the function evaluation for the specified eval type.
  */
  void finalEvaluation( EvaluationType ftype );

  /**
     Get the value of the function
  */
  TacsScalar getFunctionValue();

  /**
     Evaluate the derivative of the function w.r.t. state variables
  */
  void getElementSVSens( int elemIndex, TACSElement *element, double time,
                         TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                         const TacsScalar Xpts[], const TacsScalar vars[],
                         const TacsScalar dvars[], const TacsScalar ddvars[],
                         TacsScalar *elemSVSens );

  /**
     Add the derivative of the function w.r.t. the design variables
  */
  void addElementDVSens( int elemIndex, TACSElement *element,
                         double time, TacsScalar scale,
                         const TacsScalar Xpts[], const TacsScalar vars[],
                         const TacsScalar dvars[], const TacsScalar ddvars[],
                         int dvLen, TacsScalar dfdx[] );

  /**
     Evaluate the derivative of the function w.r.t. the node locations
  */
  void getElementXptSens( int elemIndex, TACSElement *element,
                          double time, TacsScalar scale,
                          const TacsScalar Xpts[], const TacsScalar vars[],
                          const TacsScalar dvars[], const TacsScalar ddvars[],
                          TacsScalar fXptSens[] );

 private:
  // The name of the function
  static const char *funcName;

  // The type of aggregation to use
  KSHeatFluxType ksType;
  KSRefType refType;

  // The weight on the ks function value
  double ksWeight;

  // The integral scaling value
  double alpha;

  // The reference HeatFlux value, the sum of exp(ksWeight*(f[i] - refHeatFlux)
  // and the value of the KS function
  TacsScalar ksHeatFluxSum, refHeatFlux;

  // Used for the case when this is used to evaluate the p-norm
  TacsScalar invPnorm;

  // List of elements and its associated surfaces
  std::map<int, int> element_to_face_key;
};

#endif // TACS_KS_HEAT_FLUX_H
