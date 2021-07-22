/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2019 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TACSKSHeatFlux.h"
#include "TACSAssembler.h"
#include "TACSElementAlgebra.h"

TACSKSHeatFlux::TACSKSHeatFlux( TACSAssembler *_assembler,
                            int *elem_index,
                            int *face_index,
                            int num_elems, 
                            double _ksWeight, 
                            double _alpha ):
 TACSFunction(_assembler,TACSFunction::ENTIRE_DOMAIN,
             TACSFunction::TWO_STAGE, 0){

  ksWeight = _ksWeight;
  alpha = _alpha;
  ksType = CONTINUOUS;
  refType = MAX;
  
  // Initialize the reference heatFlux value and KS sum to default values
  // that will be overwritten later.
  refHeatFlux = -1e20;
  ksHeatFluxSum = 0.0;
  invPnorm = 0.0;

  // Set the element indices
  for ( int i = 0; i < num_elems; i++ ){ 
    int value = 0;
    if (face_index[i] >= 0 && face_index[i] < MAX_SURFACE_INDEX){
      value = 1 << face_index[i];
    }

    std::map<int, int>::iterator it = element_to_face_key.find(elem_index[i]);

    // Perform a bitwise or operation if the element has already
    // been added to the domain
    if (it != element_to_face_key.end()){
      value = it->second | value;
    }

    // Set the value in the map
    element_to_face_key[elem_index[i]] = value;
  }

  // Set the domain of the function
  setDomain(num_elems, elem_index);
}

TACSKSHeatFlux::~TACSKSHeatFlux(){

}

/*
  KSHeatFlux function name
*/
const char* TACSKSHeatFlux::funcName = "TACSKSHeatFlux";


/*
  Set the KS aggregation type
*/
void TACSKSHeatFlux::setKSHeatFluxType( enum KSHeatFluxType type ){
  ksType = type;
}

void TACSKSHeatFlux::setKSRefType( enum KSRefType type ){
  refType = type;
}

/*
  Retrieve the KS aggregation weight
*/
double TACSKSHeatFlux::getParameter(){
  return ksWeight;
}

/*
  Set the KS aggregation parameter
*/
void TACSKSHeatFlux::setParameter( double _ksWeight ){
  ksWeight = _ksWeight;
}
/*
  Return the function name
*/
const char *TACSKSHeatFlux::getObjectName(){
  return funcName;
}
/*
  Retrieve the function value
*/
TacsScalar TACSKSHeatFlux::getFunctionValue(){
  // Compute the final value of the KS function on all processors
    TacsScalar ksHeatFlux;
    if (refType == MAX){
      ksHeatFlux = refHeatFlux + log(ksHeatFluxSum/alpha)/ksWeight;
    }
    else if (refType == MIN){
      ksHeatFlux = refHeatFlux - log(ksHeatFluxSum/alpha)/ksWeight;
    }

  return ksHeatFlux;
}

/*
  Retrieve the reference value
*/
TacsScalar TACSKSHeatFlux::getReferenceHeatFlux(){
  return refHeatFlux;
}


/*
  Initialize the internal values stored within the KS function
*/
void TACSKSHeatFlux::initEvaluation( EvaluationType ftype ){
  if (ftype == TACSFunction::INITIALIZE){
    if (refType == MAX){
      refHeatFlux = -1e20;
    }
    else if (refType == MIN){
      refHeatFlux = 1e20;
    }

  
  }
  else if (ftype == TACSFunction::INTEGRATE){
    ksHeatFluxSum = 0.0;
  }
}

/*
  Reduce the function values across all MPI processes
*/
void TACSKSHeatFlux::finalEvaluation( EvaluationType ftype ){
  if (ftype == TACSFunction::INITIALIZE){
    // Distribute the values of the KS function computed on this domain
    TacsScalar HeatFlux = refHeatFlux;

    if (refType == MAX){
      MPI_Allreduce(&HeatFlux, &refHeatFlux, 1, TACS_MPI_TYPE,
                    TACS_MPI_MAX, assembler->getMPIComm());
    }
    else if (refType == MIN){
      MPI_Allreduce(&HeatFlux, &refHeatFlux, 1, TACS_MPI_TYPE,
                    TACS_MPI_MIN, assembler->getMPIComm());
    }

  }
  else {
    // Find the sum of the ks contributions from all processes
    TacsScalar HeatFlux = ksHeatFluxSum;
    MPI_Allreduce(&HeatFlux, &ksHeatFluxSum, 1, TACS_MPI_TYPE,
                  MPI_SUM, assembler->getMPIComm());

    // Compute the P-norm quantity if needed
    invPnorm = 0.0;
    if (ksType == PNORM_DISCRETE || ksType == PNORM_CONTINUOUS){
      if (ksHeatFluxSum != 0.0){
        invPnorm = pow(ksHeatFluxSum, (1.0 - ksWeight)/ksWeight);
      }
    }
  }
}


/*
  Perform the element-wise evaluation of the TACSDisplacementIntegral function.
*/
void TACSKSHeatFlux::elementWiseEval( EvaluationType ftype,
                                    int elemIndex,
                                    TACSElement *element,
                                    double time,
                                    TacsScalar scale,
                                    const TacsScalar Xpts[],
                                    const TacsScalar vars[],
                                    const TacsScalar dvars[],
                                    const TacsScalar ddvars[] ){
  // Retrieve the number of stress components for this element
  TACSElementBasis *basis = element->getElementBasis();

  
  if (basis){
    // Get the surface index
    int face_key = element_to_face_key[elemIndex];

    for ( int face = 0; (face_key && face < MAX_SURFACE_INDEX); face++ ){
      // Check if this is a surface that we need to integrate
      if (1 << face & face_key){
        // Clear the bit from the surface index we just integrated
        face_key &= ~(1 << face);

        // integrate over the quadrature points to get the value in each element
        for ( int i = 0; i < basis->getNumFaceQuadraturePoints(face); i++ ){
          double pt[3], tangents[6];
          double weight = basis->getFaceQuadraturePoint(face, i, pt, tangents);

          // Evaluate the heat flux at the quadrature point
          TacsScalar flux[3];
          const int not_a_quadrature_pt = -1;
          int count = element->evalPointQuantity(elemIndex, TACS_HEAT_FLUX,
                                                 time, not_a_quadrature_pt, pt,
                                                 Xpts, vars, dvars, ddvars,
                                                 flux);

          // Compute the component of the flux normal to the surface
          if (count > 0){
            TacsScalar X[3], Xd[9], normal[3], heatFlux;
            TacsScalar Area = basis->getFaceNormal(face, i, Xpts, X, Xd, normal);
            
            if (count == 2){
              heatFlux = Area*vec2Dot(flux, normal); // this is actually the heat transfer rate, but don't tell Graeme
            }
            else if (count == 3){
              heatFlux = Area*vec3Dot(flux, normal); 
            }
            
            if (ftype == TACSFunction::INITIALIZE){
            // Set the reference HeatFlux
            
                if (refType == MAX){
                  if (TacsRealPart(heatFlux) > TacsRealPart(refHeatFlux)){
                      refHeatFlux = heatFlux;
                  }
                }
                else if (refType == MIN){
                  if (TacsRealPart(heatFlux) < TacsRealPart(refHeatFlux)){
                      refHeatFlux = heatFlux;
                      // fprintf(stderr, "elem:%d i%d heatflux:%f X: %f, %f, %f,\n", elemIndex, i,  heatFlux, X[0], X[1], X[2]);
                  }
                }
            }
            else {
                // Evaluate the determinant of the Jacobian
                TacsScalar Xd[9], J[9];
                TacsScalar detJ = basis->getJacobianTransform(i, pt, Xpts, Xd, J);

                // fprintf(stderr, "elem:%d i%d heatflux:%f X: %f, %f, %f,\n", elemIndex, i,  heatFlux, X[0], X[1], X[2]);
                // Add the heatFlux to the sum
                if (ksType == DISCRETE){
              
                    TacsScalar fexp;
                    
                    if (refType == MAX){
                      fexp = exp(ksWeight*(heatFlux - refHeatFlux));
              
                    }
                    else if (refType == MIN){
                      fexp = exp(ksWeight*(refHeatFlux - heatFlux));
                    }
                 
                    ksHeatFluxSum += scale*fexp;
                 
                }
                else if (ksType == CONTINUOUS){
                    
                    TacsScalar fexp;
                    
                    if (refType == MAX){
                      fexp = exp(ksWeight*(heatFlux - refHeatFlux));
                    }
                    else if (refType == MIN){
                      fexp = exp(ksWeight*(refHeatFlux - heatFlux));
                    }
                    ksHeatFluxSum += scale*weight*detJ*fexp;
                }
                else if (ksType == PNORM_DISCRETE){
                    TacsScalar fpow = pow(fabs(TacsRealPart(heatFlux/refHeatFlux)), ksWeight);
                    ksHeatFluxSum += scale*fpow;
                }
                else if (ksType == PNORM_CONTINUOUS){
                    TacsScalar fpow = pow(fabs(TacsRealPart(heatFlux/refHeatFlux)), ksWeight);
                    ksHeatFluxSum += scale*weight*detJ*fpow;
                }
            }
          }
        }
      }
    }
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.
*/
void TACSKSHeatFlux::getElementSVSens( int elemIndex, TACSElement *element,
                                     double time,
                                     TacsScalar alpha,
                                     TacsScalar beta,
                                     TacsScalar gamma,
                                     const TacsScalar Xpts[],
                                     const TacsScalar vars[],
                                     const TacsScalar dvars[],
                                     const TacsScalar ddvars[],
                                     TacsScalar dfdu[] ){
  // Zero the derivative of the function w.r.t. the element state
  // variables
  int numVars = element->getNumVariables();
  memset(dfdu, 0, numVars*sizeof(TacsScalar));

  // Get the element basis class
  TACSElementBasis *basis = element->getElementBasis();

  if (basis){
    // Get the surface index
    int face_key = element_to_face_key[elemIndex];

    for ( int face = 0; (face_key && face < MAX_SURFACE_INDEX); face++ ){
      // Check if this is a surface that we need to integrate
      if (1 << face & face_key){
        // Clear the bit from the surface index we just integrated
        face_key &= ~(1 << face);

        for ( int i = 0; i < basis->getNumFaceQuadraturePoints(face); i++ ){
          double pt[3], tangents[6];
          double weight = basis->getFaceQuadraturePoint(face, i, pt, tangents);

          // Evaluate the heat flux at the quadrature point
          TacsScalar flux[3];
          const int not_a_quadrature_pt = -1;
          int count = element->evalPointQuantity(elemIndex, TACS_HEAT_FLUX,
                                                 time, not_a_quadrature_pt, pt,
                                                 Xpts, vars, dvars, ddvars,
                                                 flux);

          if (count > 0){
            TacsScalar X[3], Xd[9], normal[3], heatFlux, J[9];
            TacsScalar Area = basis->getFaceNormal(face, i, Xpts, X, Xd, normal);
            TacsScalar detJ = basis->getJacobianTransform(i, pt, Xpts, Xd, J);

            if (count == 2){
                // this is actually the heat transfer rate, but don't tell Graeme
                heatFlux = Area*vec2Dot(flux, normal);
            }
            else if (count == 3){
                heatFlux = Area*vec3Dot(flux, normal); 
            }
            
            
            TacsScalar dfdq[3] = {0.0, 0.0, 0.0};
                
            // Compute the sensitivity contribution
            TacsScalar dKSHeatFlux_dHeatFlux = 0.0;
            
            if (ksType == DISCRETE){
                // d(log(ksHeatFlux))/dx = 1/(ksHeatFlux)*d(Heatflux)/dx
                                    
                if (refType == MAX){
                  dKSHeatFlux_dHeatFlux = exp(ksWeight*(heatFlux - refHeatFlux))/ksHeatFluxSum;
                }
                else if (refType == MIN){
                  dKSHeatFlux_dHeatFlux = exp(ksWeight*(refHeatFlux - heatFlux))/ksHeatFluxSum;
                }
                
            }
            else if (ksType == CONTINUOUS){
                if (refType == MAX){
                  dKSHeatFlux_dHeatFlux = exp(ksWeight*(heatFlux - refHeatFlux))/ksHeatFluxSum;
                }
                else if (refType == MIN){
                  dKSHeatFlux_dHeatFlux = exp(ksWeight*(refHeatFlux - heatFlux))/ksHeatFluxSum;
                }
                
                dKSHeatFlux_dHeatFlux *= weight*detJ;
            }
            else if (ksType == PNORM_DISCRETE){
                TacsScalar fpow = pow(fabs(TacsRealPart(heatFlux/refHeatFlux)), ksWeight-2.0);
                dKSHeatFlux_dHeatFlux = heatFlux*fpow*invPnorm;
            }
            else if (ksType == PNORM_CONTINUOUS){
                // Get the determinant of the Jacobian
                TacsScalar fpow = pow(fabs(TacsRealPart(heatFlux/refHeatFlux)), ksWeight-2.0);
                dKSHeatFlux_dHeatFlux = heatFlux*fpow*invPnorm;
                dKSHeatFlux_dHeatFlux *= weight*detJ;
            }
                
            // dfdq = dKSHeatFlux_dHeatFlux*dHeatFlux_dFlux
            // f = ksHeatFlux
            // q = flux
            
            if (count == 2){
              dfdq[0] = dKSHeatFlux_dHeatFlux*Area*normal[0];
              dfdq[1] = dKSHeatFlux_dHeatFlux*Area*normal[1];
            }
            else if (count == 3){
              dfdq[0] = dKSHeatFlux_dHeatFlux*Area*normal[0];
              dfdq[1] = dKSHeatFlux_dHeatFlux*Area*normal[1];
              dfdq[2] = dKSHeatFlux_dHeatFlux*Area*normal[2];
            }

            element->addPointQuantitySVSens(elemIndex, TACS_HEAT_FLUX, time,
                                            alpha, beta, gamma, not_a_quadrature_pt, pt,
                                            Xpts, vars, dvars, ddvars,
                                            dfdq, dfdu);
          }
        }
      }
    }
  }
}


/*
  Determine the derivative of the function with respect to
  the element nodal locations
*/
void TACSKSHeatFlux::getElementXptSens( int elemIndex,
                                      TACSElement *element,
                                      double time,
                                      TacsScalar scale,
                                      const TacsScalar Xpts[],
                                      const TacsScalar vars[],
                                      const TacsScalar dvars[],
                                      const TacsScalar ddvars[],
                                      TacsScalar dfdXpts[] ){
  // Zero the derivative of the function w.r.t. the element node
  // locations
  int numNodes = element->getNumNodes();
  memset(dfdXpts, 0, 3*numNodes*sizeof(TacsScalar));

  // Get the element basis class
  TACSElementBasis *basis = element->getElementBasis();

  if (basis){
    // Get the surface index
    int face_key = element_to_face_key[elemIndex];

    for ( int face = 0; (face_key && face < MAX_SURFACE_INDEX); face++ ){
      // Check if this is a surface that we need to integrate
      if (1 << face & face_key){
        // Clear the bit from the surface index we just integrated
        face_key &= ~(1 << face);

        for ( int i = 0; i < basis->getNumFaceQuadraturePoints(face); i++ ){
          double pt[3], tangents[6];
          double weight = basis->getFaceQuadraturePoint(face, i, pt, tangents);

          // Evaluate the heat flux at the quadrature point
          TacsScalar flux[3];
          const int not_a_quadrature_pt = -1;
          int count = element->evalPointQuantity(elemIndex, TACS_HEAT_FLUX,
                                                 time, not_a_quadrature_pt, pt,
                                                 Xpts, vars, dvars, ddvars,
                                                 flux);

          if (count > 0){
            TacsScalar X[3], Xd[9], normal[3], heatFlux, J[9];
            TacsScalar Area = basis->getFaceNormal(face, i, Xpts, X, Xd, normal);
            TacsScalar detJ = basis->getJacobianTransform(i, pt, Xpts, Xd, J);

            if (count == 2){
                // this is actually the heat transfer rate, but don't tell Graeme
                heatFlux = Area*vec2Dot(flux, normal);
            }
            else if (count == 3){
                heatFlux = Area*vec3Dot(flux, normal); 
            }
   

            TacsScalar dfdq[3] = {0.0, 0.0, 0.0};
            TacsScalar dfdn[3] = {0.0, 0.0, 0.0};
            TacsScalar dfdA = 0.0;
                
            // Compute the sensitivity contribution
            TacsScalar dKSHeatFlux_dHeatFlux = 0.0;
            
            if (ksType == DISCRETE){
                // d(log(ksHeatFlux))/dx = 1/(ksHeatFlux)*d(Heatflux)/dx
                if (refType == MAX){
                  dKSHeatFlux_dHeatFlux = exp(ksWeight*(heatFlux - refHeatFlux))/ksHeatFluxSum;
                }
                else if (refType == MIN){
                  dKSHeatFlux_dHeatFlux = exp(ksWeight*(refHeatFlux - heatFlux))/ksHeatFluxSum;
                }
                
            }
            else if (ksType == CONTINUOUS){
                if (refType == MAX){
                  dKSHeatFlux_dHeatFlux = exp(ksWeight*(heatFlux - refHeatFlux))/ksHeatFluxSum;
                }
                else if (refType == MIN){
                  dKSHeatFlux_dHeatFlux = exp(ksWeight*(refHeatFlux - heatFlux))/ksHeatFluxSum;
                }
                dKSHeatFlux_dHeatFlux *= weight*detJ;
            }
            else if (ksType == PNORM_DISCRETE){
                TacsScalar fpow = pow(fabs(TacsRealPart(heatFlux/refHeatFlux)), ksWeight-2.0);
                dKSHeatFlux_dHeatFlux = heatFlux*fpow*invPnorm;
            }
            else if (ksType == PNORM_CONTINUOUS){
                // Get the determinant of the Jacobian
                TacsScalar fpow = pow(fabs(TacsRealPart(heatFlux/refHeatFlux)), ksWeight-2.0);
                dKSHeatFlux_dHeatFlux = heatFlux*fpow*invPnorm;
                dKSHeatFlux_dHeatFlux *= weight*detJ;
            }
                
            // dfdq = dKSHeatFlux_dHeatFlux*dHeatFlux_dFlux
            // f = ksHeatFlux
            // q = flux

            if (count == 2){
              dfdq[0] = dKSHeatFlux_dHeatFlux*Area*normal[0];
              dfdq[1] = dKSHeatFlux_dHeatFlux*Area*normal[1];

              dfdn[0] = dKSHeatFlux_dHeatFlux*scale*Area*flux[0];
              dfdn[1] = dKSHeatFlux_dHeatFlux*scale*Area*flux[1];

              dfdA = dKSHeatFlux_dHeatFlux*scale*vec2Dot(flux, normal);
            }
            else if (count == 3){
              dfdq[0] = dKSHeatFlux_dHeatFlux*Area*normal[0];
              dfdq[1] = dKSHeatFlux_dHeatFlux*Area*normal[1];
              dfdq[2] = dKSHeatFlux_dHeatFlux*Area*normal[2];

              dfdn[0] = dKSHeatFlux_dHeatFlux*scale*Area*flux[0];
              dfdn[1] = dKSHeatFlux_dHeatFlux*scale*Area*flux[1];
              dfdn[2] = dKSHeatFlux_dHeatFlux*scale*Area*flux[2];

              dfdA = dKSHeatFlux_dHeatFlux*scale*vec3Dot(flux, normal);
            }
            element->addPointQuantityXptSens(elemIndex, TACS_HEAT_FLUX, time,
                                             scale, not_a_quadrature_pt, pt,
                                             Xpts, vars, dvars, ddvars,
                                             dfdq, dfdXpts);

            basis->addFaceNormalXptSens(face, i, Area, Xd, normal,
                                        dfdA, NULL, NULL, dfdn, dfdXpts);
          }
        }
      }
    }
  }
}

/*
  Determine the derivative of the function with respect to
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void TACSKSHeatFlux::addElementDVSens( int elemIndex,
                                     TACSElement *element,
                                     double time,
                                     TacsScalar scale,
                                     const TacsScalar Xpts[],
                                     const TacsScalar vars[],
                                     const TacsScalar dvars[],
                                     const TacsScalar ddvars[],
                                     int dvLen, TacsScalar dfdx[] ){
  // Get the element basis class
  TACSElementBasis *basis = element->getElementBasis();

  if (basis){
    // Get the surface index
    int face_key = element_to_face_key[elemIndex];

    for ( int face = 0; (face_key && face < MAX_SURFACE_INDEX); face++ ){
      // Check if this is a surface that we need to integrate
      if (1 << face & face_key){
        // Clear the bit from the surface index we just integrated
        face_key &= ~(1 << face);

        for ( int i = 0; i < basis->getNumFaceQuadraturePoints(face); i++ ){
          double pt[3], tangents[6];
          double weight = basis->getFaceQuadraturePoint(face, i, pt, tangents);

          // Evaluate the heat flux at the quadrature point
          TacsScalar flux[3];
          const int not_a_quadrature_pt = -1;
          int count = element->evalPointQuantity(elemIndex, TACS_HEAT_FLUX,
                                                 time, not_a_quadrature_pt, pt,
                                                 Xpts, vars, dvars, ddvars,
                                                 flux);

          if (count > 0){
            TacsScalar X[3], Xd[9], normal[3], heatFlux, J[9];
            TacsScalar Area = basis->getFaceNormal(face, i, Xpts, X, Xd, normal);
            TacsScalar detJ = basis->getJacobianTransform(i, pt, Xpts, Xd, J);

            if (count == 2){
                // this is actually the heat transfer rate, but don't tell Graeme
                heatFlux = Area*vec2Dot(flux, normal);
            }
            else if (count == 3){
                heatFlux = Area*vec3Dot(flux, normal); 
            }
            
            
            TacsScalar dfdq[3] = {0.0, 0.0, 0.0};
                
            // Compute the sensitivity contribution
            TacsScalar dKSHeatFlux_dHeatFlux = 0.0;
            
            if (ksType == DISCRETE){
                // d(log(ksHeatFlux))/dx = 1/(ksHeatFlux)*d(Heatflux)/dx
                dKSHeatFlux_dHeatFlux = exp(ksWeight*(heatFlux - refHeatFlux))/ksHeatFluxSum;
            }
            else if (ksType == CONTINUOUS){
                dKSHeatFlux_dHeatFlux = exp(ksWeight*(heatFlux - refHeatFlux))/ksHeatFluxSum;
                dKSHeatFlux_dHeatFlux *= weight*detJ;
            }
            else if (ksType == PNORM_DISCRETE){
                TacsScalar fpow = pow(fabs(TacsRealPart(heatFlux/refHeatFlux)), ksWeight-2.0);
                dKSHeatFlux_dHeatFlux = heatFlux*fpow*invPnorm;
            }
            else if (ksType == PNORM_CONTINUOUS){
                // Get the determinant of the Jacobian
                TacsScalar fpow = pow(fabs(TacsRealPart(heatFlux/refHeatFlux)), ksWeight-2.0);
                dKSHeatFlux_dHeatFlux = heatFlux*fpow*invPnorm;
                dKSHeatFlux_dHeatFlux *= weight*detJ;
            }
                
            // dfdq = dKSHeatFlux_dHeatFlux*dHeatFlux_dFlux
            // f = ksHeatFlux
            // q = flux
            
            if (count == 2){
              dfdq[0] = dKSHeatFlux_dHeatFlux*Area*normal[0];
              dfdq[1] = dKSHeatFlux_dHeatFlux*Area*normal[1];
            }
            else if (count == 3){
              dfdq[0] = dKSHeatFlux_dHeatFlux*Area*normal[0];
              dfdq[1] = dKSHeatFlux_dHeatFlux*Area*normal[1];
              dfdq[2] = dKSHeatFlux_dHeatFlux*Area*normal[2];

            element->addPointQuantityDVSens(elemIndex, TACS_HEAT_FLUX,
                                            time, scale, not_a_quadrature_pt, pt,
                                            Xpts, vars, dvars, ddvars,
                                            dfdq, dvLen, dfdx);
          }
        }
        }
      }
    }
  } 
}