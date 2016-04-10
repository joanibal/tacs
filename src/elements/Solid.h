#ifndef TACS_SOLID_H
#define TACS_SOLID_H

/*
  Solid element implementation

  Copyright (c) 2010-2015 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.

  The following code uses templates to allow for arbitrary order elements.
*/

#include "TACS3DElement.h"
#include "FElibrary.h"

template <int order>
class Solid : public TACS3DElement<order*order*order> {
 public:
  enum SolidElementType { LINEAR, NONLINEAR };
  
  Solid( SolidStiffness * _stiff, 
	 SolidElementType type = LINEAR, 
	 int _componentNum = 0 ); 
  ~Solid();

  // Return the name of this element
  // -------------------------------
  const char * elementName(){ return elemName; }

  // Retrieve the shape functions
  // ----------------------------
  void getShapeFunctions( const double pt[], double N[] );
  void getShapeFunctions( const double pt[], double N[],
			  double Na[], double Nb[], double Nc[] );

  // Retrieve the Gauss points/weights
  // ---------------------------------
  int getNumGaussPts();
  double getGaussWtsPts( const int num, double pt[] ); 

  // Functions for post-processing
  // -----------------------------
  void addOutputCount( int * nelems, int * nnodes, int * ncsr );
  void getOutputData( unsigned int out_type, 
		      double * data, int ld_data, 
		      const TacsScalar Xpts[],
		      const TacsScalar vars[] );
  void getOutputConnectivity( int * con, int node );

 private:
  // The number of nodes in the element
  static const int NUM_NODES = order*order*order;

  // The Gauss quadrature scheme
  int numGauss;
  const double *gaussWts, *gaussPts;

  // Store the name of the element
  static const char * elemName;
};

template <int order>
Solid<order>::Solid( SolidStiffness * _stiff, 
		     SolidElementType type, 
		     int _componentNum ):
TACS3DElement<order*order*order>(_stiff, (type == LINEAR), _componentNum){
  numGauss = FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);
} 

template <int order>
Solid<order>::~Solid(){}

template <int order>
const char * Solid<order>::elemName = "Solid";

/*
  Get the number of Gauss points in the scheme
*/
template <int order>
int Solid<order>::getNumGaussPts(){
  return numGauss*numGauss*numGauss;
}

/*
  Get the Gauss points
*/
template <int order>
double Solid<order>::getGaussWtsPts( int npoint, double pt[] ){
  // Compute the n/m/p indices of the Gauss quadrature scheme
  int p = (int)((npoint)/(numGauss*numGauss));
  int m = (int)((npoint - numGauss*numGauss*p)/numGauss);
  int n = npoint - numGauss*m - numGauss*numGauss*p;
  
  pt[0] = gaussPts[n];
  pt[1] = gaussPts[m];    
  pt[2] = gaussPts[p];
  
  return gaussWts[n]*gaussWts[m]*gaussWts[p];
}

/*
  Evaluate the shape functions and their derivatives
*/
template <int order>
void Solid<order>::getShapeFunctions( const double pt[],
				      double N[] ){
  FElibrary::triLagrangeSF(N, pt, order);
}

/*
  Compute the shape functions and their derivatives w.r.t. the
  parametric element location 
*/
template <int order>
void Solid<order>::getShapeFunctions( const double pt[], double N[],
				      double Na[], double Nb[],
				      double Nc[] ){
  FElibrary::triLagrangeSF(N, Na, Nb, Nc, pt, order);
}

/*
  Get the number of elemens/nodes and CSR size of the contributed by
  this element.  
*/
template <int order>
void Solid<order>::addOutputCount( int *nelems, 
				   int *nnodes, int *ncsr ){
  *nelems += (order-1)*(order-1)*(order-1);
  *nnodes += order*order*order;
  *ncsr += 8*(order-1)*(order-1)*(order-1);
}

/*
  Get the output data from this element and place it in a real
  array for visualization later. The values generated for visualization
  are determined by a bit-wise selection variable 'out_type' which is 
  can be used to simultaneously write out different data. Note that this
  is why the bitwise operation & is used below. 

  The output may consist of the following:
  - the nodal locations
  - the displacements and rotations
  - the strains or strains within the element
  - extra variables that are used for optimization
  
  output:
  data:     the data to write to the file (eventually)

  input:
  out_type: the bit-wise variable used to specify what data to generate
  vars:     the element variables
  Xpts:     the element nodal locations
*/
template <int order>
void Solid<order>::getOutputData( unsigned int out_type, 
				  double *data, int ld_data,
				  const TacsScalar Xpts[],
				  const TacsScalar vars[] ){
  for ( int p = 0; p < order; p++ ){
    for ( int m = 0; m < order; m++ ){
      for ( int n = 0; n < order; n++ ){
	int node = n + m*order + p*order*order;
	int index = 0;
	if (out_type & TACSElement::OUTPUT_NODES){
	  for ( int k = 0; k < 3; k++ ){
	    data[index+k] = RealPart(Xpts[3*node+k]);
	  }
          index += 3;
	}
	if (out_type & TACSElement::OUTPUT_DISPLACEMENTS){
	  for ( int k = 0; k < 3; k++ ){
	    data[index+k] = RealPart(vars[3*node+k]);
	  }
	  index += 3;
	}

	// Set the parametric point where to evaluate the 
	// stresses/strains
        double pt[3];
        pt[0] = -1.0 + 2.0/(order-1.0)*n;
        pt[1] = -1.0 + 2.0/(order-1.0)*m;
        pt[2] = -1.0 + 2.0/(order-1.0)*p;
        
	// Compute the shape functions
	double N[NUM_NODES];
	double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
	getShapeFunctions(pt, N, Na, Nb, Nc);

	// Compute the derivative of X with respect to the
	// coordinate directions
	TacsScalar X[3], Xa[9];
	this->solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

	// Compute the determinant of Xa and the transformation
	TacsScalar J[9];
	FElibrary::jacobian3d(Xa, J);

	// Compute the strain
	TacsScalar strain[6];
	this->evalStrain(strain, J, Na, Nb, Nc, vars);
	
        if (out_type & TACSElement::OUTPUT_STRAINS){
          for ( int k = 0; k < 6; k++ ){
            data[index+k] = RealPart(strain[k]);
          }
          index += 6;
        }
        if (out_type & TACSElement::OUTPUT_STRESSES){
	  TacsScalar stress[6];
	  this->stiff->calculateStress(pt, strain, stress);
          
          for ( int k = 0; k < 6; k++ ){
            data[index+k] = RealPart(stress[k]);
          }
          index += 6;
        }
        if (out_type & TACSElement::OUTPUT_EXTRAS){
	  // Compute the failure value
	  TacsScalar lambda;
	  this->stiff->failure(pt, strain, &lambda);
	  data[index] = RealPart(lambda);

	  this->stiff->buckling(strain, &lambda);
	  data[index+1] = RealPart(lambda);

	  data[index+2] = RealPart(this->stiff->getDVOutputValue(0, pt));
	  data[index+3] = RealPart(this->stiff->getDVOutputValue(1, pt));
	  
	  index += this->NUM_EXTRAS;
        }

	data += ld_data;
      }
    }
  }
}

/*
  Get the element connectivity for visualization purposes. Since each
  element consists of a series of sub-elements used for visualization,
  we also need the connectivity of these visualization elements.

  output:
  con:  the connectivity of the local visualization elements contributed
  by this finite-element

  input:
  node:  the node offset number - so that this connectivity is more or 
  less global
*/
template <int order>
void Solid<order>::getOutputConnectivity( int *con, int node ){
  int j = 0;
  for ( int p = 0; p < order-1; p++ ){
    for ( int m = 0; m < order-1; m++ ){
      for ( int n = 0; n < order-1; n++ ){
	con[8*j]   = node + n   +     m*order +     p*order*order;
	con[8*j+1] = node + n+1 +     m*order +     p*order*order;
	con[8*j+2] = node + n+1 + (m+1)*order +     p*order*order;
	con[8*j+3] = node + n   + (m+1)*order +     p*order*order;
	con[8*j+4] = node + n   +     m*order + (p+1)*order*order;
	con[8*j+5] = node + n+1 +     m*order + (p+1)*order*order;
	con[8*j+6] = node + n+1 + (m+1)*order + (p+1)*order*order;
	con[8*j+7] = node + n   + (m+1)*order + (p+1)*order*order;
	j++;
      }
    }
  }
}

#endif // TACS_SOLID_H