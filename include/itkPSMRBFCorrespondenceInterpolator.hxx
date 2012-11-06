/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkPSMRBFCorrespondenceInterpolator_hxx
#define __itkPSMRBFCorrespondenceInterpolator_hxx
#include "itkPSMRBFCorrespondenceInterpolator.h"

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

namespace itk
{

template <unsigned int VDimension>
PSMRBFCorrespondenceInterpolator<VDimension>
::PSMRBFCorrespondenceInterpolator()
{

  //
  // TODO:  NOTE ! This function is currently only implemented for 3D
  //
  if (VDimension != 3)
    {
      itkExceptionMacro("This function is currently only implemented for 3D point sets.");
    }

}

template <unsigned int VDimension>
void
PSMRBFCorrespondenceInterpolator<VDimension>
::Initialize()
{ 
  // Make sure we have points on the input
  if (m_PointSetA.size() == 0 || m_PointSetB.size() == 0)
    {
      itkExceptionMacro("Surface point sets have not been specified.");
    }
  
  // Make sure the point set sizes match
  if (m_PointSetA.size() != m_PointSetB.size())
    {
      itkExceptionMacro("Surface point sets A and B are not the same size.");
    }

  // N is the number of points
  const unsigned int N = m_PointSetA.size();

  // TODO: Need algorithm references
  // TODO: Need to implement for at least 2D

  //  INSERT CODE HERE 

  typedef vnl_vector<double> vectype;
  typedef vnl_matrix<double> matrixtype;
  
  // Transforms are from 1 to 2
  vectype b(N+4); // point data from file 2
  vectype x(N+4); // parameters to be solved for
  matrixtype A(N+4,N+4);
  matrixtype Phi(N,N);
  
  // Compute Phi values
  vectype Xi(3);
  vectype Xj(3);
  for (unsigned int i = 0; i < N; i++)
    {
    Xi[0] = (m_PointSetA[i])[0];
    Xi[1] = (m_PointSetA[i])[1];
    Xi[2] = (m_PointSetA[i])[2];
    
    for (unsigned int j = 0; j < N; j++)
      {
      Xj[0] = (m_PointSetB[j])[0];
      Xj[1] = (m_PointSetB[j])[1];
      Xj[2] = (m_PointSetB[j])[2];
      
      // Biharmonic Radial Basis Function
      Phi(i,j) = (Xi - Xj).magnitude();

      // Thin plate spline basis -- probably need some logic to switch on dimensionality
      // Phi(i,j) = dot_product(Xi - Xj,Xi - Xj) * log((Xi - Xj).magnitude() + 1.0e-6);
      }
    }
 
  // Construct the coefficient matrix A from PointSetA and Phi's
  for (unsigned int i = 0; i < 4; i++)
    {
    for (unsigned int j = 0; j < N +4; j++)
      {
      if      (j >= N) A(i,j) = 0.0;
      else if (i == 3) A(i,j) = 1.0;
      else             A(i,j) = (m_PointSetA[j])[i];
      }    
    }

  // Rest of the matrix A
  for (unsigned int i = 0; i < N; i++)
    {
    for (unsigned int j = 0; j < N +4; j++)
      {
      if      (j < N)    A(i+4,j) = Phi(i,j);
      else if (j == N+3) A(i+4,j) = 1.0;
      else               A(i+4,j) = (m_PointSetA[i])[2-(j-N)];
      }
    }

  // To keep things simple, solve in each dimension separately
  matrixtype P(4,3);
  matrixtype C(N,3);
  for (unsigned int D = 0; D < 3; D++)
    {
    // Construct the position data vector b from file 2 
    // [0 0 0 0 x1' x2' x3' ... xN' ]^T   
    b[0] = b[1] = b[2] = b[3] = 0.0;
    for(unsigned int i = 0; i < N; i++)
      {
      b[i+4] = (m_PointSetB[i])[D];
      }

    // Solve for parameters
    x = vnl_matrix_inverse<double>(A) * b;

    // Store results
    for (unsigned int i = 0; i < N; i++)
      {
      C(i,D) = x[i];
      }
    
    P(0,D) = x[N+3];
    P(1,D) = x[N+2];
    P(2,D) = x[N+1];
    P(3,D) = x[N];
    
    } // end D

  m_Initialized = true;
}
  
template <unsigned int VDimension>
typename PSMRBFCorrespondenceInterpolator<VDimension>::VectorType
PSMRBFCorrespondenceInterpolator<VDimension>
::Evaluate() const
{
  if (m_Initialized == false)
    {
      itkExceptionMacro("Function has not been initialized");
    }
  

}

} // end namespace itk
#endif
