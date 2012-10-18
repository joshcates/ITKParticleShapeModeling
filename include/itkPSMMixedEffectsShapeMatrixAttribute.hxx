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
#ifndef __itkPSMMixedEffectsShapeMatrixAttribute_hxx
#define __itkPSMMixedEffectsShapeMatrixAttribute_hxx
#include "itkPSMMixedEffectsShapeMatrixAttribute.h"

namespace itk
{

template <class T, unsigned int VDimension>
void PSMMixedEffectsShapeMatrixAttribute<T,VDimension>
::UpdateMeanMatrix()
{
  vnl_vector<double> tempvect;
  tempvect.set_size(this->GetMeanMatrix().rows());
  const vnl_vector<double> expl = this->GetExplanatory();

  // For each sample, which are in the columns, ...
  for (unsigned int i = 0; i < this->GetMeanMatrix().cols(); i++)
    {
      int group_indx = i / m_TimePointsPerIndividual;
      tempvect = this->GetIntercept() + this->GetSlope() * expl(i);
      tempvect = tempvect + m_InterceptRand.get_row(group_indx);
      tempvect = tempvect + m_SlopeRand.get_row(group_indx) * expl(i);
      // ... compute the mean.
      this->GetMeanMatrix().set_column(i, tempvect);
    }
}
  
template <class T, unsigned int VDimension>
void PSMMixedEffectsShapeMatrixAttribute<T,VDimension>::
EstimateParameters()
{
  vnl_matrix<double> X = *this + this->GetMeanMatrix();
  const vnl_vector<double> expl = this->GetExplanatory();
    
  // Number of samples
  int num_shapes = static_cast<double>(X.cols());
  this->m_NumberOfIndividuals = num_shapes / (this->GetDomainsPerShape()*this->GetTimePointsPerIndividual());
  int nr = X.rows(); //number of points*3
  
  //set the sizes of random slope and intercept matrix
  m_SlopeRand.set_size(m_NumberOfIndividuals, nr); //num_groups X num_points*3
  m_InterceptRand.set_size(m_NumberOfIndividuals, nr); //num_groups X num_points*3
  
  vnl_matrix<double> fixed; //slopes + intercepts for all points
  vnl_matrix<double> random; //slopes + intercepts for all groups, for all points
  fixed.set_size(2, nr);
  random.set_size(2, nr*m_NumberOfIndividuals);
  vnl_matrix<double> Ds(2,2); //covariance matrix of random parameters (2x2)
  Ds.set_identity();  //initialize to identity
  double sigma2s = 1;  //variance of error
  vnl_matrix<double> identity_n;
  identity_n.set_size(m_TimePointsPerIndividual, m_TimePointsPerIndividual);
  identity_n.set_identity();
  vnl_matrix<double> identity_2;
  identity_2.set_size(2,2);
  identity_2.set_identity();
  vnl_matrix<double> *Ws=NULL, *Vs=NULL;
  Ws = new vnl_matrix<double>[m_NumberOfIndividuals];
  Vs = new vnl_matrix<double>[m_NumberOfIndividuals];
  for (int i = 0; i < m_NumberOfIndividuals; i++)
    {
      Vs[i].set_size(m_TimePointsPerIndividual, m_TimePointsPerIndividual); 
      Ws[i].set_size(m_TimePointsPerIndividual, m_TimePointsPerIndividual);
    }
  
  vnl_matrix<double> sum_mat1(2,2,0);
  vnl_vector<double> sum_mat2(2); 
  sum_mat2.fill(0.0);
  
  vnl_vector<double> residual; 
  residual.set_size(m_TimePointsPerIndividual); 
  residual.fill(0.0);
  
  double ecorr = 0.0;
  double tracevar = 0.0;
  
  vnl_matrix<double> bscorr(2,2,0.0); 
  vnl_matrix<double> bsvar(2,2,0.0);
  vnl_matrix<double> Xp; 
  Xp.set_size(m_TimePointsPerIndividual, 2);
  vnl_vector<double> y; 
  y.set_size(m_TimePointsPerIndividual);
  vnl_vector<double> tempvect; 
  tempvect.set_size(2);
  
  for (int i = 0; i < nr; i++) //for all points (x,y,z coordinates)
    {
      sigma2s = 1.0;
      Ds.set_identity();
      for (int j = 0; j < 50; j++) //EM iterations
        {
          sum_mat1.fill(0.0); 
          sum_mat2.fill(0.0);
          residual.fill(0.0);
          ecorr = 0.0; 
          tracevar = 0.0;
          bscorr.fill(0.0);
          bsvar.fill(0.0);
          
          for (int k = 0; k < m_NumberOfIndividuals; k++)
            {
              for (int l = 0; l < m_TimePointsPerIndividual; l++)
                {
                  Xp(l,0) = expl(k*m_TimePointsPerIndividual + l);
                  Xp(l,1) = 1;
                  y(l) = X(i, k*m_TimePointsPerIndividual + l);
                }
              Vs[k] = (identity_n * sigma2s) + Xp * Ds * vnl_transpose(Xp);
              //Ws = static_cast<vnl_matrix> (vnl_matrix_inverse<double>(Vs));
              Ws[k] = vnl_inverse(Vs[k]);
              sum_mat1 = sum_mat1 + vnl_transpose(Xp) * Ws[k] * Xp;
              sum_mat2 = sum_mat2 + vnl_transpose(Xp) * Ws[k] * y;
            }
          tempvect = vnl_inverse(sum_mat1) * sum_mat2; 
          fixed.set_column(i, tempvect);
          for (int k = 0; k < m_NumberOfIndividuals; k++)
            {
              for (int l = 0; l < m_TimePointsPerIndividual; l++)
                {
                  Xp(l,0) = expl(k*m_TimePointsPerIndividual + l);
                  Xp(l,1) = 1;
                  y(l) = X(i, k*m_TimePointsPerIndividual + l);
                }
              tempvect = Ds * vnl_transpose(Xp) * Ws[k] * (y - (Xp * fixed.get_column(i)));
              random.set_column(i*m_NumberOfIndividuals + k, tempvect);
              residual = y - (Xp * fixed.get_column(i)) - (Xp * random.get_column(i*m_NumberOfIndividuals + k));
              ecorr = ecorr + dot_product(residual, residual);
              tracevar = tracevar + (m_TimePointsPerIndividual - sigma2s * vnl_trace(Ws[k]));
              bscorr = bscorr + outer_product(random.get_column(i*m_NumberOfIndividuals + k), random.get_column(i*m_NumberOfIndividuals + k));
              bsvar = bsvar + (identity_2 - (vnl_transpose(Xp) * Ws[k] * Xp * Ds));
            }
          sigma2s = (ecorr + sigma2s * tracevar) / (num_shapes);
          Ds = (bscorr + Ds * bsvar) / m_NumberOfIndividuals;
        }//endfor EM iterations
      //printf ("point #%d\n", i);
    }//endfor all points on shape (x,y & z)
  
  this->SetSlope(fixed.get_row(0));
  this->SetIntercept(fixed.get_row(1));
  for (int i = 0; i < m_NumberOfIndividuals; i++)
    {
      for (int j = 0; j < nr; j++) //for all points * 3
        {
          m_SlopeRand(i,j) = random(0, j*m_NumberOfIndividuals+i);
          m_InterceptRand(i,j) = random(1, j*m_NumberOfIndividuals+i);
        }
    }
  delete [] Vs;
  delete [] Ws;

  //printf ("points:\n");
  //for (int k = 0; k < m_NumberOfIndividuals; k++)
  //	for (int l = 0; l < m_TimePointsPerIndividual; l++)
  //		printf ("%g   %g\n", X(0,k*m_TimePointsPerIndividual + l), expl(k*m_TimePointsPerIndividual + l));
  
  //printf ("fixed: slope %g, intercept %g", m_Slope(0), this->GetIntercept()(0));
  //printf ("random: slopes %g %g, intercepts %g %g", m_SlopeRand(0,0), m_SlopeRand(1,0), m_InterceptRand(0,0), m_InterceptRand(1,0));  
}

} // end namespace itk

#endif
