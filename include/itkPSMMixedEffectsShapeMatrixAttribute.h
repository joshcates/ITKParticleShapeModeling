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
#ifndef __itkPSMMixedEffectsShapeMatrixAttribute_h
#define __itkPSMMixedEffectsShapeMatrixAttribute_h

#include "itkPSMRegressionShapeMatrixAttribute.h"
#include "vnl/vnl_trace.h"

namespace itk
{
/** \class PSMMixedEffectsShapeMatrixAttribute
 *
 *
 *
 */
template <class T, unsigned int VDimension>
class ITK_EXPORT PSMMixedEffectsShapeMatrixAttribute
  : public PSMRegressionShapeMatrixAttribute<T,VDimension>
{
public:
  /** Standard class typedefs */
  typedef T DataType;
  typedef PSMMixedEffectsShapeMatrixAttribute Self;
  typedef PSMRegressionShapeMatrixAttribute<T,VDimension> Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  typedef WeakPointer<const Self>  ConstWeakPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(PSMMixedEffectsShapeMatrixAttribute, PSMRegressionShapeMatrixAttribute);
  
  /** */
  virtual void UpdateMeanMatrix();
  
  // CAN BE INHERITED -- Is this the correct mean?
  //  inline vnl_vector<double> ComputeMean(double k) const
  // {
  //   return m_Intercept + m_Slope * k;    
  // }
  
  /** Set / Get the number of timepoints per individual in the cohort.
      A current limitation of this implementation is that all
      individuals must have the same number of time points.*/
  itkSetMacro(TimePointsPerIndividual, int);
  itkGetConstMacro(TimePointsPerIndividual, int);
  
  /** */
  const vnl_matrix<double> &GetSlopeRandom() const
  { return m_SlopeRand; }

  /** */
  const vnl_matrix<double> &GetInterceptRandom() const
  { return m_InterceptRand; }
    
  /** */
  virtual void EstimateParameters();
  
  /** */
  virtual void Initialize()
  {
    Superclass::Initialize();

    // Initialize variables specific to this subclass
    m_SlopeRand.fill(0.0);
    m_InterceptRand.fill(0.0);    
  }

 protected:
  PSMMixedEffectsShapeMatrixAttribute() 
    {
      m_NumberOfIndividuals     = 0;
      m_TimePointsPerIndividual = 1;
    }
  virtual ~PSMMixedEffectsShapeMatrixAttribute() {};
  
  void PrintSelf(std::ostream& os, Indent indent) const
  { Superclass::PrintSelf(os,indent);  }
  
 private:
  PSMMixedEffectsShapeMatrixAttribute(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  vnl_matrix<double> m_InterceptRand; //added: AK , random intercepts for each group
  vnl_matrix<double> m_SlopeRand; //added: AK , random slopes for each group
  int m_NumberOfIndividuals;
  int m_TimePointsPerIndividual;
 };
 
} // end namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPSMMixedEffectsShapeMatrixAttribute.hxx"
#endif


#endif
