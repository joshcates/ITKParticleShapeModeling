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
#ifndef __itkPSMProcrustesFunction_h
#define __itkPSMProcrustesFunction_h

#include <vector>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkObject.h"
#include "itkObjectFactory.h"

namespace itk
{

/**
   JOSH -- need Doxygen documentation here describing this struct.  Also, this
   struct should probably just be embedded in the PSMProcrustes function
   itself, since this class is only used within PSMProcrustesFunction.
*/
struct SimilarityTransform3D
{
  vnl_matrix_fixed<double, 3, 3> rotation;
  double scale;
  vnl_vector_fixed<double, 3> translation;
};

/**
   JOSH -- need Doxygen documentation here. See other classes for examples.
   
*/
class PSMProcrustesFunction : public itk::Object
{
public:
  typedef double                                RealType;
  typedef vnl_vector_fixed<double, 3>           PointType;
  typedef std::vector<PointType>                ShapeType;
  typedef ShapeType::iterator                   ShapeIteratorType;
  typedef std::vector<ShapeType>                ShapeListType;
  typedef ShapeListType::iterator               ShapeListIteratorType;
  typedef std::vector<SimilarityTransform3D>    SimilarityTransformListType;
  typedef SimilarityTransformListType::iterator SimilarityTransformListIteratorType;
  
  /** Standard class typedefs. */
  typedef PSMProcrustesFunction                 Self;
  typedef Object                                Superclass;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PSMProcrustesFunction, Object);

  // JOSH -- changed comment syntax for consistency
  /** Align a list of shapes using Generalized Procrustes Analysis */
  void RunGeneralizedProcrustes(SimilarityTransformListType & transform,
                                ShapeListType & shapes);

  /** JOSH -- need description here */
  RealType ComputeSumOfSquares(ShapeListType & shapes);

  /** Helper function to transform a shape by a similarity transform */
  ShapeType TransformShape(ShapeType shape, SimilarityTransform3D & transform);
  
private:        
  /** Check if shapes are the same */
  bool CheckDegenerateCase(PointType ssqShape, PointType ssqMean, PointType muShape,
                           PointType muMean, int rows);

  /** Compute mean of all shapes except one */
  void LeaveOneOutMean(ShapeType & mean, ShapeListType & shapeList,
                       ShapeListIteratorType & leaveOutIt);

  /** Align two shapes (translation, rotation & scale) using ordinary Procrustes
      Analysis */
  ShapeType RunProcrustes(SimilarityTransform3D & transform, ShapeType mean,
                          ShapeListIteratorType & leaveOutIt);

  // TODO: Template the class to allow for N-D
};

} // end namespace itk
#endif
