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

#include "itkPSMProcrustesFunction.h"
#include <iostream>
#include <vnl/algo/vnl_svd.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iostream>
#include <fstream>

namespace itk
{

void PSMProcrustesFunction
::RunGeneralizedProcrustes(SimilarityTransformListType & transforms,
                           ShapeListType & shapes)
{
  ShapeListIteratorType leaveOutIt;
  SimilarityTransformListIteratorType transformIt;
  ShapeIteratorType shapeIt, meanIt;
  ShapeType shape, mean;
  SimilarityTransform3D transform;
  PointType ssqShape, ssqMean;

  const RealType SOS_EPSILON = 1.0e-8;

  int numOfShapes = shapes.size();
  
  // Initialize transform structure
  transform.rotation.set_identity();
  transform.scale = 1.0;
  transform.translation.fill(0.0);
  
  // Initialize transforms vector
  transforms.clear();
  transforms.reserve(shapes.size());
  for(int i = 0; i<numOfShapes; i++)
    {
    transforms.push_back(transform);
    }
  
  RealType sumOfSquares = ComputeSumOfSquares(shapes);
  RealType newSumOfSquares, diff = 1e10;
  
  int counter = 1;
  while(diff > SOS_EPSILON)
    {
    // Initialize ssqShape vector
    ssqShape.fill(0.0);
    // Initialize ssqMean vector
    ssqMean.fill(0.0);
    transformIt = transforms.begin();
    
    for(leaveOutIt = shapes.begin(); leaveOutIt != shapes.end(); leaveOutIt++)
      {
      // Calculate mean of all shapes but one
      LeaveOneOutMean(mean, shapes, leaveOutIt);
      (*leaveOutIt) = RunProcrustes((*transformIt), mean, leaveOutIt);
      transformIt++;
      } // End shape list iteration
    
    // Fix scalings so geometric average = 1
    RealType scaleAve = 0.0;
    for(transformIt = transforms.begin(); transformIt != transforms.end(); transformIt++)
      scaleAve += log((*transformIt).scale);
    
    scaleAve = exp(scaleAve / static_cast<RealType>(transforms.size()));
    
    SimilarityTransform3D scaleSim;
    scaleSim.rotation.set_identity();
    scaleSim.translation.fill(0.0);
    scaleSim.scale = 1.0 / scaleAve;
    
    ShapeListIteratorType shapeListIt = shapes.begin();
    transformIt = transforms.begin();
    while(shapeListIt != shapes.end())
      {
      TransformShape((*shapeListIt), scaleSim);
      (*transformIt).scale /= scaleAve;
      
      shapeListIt++;
      transformIt++;
      }
    
    newSumOfSquares = ComputeSumOfSquares(shapes);
    diff = sumOfSquares - newSumOfSquares;
    
    sumOfSquares = newSumOfSquares;
    counter++;

    // JOSH -- should throw an exception here
    if(counter > 1000)
      {
      std::cout << "Number of iterations on shapes is too high. " << std::endl;
      break;
      }
    } // End while loop
  //std::cout <<  "counter: " << counter << std::endl;
}

PSMProcrustesFunction::ShapeType PSMProcrustesFunction
::RunProcrustes(SimilarityTransform3D & transform, ShapeType mean,
                ShapeListIteratorType & leaveOutIt)
{
  ShapeIteratorType shapeIt1, shapeIt2;
  SimilarityTransform3D newTransform;
  ShapeType shapeScaled, meanScaled;
  PointType colMeanShape, colMeanMean, ssqShape, ssqMean;
  double normMean, normShape;
  vnl_matrix_fixed<RealType, 3, 3> shapeMat;
  shapeMat.fill(0.0);
  
  int numPoints = (*leaveOutIt).size();
  int dim = (*leaveOutIt)[0].size();

  // This needs to be dynamically allocated since "dim" is not a constant
  // double meanScaledTranspose[dim][numPoints];
  vnl_matrix<double> meanScaledTranspose(dim,numPoints);

  // Initialize variables
  colMeanShape.fill(0.0);
  colMeanMean.fill(0.0);
  normMean = 0;
  normShape = 0;
  ssqShape.fill(0.0);
  ssqMean.fill(0.0);

  // Centering the shapes at the origin
  // First calculate mean along columns
  for(int j = 0; j<dim; j++)
    {
    for(int i = 0; i<numPoints; i++)
      {
      colMeanShape[j] += (*leaveOutIt)[i][j];
      colMeanMean[j] += mean[i][j];
      }
    colMeanShape[j] = colMeanShape[j]/numPoints;
    colMeanMean[j] = colMeanMean[j]/numPoints;
    }
  
  // Repeat rows to create new vector
  ShapeIteratorType it1 = (*leaveOutIt).begin();
  ShapeIteratorType it2 = mean.begin();
  
  while(it1 != (*leaveOutIt).end())
    {
    shapeScaled.push_back((*it1) - colMeanShape);
    meanScaled.push_back((*it2) - colMeanMean);
    it1++;
    it2++;
    }
  
  // Calculate sum of squared elements of shapeScaled and meanScaled vectors along columns
  for(int j = 0; j<dim; j++)
    {
    for(int i = 0; i<numPoints; i++)
      {
      ssqShape[j] += (shapeScaled[i][j] * shapeScaled[i][j]);
      ssqMean[j] += (meanScaled[i][j] * meanScaled[i][j]);
      }
    }
  // TODO: Check if dimensions match?
  // Check if shapes are the same
  bool constShape = CheckDegenerateCase(ssqShape, ssqMean, colMeanShape, colMeanMean, numPoints);
  
  // Continue iterations
  if(constShape)
    {
    // Calculate scale normalizing value
    for(int j = 0; j<dim; j++)
      {
      normShape += ssqShape[j];
      normMean += ssqMean[j];
      }
    normShape = sqrt(normShape);
    normMean = sqrt(normMean);
    // Scale shapes to equal (unit) norm
    ShapeIteratorType shapeScaledIt = shapeScaled.begin();
    ShapeIteratorType meanScaledIt = meanScaled.begin();
    while(shapeScaledIt != shapeScaled.end())
      {
      (*shapeScaledIt) = (*shapeScaledIt)/normShape;
      (*meanScaledIt) = (*meanScaledIt)/normMean;
      shapeScaledIt++;
      meanScaledIt++;
      }
    }
  // The degenerate cases: both shapes are the same
  else
    {
    ShapeType output;
    for(int i = 0; i < numPoints; i++)
      output.push_back(vnl_vector_fixed<RealType, 3>(0.0, 0.0, 0.0));
    ShapeIteratorType outputIt = output.begin();
    while(outputIt != output.end())
      {
      (*outputIt) = colMeanShape;
      outputIt++;
      }
    transform.scale = 1.0;
    transform.rotation.set_identity();
    outputIt = output.begin();
    transform.translation = (*outputIt);
    }
  
  for(int j = 0; j<dim; j++)
    {
    for(int i = 0; i<numPoints; i++)
      {
      meanScaledTranspose[j][i] = meanScaled[i][j];
      }
    }
  
  // Build shapeMat = meanScaledTranspose * shapeScaled
  for(int i = 0; i<dim; i++)
    {
    for(int j = 0; j<dim; j++)
      {
      for(int k = 0; k<numPoints; k++)
        {
        shapeMat(i, j) += meanScaledTranspose[i][k] * shapeScaled[k][j];
        }
      }
    }
  
  // Calculate SVD
  vnl_svd<RealType> svd(shapeMat);
  
  newTransform.rotation = svd.V() * svd.U().transpose();
  // Cumulatively multiply rotation values
  transform.rotation = newTransform.rotation * transform.rotation;
  // TODO: Calculate standardized distance between mean of shapes and registered shape?
  // Calculate scale: Sum up elements of diagonal matrix
  double trsqrt = 0;
  for(int j = 0; j<dim; j++)
    {
    trsqrt += svd.W()(j);
    }
  newTransform.scale = trsqrt * (normMean / normShape);
  
  if(newTransform.scale == 0)
    newTransform.scale = 1.0;
  
  // Cumulatively multiply scale values
  transform.scale *= newTransform.scale;
  
  // Calculate translation
  PointType mult1 = newTransform.scale * colMeanShape;
  
  PointType mult2;
  mult2.fill(0.0);
  
  for(int i = 0; i<dim; i++)
    {
    for(int j = 0; j<dim; j++)
      {
      mult2[i] += mult1[j] * newTransform.rotation(j,i);
      }
    }
  
  PointType sub = colMeanMean - mult2;
  newTransform.translation = sub;
  // Cumulatively add translation values
  transform.translation += newTransform.translation;
  // Transform the shape
  ShapeType outputShape = TransformShape((*leaveOutIt), newTransform);
  // Re-initialize variables
  colMeanShape.fill(0.0);
  colMeanMean.fill(0.0);
  
  shapeScaled.clear();
  meanScaled.clear();
  
  ssqShape.fill(0.0);
  ssqMean.fill(0.0);
  
  normMean = 0;
  normShape = 0;
  
  return outputShape;
}

PSMProcrustesFunction::ShapeType PSMProcrustesFunction
::TransformShape(ShapeType shape, SimilarityTransform3D & transform)
{
  int numPoints = shape.size();
  int dim = shape[0].size();
  ShapeIteratorType shapeIt;
  shapeIt = shape.begin();
  
  // Multiply by scale
  while(shapeIt != shape.end())
    {
    PointType & point = *shapeIt;
    (*shapeIt) = transform.scale * point;
    shapeIt++;
    }
  
  ShapeType transformedShape;
  for(int i = 0; i < numPoints; i++)
    transformedShape.push_back(vnl_vector_fixed<RealType, 3>(0.0, 0.0, 0.0));
  
  // Multiply by rotation
  for(int i = 0; i<numPoints; i++)
    {
    for(int j = 0; j<dim; j++)
      {
      for(int k = 0; k<dim; k++)
        {
        transformedShape[i][j] += shape[i][k] * transform.rotation[k][j];
        }
      }
    }
  
  shapeIt = transformedShape.begin();
  
  // Add translation
  while(shapeIt != transformedShape.end())
    {
    PointType & point = (*shapeIt);
    point += transform.translation;
    shapeIt++;
    }
  return transformedShape;
}

PSMProcrustesFunction::RealType PSMProcrustesFunction
::ComputeSumOfSquares(ShapeListType & shapes)
{
  ShapeListIteratorType shapeIt1, shapeIt2;
  ShapeIteratorType pointIt1, pointIt2;
  
  RealType sum = 0.0;
  
  for(shapeIt1 = shapes.begin(); shapeIt1 != shapes.end(); shapeIt1++)
    {
    for(shapeIt2 = shapes.begin(); shapeIt2 != shapes.end(); shapeIt2++)
      {
      ShapeType & shape1 = (*shapeIt1);
      ShapeType & shape2 = (*shapeIt2);
      
      pointIt1 = shape1.begin();
      pointIt2 = shape2.begin();
      while(pointIt1 != shape1.end() && pointIt2 != shape2.end())
        {
        sum += ((*pointIt1) - (*pointIt2)).squared_magnitude();
        pointIt1++;
        pointIt2++;
        }
      }
    }
  return sum / static_cast<RealType>(shapes.size() * shapes[0].size());
}

bool PSMProcrustesFunction
::CheckDegenerateCase(PointType ssqShape, PointType ssqMean,
                      PointType colMeanShape, PointType colMeanMean, int rows)
{
  // TODO: Calculate standardized distance between mean of shapes and
  // registered shape?
  PointType valueShape, valueMean;
  int dim = ssqShape.size();
  for(int i = 0; i<dim; i++)
    {
    valueShape[i] = 2.22e-16 * rows * colMeanShape[i];
    valueShape[i] = valueShape[i] * valueShape[i];
    
    valueMean[i] = 2.22e-16 * rows * colMeanMean[i];
    valueMean[i] = valueMean[i] * valueMean[i];
    }
  
  // Check if any element in ssqShape and ssqMean is less than any element in
  // valueShape and valueMean resp.
  for(int j = 0; j<dim; j++)
    {
    if(ssqShape[j] <= valueShape.min_value() && ssqMean[j] <= valueMean.min_value())
      return false;
    }
  
  return true;
}	

void PSMProcrustesFunction
::LeaveOneOutMean(ShapeType & mean, ShapeListType & shapeList, ShapeListIteratorType & leaveOutIt)
{
  ShapeListIteratorType shapeListIt;
  ShapeIteratorType shapeIt, meanIt;
  
  int i, numPoints = shapeList[0].size();
  
  mean.clear();
  mean.reserve(numPoints);
  for(i = 0; i < numPoints; i++)
    {
    mean.push_back(vnl_vector_fixed<RealType, 3>(0.0, 0.0, 0.0));
    }
    
  for(shapeListIt = shapeList.begin(); shapeListIt != shapeList.end(); shapeListIt++)
    {
    if(shapeListIt != leaveOutIt)
      {
      ShapeType & shape = (*shapeListIt);
      shapeIt = shape.begin();
      meanIt = mean.begin();
      while(shapeIt != shape.end())
        {
        (*meanIt) += (*shapeIt);
        
        shapeIt++;
        meanIt++;
        }
      }
    }
  
  for(meanIt = mean.begin(); meanIt != mean.end(); meanIt++)
    {
    (*meanIt) /= static_cast<RealType>(shapeList.size() - 1);
    }
}

} // end namespace
