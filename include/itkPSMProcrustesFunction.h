/*=========================================================================
 Program:   ShapeWorks: Particle-based Shape Correspondence & Visualization
 Module:    $RCSfile: itkPSMProcrustesFunction.h,v $
 Date:      $Date: 2013/21/08 15:06:15 $
 Version:   $Revision: 1.1.1.1 $
 Author:    $Author: bengali $
 =========================================================================*/

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
    struct SimilarityTransform3D
    {
        vnl_matrix_fixed<double, 3, 3> rotation;
        double scale;
        vnl_vector_fixed<double, 3> translation;
    };
    
    class PSMProcrustesFunction : public itk::Object
    {
    public:
        typedef double                                  RealType;
        typedef vnl_vector_fixed<double, 3>             PointType;
        typedef std::vector<PointType>                  ShapeType;
        typedef ShapeType::iterator                     ShapeIteratorType;
        typedef std::vector<ShapeType>                  ShapeListType;
        typedef ShapeListType::iterator                 ShapeListIteratorType;
        typedef std::vector<SimilarityTransform3D>      SimilarityTransformListType;
        typedef SimilarityTransformListType::iterator   SimilarityTransformListIteratorType;
        /** Standard class typedefs. */
        typedef PSMProcrustesFunction                   Self;
        typedef Object                                  Superclass;
        typedef SmartPointer< Self >                    Pointer;
        typedef SmartPointer< const Self >              ConstPointer;
        
        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        /** Run-time type information (and related methods). */
        itkTypeMacro(PSMProcrustesFunction, Object);
        
        // Align a list of shapes using Generalized Procrustes Analysis
        void RunGeneralizedProcrustes(SimilarityTransformListType & transform, ShapeListType & shapes);
		
        RealType ComputeSumOfSquares(ShapeListType & shapes);
        // Helper function to transform a shape by a similarity transform
        ShapeType TransformShape(ShapeType shape, SimilarityTransform3D & transform);
        
    private:        
        // Check if shapes are the same
        bool CheckDegenerateCase(PointType ssqShape, PointType ssqMean, PointType muShape, PointType muMean, int rows);
        // Compute mean of all shapes except one
        void LeaveOneOutMean(ShapeType & mean, ShapeListType & shapeList, ShapeListIteratorType & leaveOutIt);
        // Align two shapes (translation, rotation & scale) using ordinary Procrustes Analysis
        ShapeType RunProcrustes(SimilarityTransform3D & transform, ShapeType mean, ShapeListIteratorType & leaveOutIt);
        // TODO: Template the class to allow for N-D
    };
    
} // end namespace itk
#endif
