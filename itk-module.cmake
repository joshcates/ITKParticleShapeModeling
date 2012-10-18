set(DOCUMENTATION "This module contains code for computing shape models of surfaces and segmentations using particle systems and optimized correspondence calculations.")

itk_module(ITKParticleShapeModeling
  DEPENDS
    ITKCommon
    ITKImageFeature
    ITKImageFunction
    ITKImageFilterBase
    ITKImageGradient
    ITKConnectedComponents
    ITKRegionGrowing
    ITKAntiAlias
    ITKLevelSets
    ITKIONRRD
    ITKIOMeta
    ITKIOGDCM
    ITKIOXML
    ITKTestKernel #to handle IO in src
    ITKIOImageBase
    ITKIOXML
  TEST_DEPENDS
    ITKTestKernel
    ITKIOXML
  DESCRIPTION
    "${DOCUMENTATION}"
)

