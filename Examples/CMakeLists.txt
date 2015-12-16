OPTION(PSM_BUILD_EXAMPLES "Build Particle Shape Modeling module applications." OFF)

IF (PSM_BUILD_EXAMPLES)
#ADD_EXECUTABLE( psmPreprocessor  psmPreprocessor.cxx)
#TARGET_LINK_LIBRARIES(psmPreprocessor ${ITKParticleShapeModeling_LIBRARIES} ${ITKIOImageBase_LIBRARIES})

ADD_EXECUTABLE(psmCommandLineTool psmCommandLineTool.cxx)
TARGET_LINK_LIBRARIES(psmCommandLineTool ITKIOMeta ITKIOGDCM ITKIONRRD ITKParticleShapeModeling ITKIOXML  ${ITKIOImageBase_LIBRARIES} ${ITKstatistics_LIBRARIES} ${ITKTransform_LIBRARIES})

ENDIF (PSM_BUILD_EXAMPLES)
