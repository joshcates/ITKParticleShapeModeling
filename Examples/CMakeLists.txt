OPTION(PSM_BUILD_EXAMPLES "Build Particle Shape Modeling module Example applications." OFF)

IF (PSM_BUILD_EXAMPLES)
ADD_EXECUTABLE( psmPreprocessor  psmPreprocessor.cxx)
TARGET_LINK_LIBRARIES(psmPreprocessor ${ITKParticleShapeModeling_LIBRARIES} ${ITKIOImageBase_LIBRARIES})
ENDIF (PSM_BUILD_EXAMPLES)
