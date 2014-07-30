set(ITKParticleShapeModeling_SRC
itkPSMCommandLineClass.cxx
itkPSMDOMNodeXMLReader.cxx
itkPSMProject.cxx
itkPSMProjectReader.cxx
itkPSMProjectWriter.cxx
itkPSMParticleSystem.cxx
itkPSMProcrustesFunction.cxx
itkPSMProcrustesRegistration.cxx
)

add_library(ITKParticleShapeModeling ${ITKParticleShapeModeling_SRC})
target_link_libraries(ITKParticleShapeModeling  ${ITKstatistics_LIBRARIES} ${ITKTransform_LIBRARIES})
itk_module_target(ITKParticleShapeModeling)

add_executable(PSMCLI itkPSMCommandLineTool.cxx)
target_link_libraries(PSMCLI ${ITKParticleShapeModeling_LIBRARIES} ITKIOXML ${ITKIOImageBase_LIBRARIES} ITKMetaIO ${ITKstatistics_LIBRARIES} ${ITKTransform_LIBRARIES})
