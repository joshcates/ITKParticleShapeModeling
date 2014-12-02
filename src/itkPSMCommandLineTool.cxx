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

#include "itkIncludeRequiredIOFactories.h"
#include <iostream>
#include "itkPSMCommandLineClass.h"
#include "itkPSMCommandLineClass.cxx"
#include "itkPSMProjectReader.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkExceptionObject.h"

int main( int argc, char *argv[] )
{
  std::string output_path = "";
  std::string input_path_prefix = "";
  std::string errstring = "";
  // Check for proper arguments
  if (argc < 2)
  {
    std::cout << "Wrong number of arguments. \nUse: "
    << "ParticleShapeModeling_CLI parameter_file [output_path] [input_path]\n"
    << "See itk::PSMParameterFileReader for documentation on the parameter file format.\n"
    << "Note that input_path will be prefixed to any file names and paths in the xml parameter file.\n"
    << std::endl;
    return EXIT_FAILURE;
  }
  
  if (argc > 2)
  {
    output_path = std::string(argv[2]);
  }
  
  if (argc > 3)
  {
    input_path_prefix = std::string(argv[3]);
  }
  
  try
  {
    // This function is called to fix an ITK runtime error where image format is not recognized.
    RegisterRequiredFactories();
    
    // The dimensions of the input images need to be checked in order to
    // correctly initialize PSMCommandLineClass.
    itk::PSMProjectReader::Pointer xmlReader = itk::PSMProjectReader::New();
    xmlReader->SetFileName(argv[1]);
    xmlReader->Update();
    
    itk::PSMProject::Pointer project = xmlReader->GetOutput();
    // Load the distance transforms
    const std::vector<std::string> &dt_files = project->GetDistanceTransforms();
    std::string fname = input_path_prefix + dt_files[0];
    // Read the first distance transform to check the image dimensions
    std::cout << "Checking input image dimensions ..." << std::endl;
    typedef itk::ImageIOBase::IOComponentType ScalarPixelType;
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
                                       fname.c_str(), itk::ImageIOFactory::ReadMode);
    imageIO->SetFileName(fname);
    imageIO->ReadImageInformation();
    const size_t numOfDimensions =  imageIO->GetNumberOfDimensions();
    std::cout << "Number of dimensions: " << numOfDimensions << std::endl;
    
    if (numOfDimensions == 2)
    {
      itk::PSMCommandLineClass<2>::Pointer psmClass = itk::PSMCommandLineClass<2>::New();
      psmClass->Run( argv[1], input_path_prefix, output_path );
    }
    else if (numOfDimensions == 3)
    {
      itk::PSMCommandLineClass<3>::Pointer psmClass = itk::PSMCommandLineClass<3>::New();
      psmClass->Run( argv[1], input_path_prefix, output_path );
    }
  }
  
  catch(itk::ExceptionObject &e)
  {
    std::cerr << "ITK exception with description: " << e.GetDescription()
    << "\n at location:" << e.GetLocation()
    << "\n in file:" << e.GetFile() << std::endl;
    return EXIT_FAILURE;
  }
  
  catch(...)
  {
    errstring = "Unknown exception thrown";
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}