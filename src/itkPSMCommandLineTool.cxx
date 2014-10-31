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
#include "itkExceptionObject.h"

int main( int argc, char *argv[] )
{
  int dim = 0;
  std::string output_path = "";
  std::string input_path_prefix = "";
  std::string errstring = "";
  // Check for proper arguments
  if (argc < 3)
  {
    std::cout << "Wrong number of arguments. \nUse: "
    << "ParticleShapeModeling_CLI parameter_file shape_dimensions [output_path] [input_path]\n"
    << "See itk::PSMParameterFileReader for documentation on the parameter file format.\n"
    << "Enter number of dimensions of the shapes (2 or 3) following parameter file name.\n"
    << "Note that input_path will be prefixed to any file names and paths in the xml parameter file.\n"
    << std::endl;
    return EXIT_FAILURE;
  }
  
  if (argc > 3)
  {
    output_path = std::string(argv[3]);
  }
  
  if (argc > 4)
  {
    input_path_prefix = std::string(argv[4]);
  }
  
  try
  {
    dim = atoi(argv[2]);
    // This function is called to fix an ITK runtime error where image format is not recognized.
    RegisterRequiredFactories();
    if (dim == 2)
    {
      itk::PSMCommandLineClass<2>::Pointer psmClass = itk::PSMCommandLineClass<2>::New();
      psmClass->Run( argv[1], input_path_prefix, output_path );
    }
    else if (dim == 3)
    {
      itk::PSMCommandLineClass<3>::Pointer psmClass = itk::PSMCommandLineClass<3>::New();
      psmClass->Run( argv[1], input_path_prefix, output_path );
    }
    else
    {
      std::cerr << "Please specify the input dimension as 2 or 3" << std::endl;
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