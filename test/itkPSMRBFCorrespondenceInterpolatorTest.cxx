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
#include <iostream>
#include "itkImage.h"
#include "itkPSMRBFCorrespondenceInterpolator.h"

/** This test exercises functionality of the base
    itkPSMRBFCorrespondenceInterpolator class */
int itkPSMRBFCorrespondenceInterpolatorTest(int argc, char* argv[] )
{
  bool passed = true;
  std::string errstring = "";
  std::string output_path = "";

  // Check for proper arguments
  if (argc < 2)
    {
      std::cout << "Wrong number of arguments. \nUse: " 
	<< "itkPSMRBFCorrespondenceInterpolatorTest PointFileA PointFileB PointListToInterpolate\n"
	<< std::endl;
      return EXIT_FAILURE;
    }
  else if (argc >2)
    {
      output_path = std::string(argv[2]);
    }

  typedef itk::Image<float, 3> ImageType;

  try
    {
     // Create the modeling filter and set up the optimization.
     itk::PSMRBFCorrespondenceInterpolator<3>::Pointer P 
       = itk::PSMRBFCorrespondenceInterpolator<3>::New();
 
     // Load the PointFileA
   //   const std::vector<std::string> &pt_files = project->GetModel(std::string("initialization"));
//      std::cout << "Reading the initial model correspondences ..." << std::endl;
//      for (unsigned int i = 0; i < pt_files.size(); i++)
//        {
//          // Read the points for this file and add as a list
//          std::vector<itk::PSMRBFCorrespondenceInterpolator<ImageType>::PointType> c;

//          int counter = 0;
//          // Open the ascii file.
//          std::ifstream in( pt_files[i].c_str() );
//          if ( !in )
//            {
//              errstring += "Could not open point file for input.";
//              passed = false;
//            }

//          // Read all of the points, one point per line.
//          while (in)
//            {
//              itk::PSMRBFCorrespondenceInterpolator<ImageType>::PointType pt;

//              for (unsigned int d = 0; d < 3; d++)
//                {
//                  in >> pt[d];
//                }
//              c.push_back(pt);
//              counter++;
//            }
//          // this algorithm pushes the last point twice
//          c.pop_back();
//          //  std::cout << "Read " << counter-1 << " points. " << std::endl;
//          in.close();

//          P->SetInputCorrespondencePoints(i,c);
                  
//          std::cout << "  " << pt_files[i] << std::endl;
//        }
//      std::cout << "Done!" << std::endl;
     

     if (passed = true)
       {

       }
     
    }
  catch(itk::ExceptionObject &e)
    {
      errstring = "ITK exception with description: " + std::string(e.GetDescription())
        + std::string("\n at location:") + std::string(e.GetLocation())
        + std::string("\n in file:") + std::string(e.GetFile());
      passed = false;
    }
  catch(...)
    {
      errstring = "Unknown exception thrown";
      passed = false;
    }

  if (passed)
    {
    std::cout << "All tests passed" << std::endl;
    return EXIT_SUCCESS;
    }
  else
    {
    std::cout << "Test failed with the following error:" << std::endl;
    std::cout << errstring << std::endl;
    return EXIT_FAILURE;
    }
}
