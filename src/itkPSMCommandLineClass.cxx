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

#include "itkPSMCommandLineClass.h"
#include <fstream>
#include <vector>

namespace itk
{
template class PSMCommandLineClass<3>;

template <unsigned int VDimension>
PSMCommandLineClass<VDimension>
::PSMCommandLineClass()
{
  this->m_ProcrustesCounter = 0;
}
  
template <unsigned int VDimension>
void PSMCommandLineClass<VDimension>
::IterateCallback(itk::Object *caller , const itk::EventObject &)
{
  // Check if the Procrustes interval is set to a value other than 0
  if (this->m_ProcrustesRegistration->GetProcrustesInterval() != 0)
  {
    this->m_ProcrustesCounter++;
    // If the counter is greater than the interval value, run Procrustes registration
    if (this->m_ProcrustesCounter >= (int)this->m_ProcrustesRegistration->GetProcrustesInterval())
    {
      // Reset the counter
      this->m_ProcrustesCounter = 0;
      this->m_ProcrustesRegistration->RunRegistration();
      std::cout << "Run Procrustes Registration" << std::endl;
    }
  }
  PSMEntropyModelFilter<PSMCommandLineClass::ImageType> *o
  = static_cast<PSMEntropyModelFilter<PSMCommandLineClass::ImageType> *>(caller);
  // Print every 10 iterations
  if (o->GetNumberOfElapsedIterations() % 10 != 0) return;
  
  std::cout << "Iteration # " << o->GetNumberOfElapsedIterations() << std::endl;
  std::cout << " Eigenmode variances: ";
  for (unsigned int i = 0; i < o->GetShapePCAVariances().size(); i++)
  {
    std::cout << o->GetShapePCAVariances()[i] << " ";
  }
  std::cout << std::endl;
  std::cout << " Regularization = " << o->GetRegularizationConstant() << std::endl;
}

template <unsigned int VDimension>
void PSMCommandLineClass<VDimension>
::ReadInputs(std::string input_path_prefix)
{
  // Read the project parameter file
  this->m_XmlReader = PSMCommandLineClass::ProjectReaderType::New();
  this->m_XmlReader->SetFileName(this->m_ProjectParameterFile);
  this->m_XmlReader->Update();
  
  // Store the project parameters
  this->m_Project = PSMCommandLineClass::ProjectType::New();
  this->m_Project = m_XmlReader->GetOutput();
  
  // Read in the distance transforms
  const std::vector<std::string> &dt_files = this->m_Project->GetDistanceTransforms();
  this->m_Filter = PSMCommandLineClass::EntropyModelFilterType::New();
  std::cout << "Reading distance transforms ..." << std::endl;
  for (unsigned int i = 0; i < dt_files.size(); i++)
  {
    itk::ImageFileReader<PSMCommandLineClass::ImageType>::Pointer reader =
    itk::ImageFileReader<PSMCommandLineClass::ImageType>::New();
    reader->SetFileName(input_path_prefix + dt_files[i]);
    reader->Update();
    
    std::cout << "  " << dt_files[i] << std::endl;
    this->m_Filter->SetInput(i, reader->GetOutput());
  }
  std::cout << "Done!" << std::endl;
  // Load the model initialization.  It should be specified as a model with a name.
  const std::vector<std::string> &pt_files = this->m_Project->GetModel(std::string("initialization"));
  std::cout << "Reading the initial model correspondences ..." << std::endl;
  for (unsigned int i = 0; i < pt_files.size(); i++)
  {
    // Read the points for this file and add as a list
    std::vector<PSMCommandLineClass::EntropyModelFilterType::PointType> c;
    
    int counter = 0;
    // Open the ascii file.
    std::ifstream in( (input_path_prefix + pt_files[i]).c_str() );
    if ( !in )
    {
      std::cerr << "Could not open point file for input." << std::endl;
      break;
    }
    
    // Read all of the points, one point per line.
    while (in)
    {
      PSMCommandLineClass::EntropyModelFilterType::PointType pt;
      
      for (unsigned int d = 0; d < 3; d++)
      {
        in >> pt[d];
      }
      c.push_back(pt);
      counter++;
    }
    // This algorithm pushes the last point twice
    c.pop_back();
    in.close();
    
    this->m_Filter->SetInputCorrespondencePoints(i,c);
    
    std::cout << "  " << pt_files[i] << std::endl;
  }
  std::cout << "Done!" << std::endl;
  this->ReadInputParameters();
}
  
template <unsigned int VDimension>
void PSMCommandLineClass<VDimension>
::ReadInputParameters()
{
  this->m_ProcrustesRegistration = PSMCommandLineClass::ProcrustesRegistrationType::New();
  
  //  Provide some default parameters
  double regularization_initial   = 100.0f;
  double regularization_final     = 5.0f;
  double regularization_decayspan = 2000.0f;
  double tolerance                = 1.0e-8;
  unsigned int maximum_iterations = 200000;
  unsigned int procrustes_interval = 10;
  
  // Read the optimization parameters and set them in the filter
  std::cout << "Optimization parameters: " << std::endl;
  if ( this->m_Project->HasOptimizationAttribute("regularization_initial") )
  {
    regularization_initial = this->m_Project->GetOptimizationAttribute("regularization_initial");
    this->m_Filter->SetRegularizationInitial(regularization_initial);
    std::cout << "    regularization_initial = " << regularization_initial << std::endl;
  }
  if ( this->m_Project->HasOptimizationAttribute("regularization_final") )
  {
    regularization_final = this->m_Project->GetOptimizationAttribute("regularization_final");
    this->m_Filter->SetRegularizationFinal(regularization_final);
    std::cout << "      regularization_final = " << regularization_final << std::endl;
  }
  if ( this->m_Project->HasOptimizationAttribute("regularization_decayspan") )
  {
    regularization_decayspan = this->m_Project->GetOptimizationAttribute("regularization_decayspan");
    this->m_Filter->SetRegularizationDecaySpan(regularization_decayspan);
    std::cout << "  regularization_decayspan = " << regularization_decayspan << std::endl;
  }
  if ( this->m_Project->HasOptimizationAttribute("tolerance") )
  {
    tolerance = this->m_Project->GetOptimizationAttribute("tolerance");
    this->m_Filter->SetTolerance(tolerance);
    std::cout << "                 tolerance = " << tolerance << std::endl;
  }
  if ( this->m_Project->HasOptimizationAttribute("maximum_iterations") )
  {
    maximum_iterations = static_cast<unsigned int>(this->m_Project->GetOptimizationAttribute("maximum_iterations"));
    this->m_Filter->SetMaximumNumberOfIterations(maximum_iterations);
    std::cout << "        maximum_iterations = " << maximum_iterations << std::endl;
  }
  // Set Procrustes interval for ProcrustesRegistration
  if ( this->m_Project->HasOptimizationAttribute("procrustes_interval") )
  {
    procrustes_interval = static_cast<unsigned int>(this->m_Project->GetOptimizationAttribute("procrustes_interval"));
    this->m_ProcrustesRegistration->SetProcrustesInterval(procrustes_interval);
    std::cout << "        procrustes_interval = " << procrustes_interval << std::endl;
  }
  // Set ParticleSystem for ProcrustesRegistration
  this->m_ProcrustesRegistration->SetPSMParticleSystem(this->m_Filter->GetParticleSystem());
}

template <unsigned int VDimension>
void PSMCommandLineClass<VDimension>
::WriteOutputs(std::string output_path)
{
  const std::vector<std::string> &out_files = this->m_Project->GetModel(std::string("optimized"));
  if (out_files.size() != this->m_Filter->GetParticleSystem()->GetNumberOfDomains())
  {
    std::cerr << "Number of output files does not match the number of particle domains (inputs)." << std::endl;
  }  
  else
  {
    for (unsigned int d = 0; d < this->m_Filter->GetParticleSystem()->GetNumberOfDomains(); d++)
    {
      // Open the output file.
      std::string fname = output_path + out_files[d];
      std::ofstream out( fname.c_str() );
      if ( !out )
      {
        std::cerr << "Could not open point file for output: " << std::endl;
        break;
      }
      else
      {
        for (unsigned int j = 0; j < this->m_Filter->GetParticleSystem()->GetNumberOfParticles(d); j++)
        {
          for (unsigned int i = 0; i < 3; i++)
          {
            out <<  this->m_Filter->GetParticleSystem()->GetPosition(j,d)[i]  << " ";
          }
          out << std::endl;
        }
      }
    } // end for number of domains
  }
}
  
template <unsigned int VDimension>
void PSMCommandLineClass<VDimension>
::Run( const char *fname, std::string input_path_prefix, std::string output_path )
{
  // Set the name of the file containing project parameters
  this->SetProjectParameterFileName(fname);
  // Read the input point sets
  this->ReadInputs(input_path_prefix);
  // Create the callback function to run Procrustes
  m_IterateCmd = itk::MemberCommand<PSMCommandLineClass>::New();
  m_IterateCmd->SetCallbackFunction(this, &PSMCommandLineClass::IterateCallback);
  // Add the observer to the filter which will run depending on what the iteration count is
  this->m_Filter->AddObserver(itk::IterationEvent(), m_IterateCmd);
  // Run the optimization filter
  this->m_Filter->Update();
  // TODO: Is this error message really required?
  if (this->m_Filter->GetNumberOfElapsedIterations() >= this->m_Filter->GetMaximumNumberOfIterations())
  {
    std::cerr << "Optimization did not converge based on tolerance criteria." << std::endl;
  }
  // Write out the optimized point sets
  this->WriteOutputs(output_path);
}

template <unsigned int VDimension>
void PSMCommandLineClass<VDimension>
::SetProjectParameterFileName( const char *fname )
{
  this->m_ProjectParameterFile = fname;
}
// Explicit instantiation
template void PSMCommandLineClass<3>::Run( const char *fname, std::string input_path_prefix, std::string output_path );
} // end namespace itk

