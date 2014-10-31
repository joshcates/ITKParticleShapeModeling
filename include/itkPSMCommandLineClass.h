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

#ifndef ____itkPSMCommandLineClass__
#define ____itkPSMCommandLineClass__

#include <iostream>
#include <string>
#include <sstream>
#include "itkImage.h"
#include "itkPSMProcrustesRegistration.h"
#include "itkImageFileReader.h"
#include "itkPSMEntropyModelFilter.h"
#include "itkPSMProject.h"
#include "itkPSMProjectReader.h"
#include "itkPSMParticleSystem.h"
#include "itkCommand.h"

namespace itk
{
/** \class PSMCommandLineClass
 *  \brief This class provides a command line tool to run the Particle Shape
 *   Modeling. It runs the optimization as well as Procrustes
 *   Registration if the option is in the project parameter file. For now, 
 *   the tool only uses the PSMEntropyModelFilter for the optimization.
 *
 *
 *  \ingroup PSM
 */

template <unsigned int VDimension>
class ITK_EXPORT PSMCommandLineClass : public DataObject
{
  public:
    /** Standard class typedefs. */
    typedef PSMCommandLineClass  Self;
    typedef DataObject Superclass;
    typedef SmartPointer<Self>   Pointer;
    typedef SmartPointer<const Self>  ConstPointer;
  
    /** Method for creation through the object factory. */
    itkNewMacro(Self);
  
    /** Run-time type information (and related methods). */
    itkTypeMacro(PSMCommandLineClass, DataObject);
  
    /** Input distance transforms image typedef */
    typedef typename itk::Image<float, VDimension> ImageType;
    /** PSM model optimization filter typedef */
    typedef PSMEntropyModelFilter<typename PSMCommandLineClass::ImageType> EntropyModelFilterType;
    /** Procrustes Registration typedef */
    typedef PSMProcrustesRegistration<VDimension> ProcrustesRegistrationType;
    /** Project Reader typedef */
    typedef PSMProjectReader ProjectReaderType;
    /** Project typedef */
    typedef PSMProject ProjectType;
  
    /** Read the distance transforms that are provided as inputs to the 
     *  optimzation filter. */
    void ReadInputs(std::string input_path_prefix);
    /** Read the input parameters that set various optimization attribute values */
    void ReadInputParameters();
    /** Write out the optimized point sets to user specified files */
    void WriteOutputs(std::string output_path);
    /** Run the steps of the optimization process */
    void Run(const char *fname, std::string input_path_prefix, std::string output_path);
  
    /** Constructor and destructor */
    PSMCommandLineClass();
    virtual ~PSMCommandLineClass() {};
  
  protected:
    /** Set the file name of the project parameter file */
    void SetProjectParameterFileName(const char *fname);
    /** Set default optimization scale parameters */
    void SetDefaultScales(unsigned int num);
    /** Callback to run Procrustes Registration on the shapes at the interval
    *  specified in the project parameter file or by default. */
    void IterateCallback(itk::Object *, const itk::EventObject &);
  
  private:
    PSMCommandLineClass(const Self&); //purposely not implemented
    void operator=(const Self&);      //purposely not implemented
    /** The project parameter file that contains information about input and output 
     *  file names and the optimization attribute values */
    const char *m_ProjectParameterFile;
    /** The PSM model optimization filter */
    typename EntropyModelFilterType::Pointer m_Filter;
    /** The Procrustes Registration filter that will register the point sets */
    typename ProcrustesRegistrationType::Pointer m_ProcrustesRegistration;
    /** The reader which will read in the project parameters */
    typename ProjectReaderType::Pointer m_XmlReader;
    /** The variable which will store the project parameters */
    typename ProjectType::Pointer m_Project;
    /** Counter to keep track of when to run Procrustes Registration */
    int m_ProcrustesCounter;
    /** This variable calls a pointer to a member function, in this case, the
     *  IterateCallback function which will run Procrustes at specified intervals
     *  during the optimization */
    typename itk::MemberCommand<PSMCommandLineClass>::Pointer m_IterateCmd;
};
} // end namespace

#endif 
