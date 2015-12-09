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

#include "itkPSMProject.h"

namespace itk
{

// Initialize reserved keywords for XML file tags
const std::string PSMProject::data_tag                = "data";
const std::string PSMProject::distance_transforms_tag = "distance_transforms";
const std::string PSMProject::model_tag               = "model";
const std::string PSMProject::name_tag                = "name";
const std::string PSMProject::optimization_tag        = "optimization";
const std::string PSMProject::number_of_scales_tag    = "number_of_scales";
const std::string PSMProject::psm_project_tag         = "psm_project";
const std::string PSMProject::preprocessing_tag       = "preprocessing";
const std::string PSMProject::procrustes_tag          = "procrustes_registration";
const std::string PSMProject::scale_tag               = "scale";
const std::string PSMProject::scale_number_tag        = "number";
const std::string PSMProject::variables_tag           = "variables";

void PSMProject::SetDOMNode(PSMDOMNode *dom)
{
  if (dom->GetName() == psm_project_tag)
    { 
      m_DOMNode = dom;
    }
  else
    {
      itkExceptionMacro( "DOM object does not appear to contain a valid PSM Project file." );
    }
}

unsigned int PSMProject
::GetNumberOfOptimizationScales() const
{
  DOMNode *opt = m_DOMNode->GetChild(optimization_tag);
  
  if (opt != 0)
    {
      if (opt->HasAttribute(number_of_scales_tag))
        {
          return static_cast<unsigned int>(atoi(opt->GetAttribute(number_of_scales_tag).c_str()));
        }
      else
        {
          itkExceptionMacro("Number of scales not specified");
        }
    }

  itkExceptionMacro("No " + optimization_tag + " was found.");
}

bool PSMProject
::HasOptimizationAttribute(const std::string& name, unsigned int s) const
{
  DOMNode *opt = m_DOMNode->GetChild(optimization_tag);
  if (opt != 0) // found the optimization element
    {
      unsigned int nscales = 1;
      
      // Are there scale elements?
      if (opt->HasAttribute(number_of_scales_tag))
        {
          nscales = static_cast<unsigned int>(atoi(opt->GetAttribute(number_of_scales_tag).c_str()));
          // If the number_of_scales_tag is the one being checked
          if(name == number_of_scales_tag) { return true; }
        }
      
      // Did the user ask for a scale that isn't specified?
      if (s+1 > nscales)  { return false; }
      
      // If only one scale, then look in element attributes for name
      if (nscales == 1) 
        {
          if (opt->HasAttribute(name)) { return true; }
        }
      
      // Finally, look for a scale element with name i
      DOMNode::ChildrenListType children;
      opt->GetAllChildren(children);

      for (unsigned i = 0; i < children.size(); i++)
        {
          if ( children[i]->GetName() == scale_tag 
               && 
               static_cast<unsigned int>(atoi(children[i]->GetAttribute(scale_number_tag).c_str())) == s )
            {
              if (children[i]->HasAttribute(name)) { return true; }
              else { return false; }
            }
        }
    }
  return false;
}

double PSMProject
::GetOptimizationAttribute(const std::string &name, unsigned int s) const
{
  DOMNode *opt = m_DOMNode->GetChild(optimization_tag);
  if (opt != 0) // found the optimization element
    {
      unsigned int nscales = 1;
      
      // Are there scale elements?
      if (opt->HasAttribute(number_of_scales_tag))
        {
          nscales = static_cast<unsigned int>(atoi(opt->GetAttribute(number_of_scales_tag).c_str()));
        }
      
      // Did the user ask for a scale that isn't specified?
      if (s+1 > nscales)
        { 
          itkExceptionMacro("Requested data for scale that was not specified.");
        }
      
      // If only one scale, then look in element attributes for name
      if (nscales == 1) 
        {
          if (opt->HasAttribute(name)) 
            { 
              return atof(opt->GetAttribute(name).c_str());
            }
        }
      
      // Finally, look for a scale element with name i
      DOMNode::ChildrenListType children;
      opt->GetAllChildren(children);

      for (unsigned i = 0; i < children.size(); i++)
        {
          if ( children[i]->GetName() == scale_tag 
               && 
               static_cast<unsigned int>(atoi(children[i]->GetAttribute(scale_number_tag).c_str())) == s )
            {
              if (children[i]->HasAttribute(name)) 
                { 
                  return atof(children[i]->GetAttribute(name).c_str()); 
                }
              else 
                { 
                  itkExceptionMacro("Found " + scale_tag + " element, but it does not have an attribute named " + name + ".");
                }
            }
        }
    }
  itkExceptionMacro("File has no " + optimization_tag + " element");
}

bool PSMProject::HasProcrustes() const
{
  DOMNode *proc = m_DOMNode->GetChild(procrustes_tag);
  if (proc != 0) // Found the procrustes element
  {
    if (proc->HasAttribute("value"))
    {
      if (atoi(proc->GetAttribute("value").c_str()) == 0)
        return false;
      else
        return true;
    }
  }
  return false;
}
  
bool PSMProject::HasVariables(const std::string &name) const
{
  const PSMDOMNode *data = this->GetDataNode();

  // Search the data tree for nodes called variables_tag
  DOMNode::ConstChildrenListType variables;
  data->GetChildren(variables_tag, variables);
  
  if (variables.size() == 0)
    {
      return false;
    }
  
  // Search the list of variables for one with the correct name
  for (unsigned int i = 0; i < variables.size(); i++)
    {
     if (variables[i]->HasAttribute(name_tag))
       {
         if (variables[i]->GetAttribute(name_tag) == name)
           {
             return true;
             //   return ((dynamic_cast<const PSMDOMNode *>(variables[i]))->GetText());
           }
       }
   }
  // itkExceptionMacro( "PSM Project file does not have any " + model_tag + " entries with name = " << name );
  return false;
}

std::vector<double> PSMProject::GetVariables(const std::string &name) const
{
  const std::vector<std::string> &vartext = this->GetVariablesText(name);
 
  std::vector<double> ans;
  for (unsigned int i = 0; i < vartext.size(); i++)
    {
      ans.push_back(atof(vartext[i].c_str()));      
    }
  return ans;
}

const std::vector<std::string > &PSMProject::GetVariablesText(const std::string &name) const
{
  const PSMDOMNode *data = this->GetDataNode();
  
  // Search the data tree for nodes called variables_tag
  DOMNode::ConstChildrenListType variables;
  data->GetChildren(variables_tag, variables);
  
  if (variables.size() == 0)
    {
      itkExceptionMacro( "PSM Project file does not have any " + variables_tag + " entries with name = " << name );
    }
  
  // Search the list of variables for one with the correct name
  for (unsigned int i = 0; i < variables.size(); i++)
    {
      if (variables[i]->HasAttribute(name_tag))
        {
          if (variables[i]->GetAttribute(name_tag) == name)
            {
              return ((dynamic_cast<const PSMDOMNode *>(variables[i]))->GetText());
            }
        }
    }
  itkExceptionMacro( "PSM Project file does not have any " + variables_tag + " entries with name = " << name );
}

const PSMDOMNode *PSMProject::GetDataNode() const
{
 // Look for the data section
  const PSMDOMNode *data 
    = static_cast<const PSMDOMNode *>(m_DOMNode->GetChild(data_tag));
  
  if (data == 0) // data is not found
    {
      itkExceptionMacro( "PSM Project file does not specify " + data_tag );
    }
  else
    {
      return data;
    }
}

const std::vector<std::string> &
PSMProject
::GetDistanceTransforms() const
{
  const PSMDOMNode *data = this->GetDataNode();

  const PSMDOMNode *dt = 0;
  dt = static_cast<const PSMDOMNode *>(data->GetChild(distance_transforms_tag));
  if (dt == 0)
    {
      itkExceptionMacro( "PSM Project file does not specify " + distance_transforms_tag );
    }
  else return dt->GetText();
}


const std::vector<std::string> &
PSMProject::GetModel(const std::string &name)
{
 const PSMDOMNode *data = this->GetDataNode();

 // Search the data tree for nodes called model_tag
 DOMNode::ConstChildrenListType models;
 data->GetChildren(model_tag, models);

 if (models.size() == 0)
   {
     itkExceptionMacro( "PSM Project file does not specify any elements called " + model_tag );
   }

 // Search the list of models for one with the correct name
 for (unsigned int i = 0; i < models.size(); i++)
   {
     if (models[i]->HasAttribute(name_tag))
       {
         if (models[i]->GetAttribute(name_tag) == name)
           {
             return ((dynamic_cast<const PSMDOMNode *>(models[i]))->GetText());
           }
       }
   }
 itkExceptionMacro( "PSM Project file does not have any " + model_tag + " entries with name = " << name );
}


bool PSMProject::HasModel(const std::string &name) const
{
  const PSMDOMNode *data = this->GetDataNode();
  
  // Search the data tree for nodes called model_tag
  DOMNode::ConstChildrenListType models;
  data->GetChildren(model_tag, models);
  
  if (models.size() == 0)
  {
    return false;
  }
  
  // Search the list of variables for one with the correct name
  for (unsigned int i = 0; i < models.size(); i++)
    {
      if (models[i]->HasAttribute(name_tag))
        {
          if (models[i]->GetAttribute(name_tag) == name)
            {
              return true;
            }
        }
    }
  return false;
}


void PSMProject::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
    //  this->PrintChildren(m_DOMNode,os,indent);

}

void PSMProject::Stream(std::ostream &os)
{
  Indent indent;

  this->StreamChildren(m_DOMNode,os,indent);
}

void PSMProject::StreamChildren(PSMDOMNode *node,std::ostream &os, Indent indent)
{
  // Print this node's name and attributes
  os << indent << "<" << node->GetName();  
  
  DOMNode::AttributesListType attr;
  node->GetAllAttributes(attr);
  
  for (DOMNode::AttributesListType::const_iterator it = attr.begin(); 
       it != attr.end(); it++)
    {
      os << " " << it->first << " = \"" << it->second << "\"";
    }
  
  
  os << ">" << std::endl;
  
  // Now print this node's text, line by line.
  for (unsigned int i = 0; i < node->GetText().size(); i++)
    {
      os << indent.GetNextIndent() << node->GetText()[i] << std::endl;
    }

  DOMNode::ChildrenListType children;
  node->GetAllChildren(children);

  // Now print all of this node's children.
  for (unsigned i = 0; i < children.size(); i++)
    {
      this->StreamChildren(static_cast<PSMDOMNode *>(children[i]),
                           os,indent.GetNextIndent());
    }

  // Print the closing tag for this node.
  os << indent << "<\\" << node->GetName() << ">" << std::endl;
}

} // end namespace itk
