//
//  XercesDOMParser.cxx
//  
//
//  Created by Salma Bengali on 9/16/14.
//
//

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>

#if defined(XERCES_NEW_IOSTREAMS)
#include <iostream>
#else
#include <iostream.h>
#endif

XERCES_CPP_NAMESPACE_USE

int main (int argc, char* args[]) {
  
  try {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch) {
    char* message = XMLString::transcode(toCatch.getMessage());
    cout << "Error during initialization! :\n"
    << message << "\n";
    XMLString::release(&message);
    return 1;
  }
  
  XercesDOMParser* parser = new XercesDOMParser();
  parser->setValidationScheme(XercesDOMParser::Val_Always);
  parser->setDoNamespaces(true);    // optional
  
  ErrorHandler* errHandler = (ErrorHandler*) new HandlerBase();
  parser->setErrorHandler(errHandler);
  
  char* xmlFile = "PSMCommandLineClassTest.xml";
  
  try {
    parser->parse(xmlFile);
  }
  catch (const XMLException& toCatch) {
    char* message = XMLString::transcode(toCatch.getMessage());
    cout << "Exception message is: \n"
    << message << "\n";
    XMLString::release(&message);
    return -1;
  }
  catch (const DOMException& toCatch) {
    char* message = XMLString::transcode(toCatch.msg);
    cout << "Exception message is: \n"
    << message << "\n";
    XMLString::release(&message);
    return -1;
  }
  catch (...) {
    cout << "Unexpected Exception \n" ;
    return -1;
  }
  
  delete parser;
  delete errHandler;
  return 0;
}

