#pragma once

//#include <GetPot>
#include "ParameterList.h"

class InputParser
{
public:
  InputParser(const char* input_file_name);
  virtual ~InputParser();

  InputParameterList& getGlobalParamList() { return *_global_ParamList; }
  std::map<std::string, InputParameterList *>& getProblemSystemParamList() { return _problemParamList_map; }
  std::map<std::string, InputParameterList *>& getFluidParamList() { return _fluidParamList_map; }

  // helper function
  void print() { _ifile->print(); }
  void checkUnusedVariables() { _ifile->checkUnvisitedVariables(); }

protected:
  void buildGlobalParamList();
  void prepareProblemParamList();

  std::string ifile_name;
  GetPot * _ifile;

  InputParameterList * _global_ParamList;
  std::map<std::string, InputParameterList *> _problemParamList_map;
  std::map<std::string, InputParameterList *> _fluidParamList_map;
};
