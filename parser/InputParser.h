#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

//#include <GetPot>
#include "ParameterList.h"

class InputParser
{
public:
  InputParser(const char* input_file_name);
  virtual ~InputParser();

  InputParameterList& getGlobalParamList() { return *global_ParamList; }
  std::map<std::string, InputParameterList *>& getProblemSystemParamList() { return problemParamList_map; }

  // helper function
  void print() { ifile->print(); }

protected:
  void buildGlobalParamList();
  void prepareProblemParamList();

  std::string ifile_name;
  GetPot * ifile;

  InputParameterList * global_ParamList;
  std::map<std::string, InputParameterList *> problemParamList_map;
};

#endif /*INPUT_PARSER_H*/
