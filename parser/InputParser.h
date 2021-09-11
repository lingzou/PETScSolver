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

  // helper function
  void print() { ptr_ifile->print(); }

protected:
  void buildGlobalParamList();

  std::string ifile_name;
  GetPot * ptr_ifile;

  InputParameterList * global_ParamList;
};

#endif /*INPUT_PARSER_H*/
