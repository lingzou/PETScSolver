#include "InputParser.h"

InputParser::InputParser(const char* input_file_name)
{
  ifile_name = std::string(input_file_name);
  ptr_ifile = new GetPot(input_file_name);
  buildGlobalParamList();
}

InputParser::~InputParser() { delete ptr_ifile; delete global_ParamList; }

void InputParser::buildGlobalParamList()
{
  global_ParamList = new InputParameterList("GlobalParamList", *ptr_ifile);
  // ptr_ifile->set_prefix("Global/");

  global_ParamList->readRequiredInputParameter<std::string>("problem");
  global_ParamList->readRequiredInputParameter<int>("n_steps");
  global_ParamList->readRequiredInputParameter<double>("dt");
  global_ParamList->readRequiredInputParameter<TimeScheme>("ts");
  global_ParamList->readOptionalInputParameter<int>("output_interval", 1);
  global_ParamList->readOptionalInputParameter<bool>("text_output", false);
  global_ParamList->AddParameter<std::string>("input_file_name", UTILS::trim_file_name(ifile_name));
}
