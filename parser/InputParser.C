#include "InputParser.h"

InputParser::InputParser(const char* input_file_name)
{
  ifile_name = std::string(input_file_name);
  ifile = new GetPot(input_file_name);
  buildGlobalParamList();
  prepareProblemParamList();
}

InputParser::~InputParser()
{
  delete ifile;
  delete global_ParamList;

  std::map<std::string, InputParameterList *>::iterator it;
  for(it = problemParamList_map.begin(); it != problemParamList_map.end(); ++it)
    delete it->second;
}

void
InputParser::buildGlobalParamList()
{
  global_ParamList = new InputParameterList("GlobalParamList", *ifile);
  global_ParamList->setGetPotPrefix("Global/");

  global_ParamList->readRequiredInputParameter<int>("n_steps");
  global_ParamList->readRequiredInputParameter<double>("dt");
  global_ParamList->readRequiredInputParameter<TimeScheme>("ts");
  global_ParamList->readOptionalInputParameter<int>("output_interval", 1);
  global_ParamList->readOptionalInputParameter<bool>("text_output", false);
  global_ParamList->AddParameter<std::string>("input_file_name", UTILS::trim_file_name(ifile_name));
}

void
InputParser::prepareProblemParamList()
{
  std::vector<std::string> input_sections = ifile->get_section_names();

  for(unsigned int i = 0; i < input_sections.size(); i++)
  {
    std::string name = input_sections[i];

    if((name.compare(0, 7, "System/") == 0) && (name.size() > 7) ) // Avoid processing the empty 'System/' section
    {
      name.erase(0, 7);   // remove prefix "System/"
      name.erase(name.size()-1, 1);  // remove "/" at the end
      problemParamList_map[name] = new InputParameterList(name, *ifile);
      problemParamList_map[name]->setGetPotPrefix(input_sections[i]);
      problemParamList_map[name]->readRequiredInputParameter<std::string>("type");
    }
  }
}
