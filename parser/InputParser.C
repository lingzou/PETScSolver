#include "InputParser.h"

InputParser::InputParser(const char* input_file_name)
{
  ifile_name = std::string(input_file_name);
  _ifile = new GetPot(input_file_name);
  buildGlobalParamList();
  prepareProblemParamList();
}

InputParser::~InputParser()
{
  delete _ifile;
  delete _global_ParamList;
  for (auto & it : _problemParamList_map)    delete it.second;
  for (auto & it : _fluidParamList_map)      delete it.second;
}

void
InputParser::buildGlobalParamList()
{
  _global_ParamList = new InputParameterList("GlobalParamList", *_ifile);
  _global_ParamList->setGetPotPrefix("Global/");

  _global_ParamList->readRequiredInputParameter<int>("n_steps");
  _global_ParamList->readRequiredInputParameter<double>("dt");
  _global_ParamList->readRequiredInputParameter<TimeScheme>("ts");
  // solver option: 0) Newton + hand-coded jacobian; 1) Matrix-free + hand-coded jacobian; 2) MF + fd coloring J; 3) MF + fd no coloring J
  _global_ParamList->readOptionalInputParameter<int>("solver_option", 2);
  _global_ParamList->readOptionalInputParameter<double>("linear_rtol", 1.e-3);
  _global_ParamList->readOptionalInputParameter<int>("linear_max_its", 30);
  _global_ParamList->readOptionalInputParameter<int>("output_interval", 1);
  _global_ParamList->readOptionalInputParameter<bool>("text_output", false);
  _global_ParamList->AddParameter<std::string>("input_file_name", UTILS::trim_file_name(ifile_name));
}

void
InputParser::prepareProblemParamList()
{
  std::vector<std::string> section_names = _ifile->get_section_names();

  for(unsigned i = 0; i < section_names.size(); i++)
  {
    std::string name = section_names[i];

    if((name.compare(0, 7, "System/") == 0) && (name.size() > 7) ) // Avoid processing the empty 'System/' section
    {
      name.erase(0, 7);   // remove prefix "System/"
      name.erase(name.size()-1, 1);  // remove "/" at the end
      _problemParamList_map[name] = new InputParameterList(name, *_ifile);
      _problemParamList_map[name]->setGetPotPrefix(section_names[i]);
      _problemParamList_map[name]->readRequiredInputParameter<std::string>("type");
    }
    else if ((name.compare(0, 7, "Fluids/") == 0) && (name.size() > 7) ) // Avoid processing the empty 'System/' section
    {
      name.erase(0, 7);   // remove prefix "Fluids/"
      name.erase(name.size()-1, 1);  // remove "/" at the end
      _fluidParamList_map[name] = new InputParameterList(name, *_ifile);
      _fluidParamList_map[name]->setGetPotPrefix(section_names[i]);
      _fluidParamList_map[name]->readRequiredInputParameter<std::string>("type");
    }
  }
}
