#include "utils.h"
#include "ParameterList.h"

InputParameterList::InputParameterList(std::string name, GetPot &file) :
  ParameterList(name), ifile(file), prefix("")
{
}

template <>
std::string InputParameterList::getValueFromInput<std::string>(std::string para_name)
{
  ifile.set_prefix(prefix.c_str());
  if (!hasInput(para_name))   sysError("Inputfile does not have parameter: " + para_name);
  return ifile(para_name.c_str(), "string");
}

template <>
int InputParameterList::getValueFromInput<int>(std::string para_name)
{
  return std::stoi(getValueFromInput<std::string>(para_name));
}

template <>
double InputParameterList::getValueFromInput<double>(std::string para_name)
{
  return std::stod(getValueFromInput<std::string>(para_name));
}

template <>
bool InputParameterList::getValueFromInput<bool>(std::string para_name)
{
  return UTILS::StringToBool(getValueFromInput<std::string>(para_name));
}

template <>
TimeScheme InputParameterList::getValueFromInput<TimeScheme>(std::string para_name)
{
  return UTILS::StringToEnum(getValueFromInput<std::string>(para_name));
}

template <>
void InputParameterList::addOptionalParamFromInput<bool>(std::string para_name, bool default_val)
{
  ifile.set_prefix(prefix.c_str());
  if (hasInput(para_name))
  {
    std::string str = ifile(para_name.c_str(), "string");
    AddParameter<bool>(para_name, UTILS::StringToBool(str));
  }
  else
    AddParameter<bool>(para_name, default_val);
}
