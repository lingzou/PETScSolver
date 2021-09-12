#include "utils.h"
#include "ParameterList.h"

InputParameterList::InputParameterList(std::string name, GetPot &file) :
  ParameterList(name), ifile(file), prefix("")
{
}

template <>
void InputParameterList::readRequiredInputParameter<int>(std::string para_name)
{
  ifile.set_prefix(prefix.c_str());
  if (hasInput(para_name))
    AddParameter<int>(para_name, ifile(para_name.c_str(), int(1)));
  else
    sysError("Inputfile does not have parameter: " + para_name);
}

template <>
void InputParameterList::readRequiredInputParameter<double>(std::string para_name)
{
  ifile.set_prefix(prefix.c_str());
  if (hasInput(para_name))
    AddParameter<double>(para_name, ifile(para_name.c_str(), double(1.0)));
  else
    sysError("Inputfile does not have parameter: " + para_name);
}

template <>
void InputParameterList::readRequiredInputParameter<std::string>(std::string para_name)
{
  ifile.set_prefix(prefix.c_str());
  if (hasInput(para_name))
    AddParameter<std::string>(para_name, ifile(para_name.c_str(), "string"));
  else
    sysError("Inputfile does not have parameter: " + para_name);
}

template <>
void InputParameterList::readRequiredInputParameter<TimeScheme>(std::string para_name)
{
  ifile.set_prefix(prefix.c_str());
  if (hasInput(para_name))
  {
    std::string ts_str = ifile(para_name.c_str(), "string");
    AddParameter<TimeScheme>(para_name, UTILS::StringToEnum(ts_str));
  }
  else
    sysError("Inputfile does not have parameter: " + para_name);
}

template <>
void InputParameterList::readRequiredInputParameter<bool>(std::string para_name)
{
  ifile.set_prefix(prefix.c_str());
  if (hasInput(para_name))
  {
    std::string str = ifile(para_name.c_str(), "string");
    AddParameter<bool>(para_name, UTILS::StringToBool(str));
  }
  else
    sysError("Inputfile does not have parameter: " + para_name);
}

template <>
void InputParameterList::readOptionalInputParameter<bool>(std::string para_name, bool default_val)
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
