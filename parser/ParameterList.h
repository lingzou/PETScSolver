#ifndef PARAMETER_LIST_H
#define PARAMETER_LIST_H

#include <map>
#include <iomanip>      // std::setw
#include <vector>
#include <GetPot>

#include "utils.h"

// Base class for Parameter
class Parameter
{
public:
  Parameter() {}
  virtual ~Parameter() {}
};

// Using template, this is the real Parameter class to deal with different data type,
template<typename T>
class ParameterT : public Parameter
{
public:
  ParameterT() {}
  ParameterT(T value) : _value(value) {}
  virtual ~ParameterT() {}

  T getValue() const { return _value; }
  void setValue(T value) { _value = value; }

protected:
  T _value;
};

// A list of key-value Parameter
class ParameterList
{
public:
  ParameterList(std::string name) : _name(name) {}
  virtual ~ParameterList()  { for(auto & it : parameter_map)  delete it.second; }

  std::string name() { return _name; }

  template<typename T>
  void AddParameter(std::string para_name, T para_value)
  {
    if(parameter_map.find(para_name) != parameter_map.end())      sysError("ERROR: parameter with name: '" + para_name + "', already exists in ParameterList: '" + _name + "'.");
    parameter_map[para_name] = new ParameterT<T>(para_value);
  }

  template<typename T>
  void setParameterValue(std::string para_name, T para_value)
  {
    if(parameter_map.find(para_name) == parameter_map.end())      sysError("ERROR: parameter with name: '" + para_name + "', does not exist in ParameterList: '" + _name + "'.");
    else
    {
      ParameterT<T> * param = dynamic_cast<ParameterT<T> * >(parameter_map[para_name]);
      if(param == NULL)  // Check if the type is right
        sysError("ERROR: parameter with name: '" + para_name + "', has a different type in ParameterList: '" + _name + "'.");
      else
        param->setValue(para_value);
    }
  }

  template<typename T>
  T getParameterValue(std::string para_name)
  {
    if(parameter_map.find(para_name) == parameter_map.end())      sysError("ERROR: parameter with name: '" + para_name + "', does not exist in ParameterList: '" + _name + "'.");
    else
    {
      ParameterT<T> * param = dynamic_cast<ParameterT<T> * >(parameter_map[para_name]);
      if(param == NULL)  // Check if the type is right
        sysError("ERROR: parameter with name: '" + para_name + "', has a different type in ParameterList: '" + _name + "'.");
      else
        return param->getValue();
    }
    // Ok, code should not reach here but compiler complains for not giving return value.
    // and we don't know what T() really does
    return T();
  }

  bool hasParameter(std::string name)
  { return (!(parameter_map.find(name) == parameter_map.end())); }

protected:
  std::string _name;
  std::map<std::string, Parameter*> parameter_map;
};

class InputParameterList : public ParameterList
{
public:
  InputParameterList(std::string name, GetPot &file);
  virtual ~InputParameterList() {}

  void setGetPotPrefix(std::string section) { prefix = section; }
  bool hasInput(std::string para_name) { return ifile.hasVariable(para_name); }

  template <typename T>
  T getValueFromInput(std::string para_name);

  template <typename T>
  void addRequiredParamFromInput(std::string para_name)
  { AddParameter<T>(para_name, getValueFromInput<T>(para_name)); }

  template <typename T>
  void addOptionalParamFromInput(std::string para_name, T default_val)
  { ifile.set_prefix(prefix.c_str());
    AddParameter<T>(para_name, ifile(para_name.c_str(), default_val)); }

protected:
  GetPot &ifile;
  std::string prefix;
};

#endif /*PARAMETER_LIST_H*/
