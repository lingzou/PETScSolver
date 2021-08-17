#ifndef PARAMETER_LIST_H
#define PARAMETER_LIST_H

#include <map>
#include <iomanip>      // std::setw
#include <vector>

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
  virtual ~ParameterList()
  {
    for(it = parameter_map.begin(); it != parameter_map.end(); ++it)
      delete it->second;
  }

  template<typename T>
  void AddParameter(std::string para_name, T para_value)
  {
    if(parameter_map.find(para_name) == parameter_map.end())  // not already existing, good to add
      parameter_map[para_name] = new ParameterT<T>(para_value);
    else
      sysError("ERROR: parameter with name: '" + para_name + "', already exists in ParameterList: '" + _name + "'.");
  }

  template<typename T>
  void setParameterValue(std::string para_name, T para_value)
  {
    it = parameter_map.find(para_name);
    if(it == parameter_map.end())
      sysError("ERROR: parameter with name: '" + para_name + "', does not exist in ParameterList: '" + _name + "'.");
    else
    {
      ParameterT<T> * param = dynamic_cast<ParameterT<T> * >(it->second);
      if(param == NULL)  // Check if the type is right
        sysError("ERROR: parameter with name: '" + para_name + "', has a different type in ParameterList: '" + _name + "'.");
      else
        param->setValue(para_value);
    }
  }

  template<typename T>
  T getParameterValue(std::string para_name)
  {
    it = parameter_map.find(para_name);
    if(it == parameter_map.end())
      sysError("ERROR: parameter with name: '" + para_name + "', does not exist in ParameterList: '" + _name + "'.");
    else
    {
      ParameterT<T> * param = dynamic_cast<ParameterT<T> * >(it->second);
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
  std::map<std::string, Parameter*>::iterator it;
};
#endif /*PARAMETER_LIST_H*/
