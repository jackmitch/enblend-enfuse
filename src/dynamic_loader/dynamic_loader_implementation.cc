#include "dynamic_loader_implementation.h"


DynamicLoaderImplementation::DynamicLoaderImplementation(const std::string& a_library_name) :
    name_(a_library_name)
{}


const std::string&
DynamicLoaderImplementation::library_name() const
{
    return name_;
}
