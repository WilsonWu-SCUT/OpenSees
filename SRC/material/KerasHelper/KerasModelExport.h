#pragma once

#ifdef KERAS_MODULE
#define DLL_KERAS_API __declspec(dllexport)
#else
#define DLL_KERAS_API __declspec(dllimport)
#pragma comment(lib, "KerasHelperCPP.lib")
#endif

#include <memory>
#include <string>
#include <vector>

#include "model.h"
#include "KerasModel.h"

using namespace keras2cpp;