#pragma once

#ifdef SECTIONANALYSIS_MODULE
#define DLL_SECTIONANALYSIS_API __declspec(dllexport)
#else
#define DLL_SECTIONANALYSIS_API __declspec(dllimport)
#pragma comment(lib, "SectionAnalysis.lib")
#endif

#include "SectionAnalysis\\MomentAxialLoad.h"
#include "SectionAnalysis\\MomentAxialLoadSet.h"

#include "SectionAnalysis\\Section.h"
