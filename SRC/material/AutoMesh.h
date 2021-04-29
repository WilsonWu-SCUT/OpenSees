#pragma once



#ifdef AUTOMESH_MODULE
#define DLL_AUTOMESH_API __declspec(dllexport)
#else
#define DLL_AUTOMESH_API __declspec(dllimport)
#pragma comment(lib, "AutoMesh.lib")
#endif



#include "AutoMesh//AutoMeshEnum.h"

#include "AutoMesh//MeshNode.h"
#include "AutoMesh//MeshElement.h"
#include "AutoMesh//MeshRegion.h"

#include "AutoMesh//FRAMSection.h"

#include "AutoMesh//MeshCpp.h"



