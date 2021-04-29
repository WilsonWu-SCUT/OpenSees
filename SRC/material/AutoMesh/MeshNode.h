#pragma once

namespace AutoMesh
{

	#define NODE_DIST_TOL 1E-6 

	class DLL_AUTOMESH_API MeshNode
	{
	public:
		MeshNode();
		MeshNode(const int& id, const double *rz);
		MeshNode(const double& x, const double& y);
		~MeshNode();

	public:
		int inline get_id() const { return this->ID; }
		double inline get_x() const { return this->X; }
		double inline get_y() const { return this->Y; }

	private:
		int ID;
		double X;
		double Y;
	};
}

