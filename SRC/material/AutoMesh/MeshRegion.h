#pragma once

class CAutomesh2DUI;

namespace AutoMesh
{
	//剖分节点
	class MeshNode;
	//剖分单元
	class MeshElement;

	//剖分子区域
	class MeshRegion
	{
	public:
		MeshRegion();
		MeshRegion(const int& id);
		MeshRegion(const int& id, std::shared_ptr<CAutomesh2DUI> amg);
		~MeshRegion();

	public:
		inline void set_mat_tag(const int& mat_type) { mat_type_ = mat_type; }
		inline int get_region_id(void) const { return region_id_; }
		inline int get_mat_type(void) const { return mat_type_; }

	public:
		inline int get_element_num(void) const
		{
			return this->element_map.size();
		}
		inline int get_node_num(void) const
		{
			return this->node_map.size();
		}
		inline double get_perimeter(void) const {
			return this->perimeter_;
		}
		inline double get_x(void) const {
			return this->centroid_sp_->get_x();
		}
		inline double get_y(void) const {
			return this->centroid_sp_->get_y();
		}

	public:
		std::shared_ptr<MeshElement> get_element_sp(const int& id);

	private:
		//区域编号
		int region_id_;
		//材料描述
		int mat_type_;
		//节点列表
		std::map<int, std::shared_ptr<MeshNode>> node_map;
		//单元列表
		std::map<int, std::shared_ptr<MeshElement>> element_map;
		//区域周长
		double perimeter_;
		//形心
		std::shared_ptr<MeshNode> centroid_sp_;
	};
}

