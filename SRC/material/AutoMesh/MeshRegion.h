#pragma once

class CAutomesh2DUI;

namespace AutoMesh
{
	//�ʷֽڵ�
	class MeshNode;
	//�ʷֵ�Ԫ
	class MeshElement;

	//�ʷ�������
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
		//������
		int region_id_;
		//��������
		int mat_type_;
		//�ڵ��б�
		std::map<int, std::shared_ptr<MeshNode>> node_map;
		//��Ԫ�б�
		std::map<int, std::shared_ptr<MeshElement>> element_map;
		//�����ܳ�
		double perimeter_;
		//����
		std::shared_ptr<MeshNode> centroid_sp_;
	};
}

