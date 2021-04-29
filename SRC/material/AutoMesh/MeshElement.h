#pragma once

namespace AutoMesh
{
	class MeshPointParam;
	class MeshNode;

	class MeshElement
	{
	public:
		MeshElement();
		MeshElement(const int& id, std::vector<std::shared_ptr<MeshNode>> node_sp_vec,
			const double& A, const double& x, const double& y);
		~MeshElement();

	public:
		inline int get_id(void) const 
		{ 
			return ID;  
		}
		inline double get_a(void) const {
			return this->a_;
		}
		double get_x(void) const;
		double get_y(void) const;
		inline int get_node_num(void) const
		{
			return this->node_sp_vec_.size();
		}
		std::shared_ptr<MeshNode> get_node_sp(const int& index);

	private:
		inline bool isTriangle(void) const {
			return this->node_sp_vec_.size() == 3;
		}

	private:
		//�ʷֵ�Ԫ���
		int ID;
		//�ڵ��б�
		std::vector<std::shared_ptr<MeshNode>> node_sp_vec_;
		//��Ԫ���
		double a_;
		//����
		std::shared_ptr<MeshNode> centriod_sp_;
	};
}

