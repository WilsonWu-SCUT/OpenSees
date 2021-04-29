#pragma once

class CAutomesh2DUI;

namespace AutoMesh
{
	//��
	class MeshLine;
	//Բ
	class MeshCircle;
	class MeshRectangular;
	//��
	class MeshNode;
	//�ʷֽ��
	class MeshResult;
	//�ʷֱ߽�
	class MeshBoundary;
	//�ʷ���
	class MeshRegion;

	class MeshCpp
	{
	public:
		MeshCpp();
		MeshCpp(const int& mesh_size,const MeshMethod& method, bool iclude_hard_point,
			bool remove_inner, bool out_mesh_file = false);
		~MeshCpp();

	public:
		//������Ϣ
		static std::vector<std::string> warning_vec;
		//������Ϣ
		static std::vector<std::string> error_vec;

	public:
		//���Բ
		void add_circle(const double& r, const double& x, const double& y,  const int& mat_tag = 1);
		//�����
		void add_line(const double& xi, const double& xj, const double& yi, const double& yj);
		//��ӵ�
		void add_point(const double& x, const double& y);
		//��Ӿ���
		void add_rect(const double& width, const double& height, 
			const double& offset_x, const double& offset_y, const int& mat_tag = 1);
		//��ӻ���
		void add_ring(const double& outer_radius, const double& inner_radius, 
			const double& offset_x, const double& offset_y, const int& mat_tag = 1);
		//��Ӳ�����ƽ���ı���
		void add_trapezoid(const double& width, const double& height, 
			const double& T, const double& F, const int& mat_tag = 1);
		//��Ӳ���
		void add_groove(const double& B, const double& H, const double& U, const double& T,
			const double& D, const double& F, const double& X, const double& Y, const int& mat_tag = 1);
		//��Ӷ����
		void add_polyLine(const std::vector<double>& x_vec, const std::vector<double>& y_vec, const int& mat_tag = 1);
		
	public:
		//�Ƿ�����ڲ�����
		inline void set_remove_inner(bool is_remove) { RemoveInner = is_remove ? 1 : 0; }
		//�Ƿ�����ڲ�����
		inline bool get_remove_inner(void) const { return RemoveInner; }
		//�Զ��ʷ�
		bool mesh(void);
		std::shared_ptr<MeshResult> auto_mesh(void);

	public:
		std::shared_ptr<MeshRegion> get_region_sp(const int& id);
		std::vector<std::shared_ptr<MeshNode>> get_node_sp_vec();
		int get_region_num(void) const;

	private:
		bool is_mesh_valid(std::shared_ptr<CAutomesh2DUI> amg) const;

	private:
		//�߽���
		std::vector<std::shared_ptr<MeshLine>> line_sp_vec;
		//�߽�Բ
		std::vector<std::shared_ptr<MeshCircle>> circle_sp_vec;
		//�߽����
		std::vector<std::shared_ptr<MeshRectangular>> rect_sp_vec;
		//�߽��
		std::vector<std::shared_ptr<MeshNode>> point_sp_vec;
		//�ʷַ���ָ��
		std::shared_ptr<CAutomesh2DUI> amg_sp;
		//�ʷֽ��
		std::shared_ptr<MeshResult> mesh_result_sp;
		//�߽�
		std::shared_ptr<MeshBoundary> bound_sp;

		//Ŀ������ߴ�
		int TargetMeshSize;
		//��������ߴ�
		int FinalMeshSize;
		//�ʷַ���
		MeshMethod method_;
		//�Ƿ���Ӳ��
		bool IncludeHardPoint;
		//�Ƿ��������ļ�
		bool OutFile;
		//�Ƿ�ɾ���ն�
		bool RemoveInner;
	};

}

