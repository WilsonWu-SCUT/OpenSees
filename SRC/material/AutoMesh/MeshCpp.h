#pragma once

class CAutomesh2DUI;

namespace AutoMesh
{
	//线
	class MeshLine;
	//圆
	class MeshCircle;
	class MeshRectangular;
	//点
	class MeshNode;
	//剖分结果
	class MeshResult;
	//剖分边界
	class MeshBoundary;
	//剖分域
	class MeshRegion;

	class MeshCpp
	{
	public:
		MeshCpp();
		MeshCpp(const int& mesh_size,const MeshMethod& method, bool iclude_hard_point,
			bool remove_inner, bool out_mesh_file = false);
		~MeshCpp();

	public:
		//警告信息
		static std::vector<std::string> warning_vec;
		//错误信息
		static std::vector<std::string> error_vec;

	public:
		//添加圆
		void add_circle(const double& r, const double& x, const double& y,  const int& mat_tag = 1);
		//添加线
		void add_line(const double& xi, const double& xj, const double& yi, const double& yj);
		//添加点
		void add_point(const double& x, const double& y);
		//添加矩形
		void add_rect(const double& width, const double& height, 
			const double& offset_x, const double& offset_y, const int& mat_tag = 1);
		//添加环形
		void add_ring(const double& outer_radius, const double& inner_radius, 
			const double& offset_x, const double& offset_y, const int& mat_tag = 1);
		//添加不规则平行四边形
		void add_trapezoid(const double& width, const double& height, 
			const double& T, const double& F, const int& mat_tag = 1);
		//添加槽型
		void add_groove(const double& B, const double& H, const double& U, const double& T,
			const double& D, const double& F, const double& X, const double& Y, const int& mat_tag = 1);
		//添加多段线
		void add_polyLine(const std::vector<double>& x_vec, const std::vector<double>& y_vec, const int& mat_tag = 1);
		
	public:
		//是否清空内部区域
		inline void set_remove_inner(bool is_remove) { RemoveInner = is_remove ? 1 : 0; }
		//是否清空内部区域
		inline bool get_remove_inner(void) const { return RemoveInner; }
		//自动剖分
		bool mesh(void);
		std::shared_ptr<MeshResult> auto_mesh(void);

	public:
		std::shared_ptr<MeshRegion> get_region_sp(const int& id);
		std::vector<std::shared_ptr<MeshNode>> get_node_sp_vec();
		int get_region_num(void) const;

	private:
		bool is_mesh_valid(std::shared_ptr<CAutomesh2DUI> amg) const;

	private:
		//边界线
		std::vector<std::shared_ptr<MeshLine>> line_sp_vec;
		//边界圆
		std::vector<std::shared_ptr<MeshCircle>> circle_sp_vec;
		//边界矩形
		std::vector<std::shared_ptr<MeshRectangular>> rect_sp_vec;
		//边界点
		std::vector<std::shared_ptr<MeshNode>> point_sp_vec;
		//剖分方法指针
		std::shared_ptr<CAutomesh2DUI> amg_sp;
		//剖分结果
		std::shared_ptr<MeshResult> mesh_result_sp;
		//边界
		std::shared_ptr<MeshBoundary> bound_sp;

		//目标网格尺寸
		int TargetMeshSize;
		//最终网格尺寸
		int FinalMeshSize;
		//剖分方法
		MeshMethod method_;
		//是否考虑硬点
		bool IncludeHardPoint;
		//是否输出相关文件
		bool OutFile;
		//是否删除空洞
		bool RemoveInner;
	};

}

