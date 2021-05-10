#pragma once

namespace AutoMesh
{
	class MeshCpp;
}

namespace SectionAnalysis
{
	class Fiber;
	class FiberSet;
	class MomentAxialLoad;
	class MomentAxialLoadSet;

	#define PM_DENSE_NUM 5
	#define SUB_MESH_AIXIAL_LOAD_TOL 0.3

	class Section
	{
	public:
		Section();
		Section(std::shared_ptr<AutoMesh::MeshCpp> mesh_sp,
			const double& fck, const double& ftk, const double& steel_fy);
		Section(std::shared_ptr<AutoMesh::MeshCpp> mesh_sp,
			const int& fcu,const double& steel_fy);
		Section(AutoMesh::MeshCpp* mesh_sp,
			const int& fcu, const double& steel_fy);
		~Section();

	public:
		void analysis(const int& circle_num);
		void add_reinforced_bar(const double& y, const double& z,
			const double& A, const double& fy);

	public:
		std::vector<std::shared_ptr<MomentAxialLoad>> get_force_sp_vec(const int& theta, bool is_pos);
		int get_moment(const double& axial_load, int theta);
		int get_moment(const double& My, const double& Mz, const double& axial_load, bool isI);
		int get_axial_load(int theta, const double& min_strain, const double& max_strain, 
			const double& min_z, const double& max_z);

	private:
		void pm_analysis(int theta);
		void semi_pm_analysis(int theta, std::shared_ptr<MomentAxialLoadSet>& force_set_sp, bool is_pos);
		void rapid_semi_pm_analysis(int theta, std::shared_ptr<MomentAxialLoadSet>& force_set_sp, bool is_pos);
		void CSI_semi_pm_analysis(int theta, std::shared_ptr<MomentAxialLoadSet>& force_set_sp, bool is_pos);

	private:
		void rapid_mesh_prt(std::shared_ptr<MomentAxialLoadSet>& force_set_sp, bool is_pos,
			std::shared_ptr<MomentAxialLoad>& force_i_sp, std::shared_ptr<MomentAxialLoad>& force_j_sp,
			const double& min_z, const double& max_z);
		void rapid_add_prt(std::shared_ptr<MomentAxialLoadSet>& force_set_sp, bool is_pos,
			std::shared_ptr<MomentAxialLoad>& force_i_sp, std::shared_ptr<MomentAxialLoad>& force_j_sp,
			const double& min_z, const double& max_z);

	public:
		double A(void) const;
		double Iy(void) const;
		double Iz(void) const;
		double E(void)const;

	private:
		//初始纤维截面
		std::shared_ptr<FiberSet> ini_fiber_set_sp_;
		//用于计算的纤维截面
		std::shared_ptr<FiberSet> target_fiber_set_sp_;
		//分析结果指针
		std::map<int, std::shared_ptr<MomentAxialLoadSet>> load_set_sp_map_;
		//快速截面分析时需要细分轴力点的轴力大小
		int sub_mesh_axial_load_;
		//获得角度值
		int delta_theta_;
	};
}
