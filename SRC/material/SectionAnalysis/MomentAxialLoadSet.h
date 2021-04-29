#pragma once

namespace SectionAnalysis
{
	#define MOMENT_AXIAL_NORM_NUM 20

	class MomentAxialLoad;

	class MomentAxialLoadSet
	{
	public:
		MomentAxialLoadSet();
		MomentAxialLoadSet(const int& rotation);
		MomentAxialLoadSet(std::shared_ptr<MomentAxialLoadSet>& load_set_sp);
		~MomentAxialLoadSet();

	public:
		void add(std::shared_ptr<MomentAxialLoad> load_sp, bool is_pos);
		void add(const double& moment, const double& axial_load, const double& strain, bool is_pos);

	public:
		std::vector<std::shared_ptr<MomentAxialLoad>> get_force_sp_vec(bool is_pos);
		std::vector<std::shared_ptr<MomentAxialLoad>> get_force_sp_vec();

	public:
		int get_max_axial_load(bool is_pos);
		int get_min_axial_load(bool is_pos);
		int get_moment(const double& axial_load, bool is_pos);

	private:
		void initial_map(std::shared_ptr<MomentAxialLoadSet>& load_set_sp, bool is_pos);

	private:
		std::map<int, std::shared_ptr<MomentAxialLoad>> load_sp_map_pos_;
		std::map<int, std::shared_ptr<MomentAxialLoad>> load_sp_map_neg_;

	private:
		int rotation_;
		int basic_axial_load_;
	};
}
