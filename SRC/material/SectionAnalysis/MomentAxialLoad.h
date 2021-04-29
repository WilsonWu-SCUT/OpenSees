#pragma once

namespace SectionAnalysis
{
	class MomentAxialLoad
	{
	public:
		MomentAxialLoad();
		MomentAxialLoad(const int& rotation, const double& moment, const double& axial_load, 
			const double& inital_strain, const double& min_strain, const double& max_strain);
		MomentAxialLoad(const std::shared_ptr<MomentAxialLoad>& load_sp,
			const int& moment_amplitude, const int& axial_load_amplitude);

		~MomentAxialLoad();

	public:
		inline int get_moment(void) const {
			return this->moment_;
		}
		inline int get_axial_load(void) const {
			return this->axial_load_;
		}
		inline double get_min_strain(void) const {
			return this->min_strain_;
		}
		inline double get_max_strain(void) const {
			return this->max_strain_;
		}
		inline double get_ini_strain(void) const {
			return this->initial_strain_;
		}

	public:
		std::shared_ptr<MomentAxialLoad> interpolation(std::shared_ptr<MomentAxialLoad>& j_load_sp,
			const int& axial_load);

	private:
		double min_strain_;
		double max_strain_;
		double initial_strain_;
		//KN M
		int moment_;
		//kN
		int axial_load_;
		int rotation_;
	};
}
