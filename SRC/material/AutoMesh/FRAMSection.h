#pragma once

namespace AutoMesh
{
	//剖分边界
	class MeshBoundary;
	//剖分
	class MeshCpp;
	//默认纵筋位置
	class MeshNode;

	class FRAMSection
	{
	public:
		FRAMSection();
		FRAMSection(const int& type, const SectionType& section_type);
		~FRAMSection();

	public:
		//设定尺寸信息
		bool set_dimension(const std::vector<double>& dimension_vec, const int& cover);
		//设定配筋信息
		bool set_As(const std::vector<int>& As_vec, bool is_beam);

	public:
		//获得边界指针
		std::shared_ptr<MeshCpp> get_mesh_sp(const int& mesh_size,
			const MeshMethod& method, bool auto_set_as);

	public:
		//设定纵筋直径（用于确定纵筋）
		inline void set_bar_diameter(const int& dia) { bar_diameter = dia; }

	public:
		std::vector<std::shared_ptr<MeshNode>>& get_as_pos_vec(void)
		{
			return this->as_pos_vec_;
		}
		std::vector<double>& get_as_A_vec(void)
		{
			return this->as_A_vec_;
		}

	private:
		//添加混凝土边界
		bool append_concrete_bound() const;
		//添加型钢边界
		bool append_steel_bound() const;

	private:
		//设定纵筋位置(Straight)
		void set_straight_bars(const double& yi, const double& yj, const double& zi, const double& zj, const int& Atol, bool includeCorner = false);
		//设定纵筋位置(Straight)
		void set_straight_bars(const double& yi, const double& yj, const double& zi, const double& zj, const int& Abasic, const int& n, bool includeCorner = false);
		//设定纵筋位置(ring)
		void set_ring_bars(const double& d, const int& Atol);
		//设定梁纵筋位置
		void set_beam_bars(const double& cover);
		//设定柱纵筋位置
		void set_column_bars(const double& cover);

	public:
		//获得截面最大高度
		double get_H(void) const;
		//获得截面最大宽度
		double get_B(void) const;
		//是否为梁
		inline bool isBeam(void) const {
			return this->is_beam_;
		};

	private:
		//获得工字钢标准数据格式
		std::vector<double> get_I_dimetion_vec(void) const;
		//获得十字钢标准数据格式
		std::vector<double> get_X_dimetion_vec(void) const;
		//获得箱型钢标准数据格式
		std::vector<double> get_Box_dimetion_vec(void) const;
		//获得环形钢标准数据格式
		std::vector<double> get_Ring_dimention_vec(void) const;
		//获得矩形截面标准数据格式
		std::vector<double> get_Rectangular_dimention_vec(void) const;
		//获得圆形截面标准数据格式
		std::vector<double> get_Circle_dimention_vec(void) const;

	private:
		//是否型钢混凝土截面
		bool is_src_section(void) const;
		//是否带矩形混凝土的截面
		bool is_rect_section_with_bar() const;
		//是否带圆形混凝土的截面
		bool is_circ_section_with_bar() const;
		//截面信息是否完备
		bool is_dimension_correct(void) const;
		//配筋信息是否完备
		bool is_As_correct(void) const;
		//是否包含某类型钢
		bool include_shape(const SectionBasicShapeType& steel_type) const;

	private:
		//保护层厚度
		int cover_;
		//截面类型
		int type_;
		//截面材料类型
		SectionType section_type_;
		//截面尺寸列表
		std::vector<double> dimention_vec_;
		//配筋信息列表
		std::vector<int> as_vec_;
		//是否梁截面
		bool is_beam_;

	private:
		//剖分求解器
		std::shared_ptr<MeshCpp> mesh_sp_;
		//节点位置
		std::vector<std::shared_ptr<MeshNode>> as_pos_vec_;
		//节点面积
		std::vector<double> as_A_vec_;

	private:
		//纵筋直径
		int bar_diameter = 20;
		//纵筋间距最小值
		int bar_dist_min = 50;
	};
}

