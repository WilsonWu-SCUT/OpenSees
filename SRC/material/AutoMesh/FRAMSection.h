#pragma once

namespace AutoMesh
{
	//�ʷֱ߽�
	class MeshBoundary;
	//�ʷ�
	class MeshCpp;
	//Ĭ���ݽ�λ��
	class MeshNode;

	class FRAMSection
	{
	public:
		FRAMSection();
		FRAMSection(const int& type, const SectionType& section_type);
		~FRAMSection();

	public:
		//�趨�ߴ���Ϣ
		bool set_dimension(const std::vector<double>& dimension_vec, const int& cover);
		//�趨�����Ϣ
		bool set_As(const std::vector<int>& As_vec, bool is_beam);

	public:
		//��ñ߽�ָ��
		std::shared_ptr<MeshCpp> get_mesh_sp(const int& mesh_size,
			const MeshMethod& method, bool auto_set_as);

	public:
		//�趨�ݽ�ֱ��������ȷ���ݽ
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
		//��ӻ������߽�
		bool append_concrete_bound() const;
		//����͸ֱ߽�
		bool append_steel_bound() const;

	private:
		//�趨�ݽ�λ��(Straight)
		void set_straight_bars(const double& yi, const double& yj, const double& zi, const double& zj, const int& Atol, bool includeCorner = false);
		//�趨�ݽ�λ��(Straight)
		void set_straight_bars(const double& yi, const double& yj, const double& zi, const double& zj, const int& Abasic, const int& n, bool includeCorner = false);
		//�趨�ݽ�λ��(ring)
		void set_ring_bars(const double& d, const int& Atol);
		//�趨���ݽ�λ��
		void set_beam_bars(const double& cover);
		//�趨���ݽ�λ��
		void set_column_bars(const double& cover);

	public:
		//��ý������߶�
		double get_H(void) const;
		//��ý��������
		double get_B(void) const;
		//�Ƿ�Ϊ��
		inline bool isBeam(void) const {
			return this->is_beam_;
		};

	private:
		//��ù��ֱָ�׼���ݸ�ʽ
		std::vector<double> get_I_dimetion_vec(void) const;
		//���ʮ�ֱָ�׼���ݸ�ʽ
		std::vector<double> get_X_dimetion_vec(void) const;
		//������͸ֱ�׼���ݸ�ʽ
		std::vector<double> get_Box_dimetion_vec(void) const;
		//��û��θֱ�׼���ݸ�ʽ
		std::vector<double> get_Ring_dimention_vec(void) const;
		//��þ��ν����׼���ݸ�ʽ
		std::vector<double> get_Rectangular_dimention_vec(void) const;
		//���Բ�ν����׼���ݸ�ʽ
		std::vector<double> get_Circle_dimention_vec(void) const;

	private:
		//�Ƿ��͸ֻ���������
		bool is_src_section(void) const;
		//�Ƿ�����λ������Ľ���
		bool is_rect_section_with_bar() const;
		//�Ƿ��Բ�λ������Ľ���
		bool is_circ_section_with_bar() const;
		//������Ϣ�Ƿ��걸
		bool is_dimension_correct(void) const;
		//�����Ϣ�Ƿ��걸
		bool is_As_correct(void) const;
		//�Ƿ����ĳ���͸�
		bool include_shape(const SectionBasicShapeType& steel_type) const;

	private:
		//��������
		int cover_;
		//��������
		int type_;
		//�����������
		SectionType section_type_;
		//����ߴ��б�
		std::vector<double> dimention_vec_;
		//�����Ϣ�б�
		std::vector<int> as_vec_;
		//�Ƿ�������
		bool is_beam_;

	private:
		//�ʷ������
		std::shared_ptr<MeshCpp> mesh_sp_;
		//�ڵ�λ��
		std::vector<std::shared_ptr<MeshNode>> as_pos_vec_;
		//�ڵ����
		std::vector<double> as_A_vec_;

	private:
		//�ݽ�ֱ��
		int bar_diameter = 20;
		//�ݽ�����Сֵ
		int bar_dist_min = 50;
	};
}

