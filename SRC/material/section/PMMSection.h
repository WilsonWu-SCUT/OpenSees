#ifndef PMMSection_h
#define PMMSection_h

#include <memory>
#include <vector>
#include <algorithm>
#include <map>

#include <SectionForceDeformation.h>


class AIDMMaterial;

namespace AutoMesh
{
	class FRAMSection;
	enum class SectionType;
}

namespace SectionAnalysis
{
	class Section;
}


class PMMSection: public SectionForceDeformation
{
public:
	PMMSection();
	PMMSection(const int& tag, const int& AIDMTag_3, const int& AIDMTag_2,
		const int& section_type, 
		const std::vector<double> dimension_vec,
		const std::vector<int> As_vec, bool is_beam,
		const int& fcu, const double& bar_fy, const double& steel_fy,
		AutoMesh::SectionType matType);
	PMMSection(const int& tag, const int& section_type,
		const std::vector<double> dimension_vec, AutoMesh::SectionType matType);
	PMMSection(const int& tag);
	~PMMSection();

public:
	/*��øն�*/
	inline double EAoverL(const double& L) const
	{
		return this->A() * this->E() / L;
	}
	inline double GJoverL(const double& L) const
	{
		return this->G() * this->Jx() / L;
	}

private:
	double A(void) const;
	double Iy(void) const;
	double Iz(void) const;
	double E(void) const;
	double Jx(void) const;
	double G(void) const;

private:
	//��ά�����ʼ����
	bool iniSection(const int& section_type, const std::vector<double> dimension_vec,
		AutoMesh::SectionType& matType, const int& fcu, const double& steel_fy);

public:
	//�趨��ʼ�ն�
	void setInitialK(const double& L, const int& factor);
	//�趨�����
	void setLammda(const Vector& force_vec, bool is_I);
	//�趨������
	void setCapacity(const Vector& force_vec, bool is_I);
	//�жϳ������Ƿ�Խ��
	bool checkCapacity(const Vector& force_vec, const int& eleTag, bool is_I);
	//�趨����
	int setTrialDeformation(const Vector& deformation_vec, bool is_I);

public:
	const Vector& getSectionDeformation(void);
	const Vector& getStressResultant(void);
	const Matrix& getSectionTangent(void);
	const Matrix& getInitialTangent(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	SectionForceDeformation* getCopy(void);
	const ID& getType(void);
	int getOrder(void) const;

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
	void Print(OPS_Stream& s, int flag);


private:
	//��ֹ���������������˻���������������µ���أ�
	double min_capacity_factor = 0.2;
	//�Ƿ���˫ƫѹ
	bool consider_double_bending = true;
	//�Ƿ�Ӳ����ʼ�ն�
	bool ensureIniK = true;

private:
	//Section Analysis
	std::shared_ptr<SectionAnalysis::Section> section_sp_;
	//������Ϣ
	std::shared_ptr<AutoMesh::FRAMSection> FRAMSection_sp_;

	//AIDM����ָ����
	int AIDMTag_2_;
	int AIDMTag_3_;
	//AIDM����ָ�� ��3��
	AIDMMaterial* AIDM_3_ptr;
	//AIDM����ָ�� ��2��
	AIDMMaterial* AIDM_2_ptr;
	//����߶�
	double sectionHeight_3_;
	double sectionHeight_2_;
	//��ʼ�������������������
	double capacity_3pos_ini_;
	double capacity_3neg_ini_;
	double capacity_2pos_ini_;
	double capacity_2neg_ini_;
	//��ʼ�ն�
	double initial_K3_;
	double initial_K2_;

private:
	static Vector s;
	static Matrix ks;
};

#endif
