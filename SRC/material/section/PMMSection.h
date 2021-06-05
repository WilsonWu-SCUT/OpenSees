#ifndef PMMSection_h
#define PMMSection_h

#include <memory>
#include <vector>
#include <algorithm>
#include <map>

#include <SectionForceDeformation.h>


class AIDMMaterial;
class Response;

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
	//ͨ�õ����Խ���
	PMMSection(const int& tag, const int& section_type, 
		const std::vector<double> dimension_vec, const std::vector<int> As_vec, 
		bool is_beam, const double& fck, const double& bar_fy, const double& steel_fy,
		AutoMesh::SectionType matType, const double& lammdaSV_y, const double& lammdaSV_z);
	//���Խ���
	PMMSection(const int& tag, const int& section_type,
		const std::vector<double> dimension_vec, const double& main_strength, AutoMesh::SectionType matType);
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
	inline double EIzoverL(const double& L, const int& factor) const
	{
		return factor * this->Iz() * this->E() / L;
	}
	inline double EIyoverL(const double& L, const int& factor) const
	{
		return factor * this->Iy() * this->E() / L;
	}

private:
	double A(void) const;
	double Iy(void) const;
	double Iz(void) const;
	double E(void) const;
	double Jx(void) const;
	double G(void) const;

private:
	//��ά�����ʼ��
	bool iniSection(const int& section_type, const std::vector<double> dimension_vec,
		AutoMesh::SectionType& matType, const double& bar_fy, bool isbeam);

public:
	//�趨����Լ���ն�
	void setARK(const double& ark);
	//�趨��ʼ�ն�
	void setInitialK(const double& L, const int& factor, bool ensureIniK);
	//�趨�����
	void setLammda(const double& shearSpan3, const double& shearSpan2);
	//�趨������
	void setCapacity(const Vector& force_vec, const Vector& deformation_vec, bool is_I);
	//�жϳ������Ƿ�Խ��
	bool checkCapacity(const Vector& force_vec, const int& eleTag, bool is_I);
	//�趨����
	int setTrialDeformation(const Vector& deformation_vec, bool is_I);
	//���Response��������
	std::vector<std::string> getResponseStrVec(bool is_I);
	std::vector<double> getResponseVec();
	//�����Ƿ�ɱ��
	bool isKill(void) const;

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
	//������ǿ�ȷ�����˫ƫѹӰ��
	double double_bending_factor = 0.1;
	

private:
	//Section Analysis
	std::shared_ptr<SectionAnalysis::Section> section_sp_;
	//������Ϣ
	std::shared_ptr<AutoMesh::FRAMSection> FRAMSection_sp_;

private:
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
	
private:
	/*ǿ��*/
	double fck_;
	double steel_fy_;
	/*���*/
	std::vector<int> as_vec_;
	//�Ƿ���˫ƫѹ
	bool consider_double_bending = true ;

private:
	static Vector s;
	static Matrix ks;
};

#endif
