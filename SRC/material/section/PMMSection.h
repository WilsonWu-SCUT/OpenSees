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
	//通用弹塑性截面
	PMMSection(const int& tag, const int& section_type, 
		const std::vector<double> dimension_vec, const std::vector<int> As_vec, 
		bool is_beam, const double& fck, const double& bar_fy, const double& steel_fy,
		AutoMesh::SectionType matType, const double& lammdaSV_y, const double& lammdaSV_z);
	//弹性截面
	PMMSection(const int& tag, const int& section_type,
		const std::vector<double> dimension_vec, const double& main_strength, AutoMesh::SectionType matType);
	PMMSection(const int& tag);
	~PMMSection();

public:
	/*获得刚度*/
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
	//纤维截面初始化
	bool iniSection(const int& section_type, const std::vector<double> dimension_vec,
		AutoMesh::SectionType& matType, const double& bar_fy, bool isbeam);

public:
	//设定轴向约束刚度
	void setARK(const double& ark);
	//设定初始刚度
	void setInitialK(const double& L, const int& factor, bool ensureIniK);
	//设定剪跨比
	void setLammda(const double& shearSpan3, const double& shearSpan2);
	//设定承载力
	void setCapacity(const Vector& force_vec, const Vector& deformation_vec, bool is_I);
	//判断承载力是否越界
	bool checkCapacity(const Vector& force_vec, const int& eleTag, bool is_I);
	//设定变形
	int setTrialDeformation(const Vector& deformation_vec, bool is_I);
	//获得Response类型描述
	std::vector<std::string> getResponseStrVec(bool is_I);
	std::vector<double> getResponseVec();
	//截面是否被杀死
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
	//防止承载力发生显著退化（相对于无轴力下的弯矩）
	double min_capacity_factor = 0.2;
	//超过该强度方考虑双偏压影响
	double double_bending_factor = 0.1;
	

private:
	//Section Analysis
	std::shared_ptr<SectionAnalysis::Section> section_sp_;
	//截面信息
	std::shared_ptr<AutoMesh::FRAMSection> FRAMSection_sp_;

private:
	//AIDM材料指针 绕3轴
	AIDMMaterial* AIDM_3_ptr;
	//AIDM材料指针 绕2轴
	AIDMMaterial* AIDM_2_ptr;
	//截面高度
	double sectionHeight_3_;
	double sectionHeight_2_;
	//初始截面承载力（无轴力）
	double capacity_3pos_ini_;
	double capacity_3neg_ini_;
	double capacity_2pos_ini_;
	double capacity_2neg_ini_;
	
private:
	/*强度*/
	double fck_;
	double steel_fy_;
	/*配筋*/
	std::vector<int> as_vec_;
	//是否考虑双偏压
	bool consider_double_bending = true ;

private:
	static Vector s;
	static Matrix ks;
};

#endif
