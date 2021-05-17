#include "PMMSection.h"
#include "AutoMesh.h"
#include "SectionAnalysisExport.h"
#include <elementAPI.h>
#include "AIDMMaterial.h"

using namespace AutoMesh;



void* OPS_PMMSectionRectBeam()
{
	//参数是不是满足
	if (OPS_GetNumRemainingInputArgs() != 9 || 
		OPS_GetNumRemainingInputArgs() != 3) {
		opserr << "insufficient arguments for PMMSectionRectBeam\n";
		return 0;
	}
	//判断是否未弹性
	bool isPlastic = OPS_GetNumRemainingInputArgs() > 3;

	// get tag
	int tag;
	int numData = 1;
	if (OPS_GetIntInput(&numData, &tag) < 0) return 0;

	//参数容器
	std::vector<double> dimension_vec;
	std::vector<int> As_vec = {0};
	int fcu;
	double fy;
	double lammdaSV_y, lammdaSV_z;

	//尺寸
	numData = 2;
	double double_data[2];
	if (OPS_GetDoubleInput(&numData, &double_data[0]) < 0) return 0;
	for (auto info : double_data)  dimension_vec.push_back(info);

	if (isPlastic)
	{
		//配筋
		numData = 2;
		int int_data[2];
		if (OPS_GetIntInput(&numData, &int_data[0]) < 0) return 0;
		for (auto info : int_data)  As_vec.push_back(info);
		//材料强度
		numData = 1;
		if (OPS_GetIntInput(&numData, &fcu) < 0) return 0;
		if (OPS_GetDoubleInput(&numData, &fy) < 0) return 0;
		if (OPS_GetDoubleInput(&numData, &lammdaSV_y) < 0) return 0;
		if (OPS_GetDoubleInput(&numData, &lammdaSV_z) < 0) return 0;
		return new PMMSection(tag, 1, dimension_vec, As_vec, true,
			fcu, fy, 345, AutoMesh::SectionType::RC, lammdaSV_y, lammdaSV_z);
	}
	else return new PMMSection(tag, 1, dimension_vec, AutoMesh::SectionType::RC);	
}

PMMSection::PMMSection()
	:SectionForceDeformation(0, SEC_TAG_PMMSection), 
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	sectionHeight_3_(0), sectionHeight_2_(0),
	capacity_3pos_ini_(0), capacity_3neg_ini_(0),
	capacity_2pos_ini_(0), capacity_2neg_ini_(0),
	fcu_(0), steel_fy_(0)
{

}

PMMSection::PMMSection(const int& tag, const int& section_type,
	const std::vector<double> dimension_vec, const std::vector<int> As_vec,
	bool is_beam, const int& fcu, const double& bar_fy, const double& steel_fy,
	AutoMesh::SectionType matType, const double& lammdaSV_y, const double& lammdaSV_z)
	: SectionForceDeformation(tag, SEC_TAG_PMMSection),
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	sectionHeight_3_(0), sectionHeight_2_(0),
	capacity_3pos_ini_(0), capacity_3neg_ini_(0),
	capacity_2pos_ini_(0), capacity_2neg_ini_(0),
	fcu_(fcu), steel_fy_(steel_fy)
{
	//初始化截面
	this->as_vec_ = As_vec;
	if (!this->iniSection(section_type, dimension_vec, matType, fcu, bar_fy, steel_fy, is_beam)) return;
	//计算受拉受压配筋比
	auto lammdaT_2_pos = 1;
	auto lammdaT_3_pos = is_beam ? As_vec[2] / As_vec[1] : 1;
	//截面分析
	this->section_sp_->analysis(4);
	//计算承载力
	this->capacity_3pos_ini_ = this->section_sp_->get_moment(0, 180) * 1E6;
	this->capacity_3neg_ini_ = this->section_sp_->get_moment(0, 0) * 1E6;
	this->capacity_2pos_ini_ = this->section_sp_->get_moment(0, 270) * 1E6;
	this->capacity_2neg_ini_ = this->section_sp_->get_moment(0, 90) * 1E6;
	//纵筋配筋特征值
	auto lammdaS = this->section_sp_->A(AutoMesh::MatType::ReinforceBar) / this->section_sp_->A()
		* bar_fy / this->get_fck(fcu);
	//创建材料指针
	this->AIDM_2_ptr = lammdaSV_y <= 0 ? new AIDMMaterial() : new AIDMMaterial(lammdaSV_y);
	this->AIDM_3_ptr = lammdaSV_z <= 0 ? new AIDMMaterial() : new AIDMMaterial(lammdaSV_z);
	//设定参数
	this->AIDM_3_ptr->setlammdaT_pos(lammdaT_2_pos);
	this->AIDM_3_ptr->setlammdaS(lammdaS);
	this->AIDM_3_ptr->setCapcacity(this->capacity_3pos_ini_, this->capacity_3neg_ini_, 0);
	this->AIDM_3_ptr->updateSkeletonParams();
	this->AIDM_2_ptr->setlammdaT_pos(lammdaT_3_pos);
	this->AIDM_2_ptr->setlammdaS(lammdaS);
	this->AIDM_2_ptr->setCapcacity(this->capacity_2pos_ini_, this->capacity_2neg_ini_, 0);
	this->AIDM_2_ptr->updateSkeletonParams();
}

PMMSection::PMMSection(const int& tag, const int& section_type,
	const std::vector<double> dimension_vec, AutoMesh::SectionType matType)
	: SectionForceDeformation(tag, SEC_TAG_PMMSection),
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	sectionHeight_3_(0), sectionHeight_2_(0),
	capacity_3pos_ini_(0), capacity_3neg_ini_(0),
	capacity_2pos_ini_(0), capacity_2neg_ini_(0),
	fcu_(0), steel_fy_(0)
{
	//初始化截面
	if (!this->iniSection(section_type, dimension_vec, matType, 30, 400, 345, true)) return;
	//初始化指针
	this->AIDM_2_ptr = new AIDMMaterial();
	this->AIDM_3_ptr = new AIDMMaterial();
}

PMMSection::PMMSection(const int& tag)
	: SectionForceDeformation(tag, SEC_TAG_PMMSection),
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	sectionHeight_3_(0), sectionHeight_2_(0),
	capacity_3pos_ini_(0), capacity_3neg_ini_(0),
	capacity_2pos_ini_(0), capacity_2neg_ini_(0),
	fcu_(0), steel_fy_(0)
{
}

PMMSection::~PMMSection()
{
	if (this->AIDM_3_ptr) delete this->AIDM_3_ptr;
	if (this->AIDM_2_ptr) delete this->AIDM_2_ptr;
}

bool PMMSection::iniSection(const int& section_type, const std::vector<double> dimension_vec,
	AutoMesh::SectionType& matType, const int& fcu, const double& bar_fy, const double& steel_fy, bool isbeam)
{
	//构造截面
	this->FRAMSection_sp_.reset(new AutoMesh::FRAMSection(section_type, matType));
	//设定截面基本参数（保护层默认20）
	if (!this->FRAMSection_sp_->set_dimension(dimension_vec, 20))
	{
		opserr << "PMMSection::PMMSection Error to set_dimension\n";
		return false;
	}
	//计算截面高度
	this->sectionHeight_2_ = this->FRAMSection_sp_->get_B();
	this->sectionHeight_3_ = this->FRAMSection_sp_->get_H();
	//判断是否为有效配筋
	if (this->as_vec_.size() != 0)
	{
		if (!this->FRAMSection_sp_->set_As(this->as_vec_, isbeam))
		{
			opserr << "PMMSection::PMMSection Error to set_As\n";
			return false;
		}
		
	}
	//剖分截面
	auto mesh_ptr = this->FRAMSection_sp_->get_mesh_sp(50, AutoMesh::MeshMethod::LoopingTri, false);
	if (!mesh_ptr->mesh())
	{
		opserr << "PMMSection::PMMSection Error to mesh\n";
		return false;
	}
	//设定截面参数
	this->section_sp_.reset(new SectionAnalysis::Section(mesh_ptr, fcu, steel_fy));
	//获得配筋及配筋面积
	if (this->as_vec_.size() != 0)
	{
		auto as_pos_vec = this->FRAMSection_sp_->get_as_pos_vec();
		auto as_A_vec = this->FRAMSection_sp_->get_as_A_vec();
		//遍历向量
		for (int i = 0; i < std::min(as_pos_vec.size(), as_A_vec.size()); i++)
			this->section_sp_->add_reinforced_bar(as_pos_vec[i]->get_x(),
				as_pos_vec[i]->get_y(), as_A_vec[i], bar_fy);
	}
	return true;
}

void PMMSection::setInitialK(const double& L, const int& factor)
{
	/*设定初始刚度*/
	auto initial_K2 = this->E() * this->Iz() * factor / L;
	auto initial_K3 = this->E() * this->Iy() * factor / L;
	//调整初始刚度
	this->AIDM_2_ptr->setInitialK(initial_K2);
	this->AIDM_3_ptr->setInitialK(initial_K3);
}

void PMMSection::setLammda(const Vector& force_vec, bool is_I)
{
	//防止非法计算
	if (force_vec(2) == 0 || force_vec(8) == 0 || force_vec(1) == 0 || force_vec(7) == 0 || 
		this->sectionHeight_3_ == 0 || this->sectionHeight_2_ == 0)
		return;
	//计算剪跨
	auto shearSpan_3 = is_I ?
		force_vec(4) / force_vec(2) : force_vec(10) / force_vec(8);
	auto shearSpan_2 = is_I ?
		force_vec(5) / force_vec(1) : force_vec(11) / force_vec(7);
	//计算剪跨比
	auto shearspanRatio_3 = std::abs(shearSpan_3) / this->sectionHeight_3_;
	auto shearspanRatio_2 = std::abs(shearSpan_2) / this->sectionHeight_2_;
	//设定剪跨比
	this->AIDM_2_ptr->setLammda(shearspanRatio_2);
	this->AIDM_3_ptr->setLammda(shearspanRatio_3);
}

void PMMSection::setCapacity(const Vector& force_vec, bool is_I)
{
	//AIDM不存在
	if (this->as_vec_.size() == 0) return;
	//计算内力
	auto my = (is_I ? force_vec(4) : force_vec(10) * -1) / 1E6;
	auto mz = (is_I ? force_vec(5) : force_vec(11) * -1) / 1E6;
	auto P = (is_I ? force_vec(0) : force_vec(6) * -1) / 1E3;
	//初始化承载力
	double myca_pos, myca_neg, mzca_pos, mzca_neg;
	//计算承载力 单偏压
	if (!this->consider_double_bending)
	{
		myca_pos = this->section_sp_->get_moment(P, 180);
		myca_neg = this->section_sp_->get_moment(P, 0);
		mzca_pos = this->section_sp_->get_moment(P, 270);
		mzca_neg = this->section_sp_->get_moment(P, 90);
	}
	//考虑双偏压
	else this->section_sp_->get_moment(my, mz, P, myca_pos, myca_neg, mzca_pos, mzca_neg);
	//更新单位
	myca_pos *= 1E6; myca_neg *= 1E6; mzca_pos *= 1E6; mzca_neg *= 1E6;
	//防止承载力过小
	myca_pos = myca_pos < this->capacity_3pos_ini_* this->min_capacity_factor ? 
		this->capacity_3pos_ini_ * this->min_capacity_factor : myca_pos;

	myca_neg = myca_neg < this->capacity_3neg_ini_* this->min_capacity_factor ?
		this->capacity_3neg_ini_ * this->min_capacity_factor : myca_neg;

	mzca_pos = mzca_pos < this->capacity_2pos_ini_* this->min_capacity_factor ?
		this->capacity_2pos_ini_ * this->min_capacity_factor : mzca_pos;

	mzca_neg = mzca_neg < this->capacity_2neg_ini_* this->min_capacity_factor ?
		this->capacity_2neg_ini_ * this->min_capacity_factor : mzca_neg;
	//计算轴压系数
	auto Ac = this->section_sp_->A(AutoMesh::MatType::Concrete) +
		this->section_sp_->A(AutoMesh::MatType::CoverConcrete);
	auto fA = this->steel_fy_ * this->section_sp_->A(AutoMesh::MatType::Steel) + Ac * 
		(P <= 0 ? this->get_fck(this->fcu_): this->get_ftk(this->fcu_));
	auto axialRatio = fA == 0 ? 0 : P / fA;
	//设定承载力
	this->AIDM_3_ptr->setCapcacity(myca_pos, myca_neg, axialRatio);
	this->AIDM_2_ptr->setCapcacity(mzca_pos, mzca_neg, axialRatio);
}

bool PMMSection::checkCapacity(const Vector& force_vec, const int& eleTag, bool is_I)
{
	//计算内力
	auto my = (is_I ? force_vec(4) : force_vec(10) * -1);
	auto mz = (is_I ? force_vec(5) : force_vec(11) * -1);
	//判断内力是否满足需求
	bool is2Bool = this->AIDM_2_ptr->checkCapacity(mz);
	bool is3Bool = this->AIDM_3_ptr->checkCapacity(my);
	if (!is2Bool || !is3Bool)
		opserr << "Warinning: the capacity of AIDMBeamColumn " << eleTag << " with PMMSection " << this->getTag() <<
		" is too weak: convert to elastic automatically." << endln;
	return (is2Bool && is3Bool);
}

int PMMSection::setTrialDeformation(const Vector& deformation_vec, bool is_I)
{
	int retVal = this->AIDM_3_ptr->setTrialStrain(
		is_I ? deformation_vec(3): -deformation_vec(4));
	retVal += this->AIDM_2_ptr->setTrialStrain(
		is_I ? deformation_vec(1) : -deformation_vec(2));
	return retVal;
}

std::vector<std::string> PMMSection::getResponseStrVec(bool is_I)
{
	std::string descp = is_I ? "I" : "J";
	return
	{
		"strain2" + descp,
		"stress2" + descp,
		"lammda2" + descp,
		"strain3" + descp,
		"stress3" + descp,
		"lammda3" + descp,
	};
}

std::vector<double> PMMSection::getResponseVec()
{
	return
	{
		this->AIDM_2_ptr->getStrain(),
		this->AIDM_2_ptr->getStress(),
		this->AIDM_2_ptr->getLammda(),
		this->AIDM_3_ptr->getStrain(),
		this->AIDM_3_ptr->getStress(),
		this->AIDM_3_ptr->getLammda(),
	};
}

SectionForceDeformation* PMMSection::getCopy(void)
{
	auto section = new PMMSection(this->getTag());
	section->section_sp_ = this->section_sp_;
	section->FRAMSection_sp_ = this->FRAMSection_sp_;
	//重构AIDM指针
	section->AIDM_2_ptr = new AIDMMaterial(this->AIDM_2_ptr->getLammdaSV());
	section->AIDM_3_ptr = new AIDMMaterial(this->AIDM_3_ptr->getLammdaSV());
	//基本变量
	section->sectionHeight_3_ = this->sectionHeight_3_;
	section->sectionHeight_2_ = this->sectionHeight_2_;
	//初始承载力
	section->capacity_3pos_ini_ = this->capacity_3pos_ini_;
	section->capacity_3neg_ini_ = this->capacity_3neg_ini_;
	section->capacity_2pos_ini_ = this->capacity_2pos_ini_;
	section->capacity_2neg_ini_ = this->capacity_2neg_ini_;
	return section;
}

#pragma region Section basic information

double PMMSection::A(void) const
{
	return this->section_sp_->A();
}

double PMMSection::Iy(void) const
{
	return this->section_sp_->Iy();
}

double PMMSection::Iz(void) const
{
	return this->section_sp_->Iz();
}

double PMMSection::E(void) const
{
	return this->section_sp_->E();
}

double PMMSection::Jx(void) const
{
	return this->Iy() + this->Iz();
}

double PMMSection::G(void) const
{
	return  0.4 * this->E();
}

#pragma endregion

#pragma region unavailable

const Vector& PMMSection::getSectionDeformation(void)
{
	//初始化
	s.Zero();
	//获得初始刚度
	s(0) = this->AIDM_2_ptr->getStrain();
	s(1) = this->AIDM_3_ptr->getStrain();
	return s;
}

const Vector& PMMSection::getStressResultant(void)
{
	//初始化
	s.Zero();
	//获得初始刚度
	s(0) = this->AIDM_2_ptr->getStress();
	s(1) = this->AIDM_3_ptr->getStress();
	return s;
}

const Matrix& PMMSection::getSectionTangent(void)
{
	//初始化
	ks.Zero();
	ks(0, 0) = this->AIDM_2_ptr->getTangent();
	ks(1, 1) = this->AIDM_3_ptr->getTangent();	
	return ks;
}

const Matrix& PMMSection::getInitialTangent(void)
{
	//初始化
	ks.Zero();
	//获得初始刚度
	ks(0, 0) = this->AIDM_2_ptr->getInitialTangent();
	ks(1, 1) = this->AIDM_3_ptr->getInitialTangent();
	return ks;
}

int PMMSection::commitState(void)
{
	//AIDM CommitState
	int retVal = this->AIDM_2_ptr->commitState();
	retVal += this->AIDM_3_ptr->commitState();
	//是否杀死单元
	if (this->AIDM_2_ptr->isKill() || this->AIDM_3_ptr->isKill())
	{
		this->AIDM_2_ptr->Kill();
		this->AIDM_3_ptr->Kill();
	}
	//是否成功
	return retVal;
}

int PMMSection::revertToLastCommit(void)
{
	int retVal = this->AIDM_2_ptr->revertToLastCommit();
	retVal += this->AIDM_3_ptr->revertToLastCommit();
	return retVal;
}

int PMMSection::revertToStart(void)
{
	this->AIDM_2_ptr->revertToStart();
	this->AIDM_3_ptr->setCapcacity(this->capacity_3pos_ini_, this->capacity_3neg_ini_, 0);
	this->AIDM_3_ptr->revertToStart();
	this->AIDM_2_ptr->setCapcacity(this->capacity_2pos_ini_, this->capacity_2neg_ini_, 0);
	return 0;
}

const ID& PMMSection::getType(void)
{
	return NULL;
}

int PMMSection::getOrder(void) const
{
	return 2;
}

int PMMSection::sendSelf(int commitTag, Channel& theChannel)
{
	opserr << "PMMSection::sendSelf() - failed to send data\n";
	return -1;
}

int PMMSection::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	opserr << "PMMSection::recvSelf() - failed to receive data\n";
	return -1;
}

void PMMSection::Print(OPS_Stream& s, int flag)
{
	opserr << "PMMSection::Print() - failed to print\n";
}

Vector PMMSection::s(2);

Matrix PMMSection::ks(2, 2);

#pragma endregion

double PMMSection::get_fck(const double& fcu) const
{
	//(混规条文说明4.1.3)
		 //混凝土强度小于C40（标准强度）
	double afa_1 = 0.0;
	double afa_2 = 0.0;
	if (fcu <= 40)
	{
		afa_1 = 0.76;
		afa_2 = 1;
	}
	//混凝土强度小于C50（标准强度）
	else if (fcu <= 50)
	{
		afa_1 = 0.76;
		afa_2 = 1 - (0.13F) / 40 * (fcu - 40);
	}
	//混凝土强度小于C80（标准强度）
	else if (fcu <= 80)
	{
		afa_1 = 0.76 + (0.06) / 30 * (fcu - 50);
		afa_2 = 1 - (0.13) / 40 * (fcu - 40);
	}
	else
	{
		afa_1 = 0.82;
		afa_2 = 0.87;
	}
	//抗压强度标准值（混规条文说明4.1.3）
	auto fck = fcu * (0.88 * afa_1 * afa_2);
	fck = (int)(fck / 0.1) * 0.1;
	return fck;
}

double PMMSection::get_ftk(const double& fcu) const
{
	if (14 < fcu && fcu < 16)
		return 1.27;
	else if (19 < fcu && fcu < 21)
		return 1.54;
	else if (24 < fcu && fcu < 26)
		return 1.78;
	else if (29 < fcu && fcu < 31)
		return 2.01;
	else if (34 < fcu && fcu < 36)
		return 2.20;
	else if (39 < fcu && fcu < 41)
		return 2.39;
	else if (44 < fcu && fcu < 46)
		return 2.51;
	else if (49 < fcu && fcu < 51)
		return 2.64;
	else if (54 < fcu && fcu < 56)
		return 2.74;
	else if (59 < fcu && fcu < 61)
		return 2.85;
	else if (64 < fcu && fcu < 66)
		return 2.93;
	else if (69 < fcu && fcu < 71)
		return 2.99;
	else if (74 < fcu && fcu < 76)
		return 3.05;
	else if (79 < fcu && fcu < 81)
		return 3.11;
	else
	{
		auto ftk = 0.88 * pow(fcu, 0.55) * 0.395 * 0.924 * 0.87;
		ftk = (int)(ftk / 0.01) * 0.01;
		return ftk;
	}
}
