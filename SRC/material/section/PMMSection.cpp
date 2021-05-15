#include "PMMSection.h"
#include "SectionAnalysisExport.h"
#include "AutoMesh.h"

using namespace AutoMesh;

#include <elementAPI.h>
#include "AIDMMaterial.h"

void* OPS_PMMSectionRectBeam()
{
	if (OPS_GetNumRemainingInputArgs() < 7) {
		opserr << "insufficient arguments for PMMSectionRectBeam\n";
		return 0;
	}

	// get tag
	int tag, AIDMtag_3, AIDMtag_2;
	int numData = 1;
	if (OPS_GetIntInput(&numData, &tag) < 0) return 0;
	if (OPS_GetIntInput(&numData, &AIDMtag_3) < 0) return 0;
	if (OPS_GetIntInput(&numData, &AIDMtag_2) < 0) return 0;

	//参数容器
	std::vector<double> dimension_vec;
	std::vector<int> As_vec = {0};
	int fcu;
	double fy;

	// get diemention
	numData = 2;
	double double_data[2];
	if (OPS_GetDoubleInput(&numData, &double_data[0]) < 0) return 0;
	for (auto info : double_data)  dimension_vec.push_back(info);

	// get As
	numData = 2;
	int int_data[2];
	if (OPS_GetIntInput(&numData, &int_data[0]) < 0) return 0;
	for (auto info : int_data)  As_vec.push_back(info);

	numData = 1;
	if (OPS_GetIntInput(&numData, &fcu) < 0) return 0;
	if (OPS_GetDoubleInput(&numData, &fy) < 0) return 0;


	return new PMMSection(tag, AIDMtag_3, AIDMtag_2, 1, 
		dimension_vec, As_vec, true, 
		fcu, fy, 345, AutoMesh::SectionType::RC);
}

PMMSection::PMMSection()
	:SectionForceDeformation(0, SEC_TAG_PMMSection), 
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	sectionHeight_3_(0), sectionHeight_2_(0),
	capacity_3pos_ini_(0), capacity_3neg_ini_(0),
	capacity_2pos_ini_(0), capacity_2neg_ini_(0),
	AIDMTag_2_(-1), AIDMTag_3_(-1)
{

}

PMMSection::PMMSection(const int& tag, const int& AIDMTag_3, const int& AIDMTag_2,
	const int& section_type,
	const std::vector<double> dimension_vec,
	const std::vector<int> As_vec, bool is_beam,
	const int& fcu, const double& bar_fy, const double& steel_fy, AutoMesh::SectionType matType)
	: SectionForceDeformation(tag, SEC_TAG_PMMSection),
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	sectionHeight_3_(0), sectionHeight_2_(0),
	capacity_3pos_ini_(0), capacity_3neg_ini_(0),
	capacity_2pos_ini_(0), capacity_2neg_ini_(0),
	AIDMTag_2_(-1), AIDMTag_3_(-1)
{
	//初始化截面
	if (!this->iniSection(section_type, dimension_vec, matType, fcu, steel_fy)) return;
	//设定配筋
	bool isSucess = this->FRAMSection_sp_->set_As(As_vec, is_beam);
	//失败直接返回
	if (!isSucess)
	{
		opserr << "PMMSection::PMMSection Error to set_As\n";
		return;
	}
	//计算受拉受压配筋比
	auto lammdaT_2_pos = 1;
	auto lammdaT_3_pos = is_beam ? As_vec[2] / As_vec[1] : 1;
	//获得配筋及配筋面积
	auto as_pos_vec = this->FRAMSection_sp_->get_as_pos_vec();
	auto as_A_vec = this->FRAMSection_sp_->get_as_A_vec();
	//遍历向量
	for (int i = 0; i < std::min(as_pos_vec.size(), as_A_vec.size()); i++)
		this->section_sp_->add_reinforced_bar(as_pos_vec[i]->get_x(),
			as_pos_vec[i]->get_y(), as_A_vec[i], bar_fy);
	//截面分析
	this->section_sp_->analysis(4);
	//计算承载力
	this->capacity_3pos_ini_ = this->section_sp_->get_moment(0, 180) * 1E6;
	this->capacity_3neg_ini_ = this->section_sp_->get_moment(0, 0) * 1E6;
	this->capacity_2pos_ini_ = this->section_sp_->get_moment(0, 270) * 1E6;
	this->capacity_2neg_ini_ = this->section_sp_->get_moment(0, 90) * 1E6;
	//获得AIDM指针
	auto material3_ptr = OPS_getUniaxialMaterial(AIDMTag_3);
	auto material2_ptr = OPS_getUniaxialMaterial(AIDMTag_2);
	//判断指针类型
	if (material3_ptr != 0)
	{
		if (material3_ptr->getClassTag() != MAT_TAG_AIDMMaterial)
		{
			opserr << "WARNING Material " << material3_ptr->getTag() <<
				"is not AIDMmaterial - section PMMSection\n";
			this->AIDM_3_ptr = new AIDMMaterial();
		}
		else
		{
			this->AIDM_3_ptr = dynamic_cast<AIDMMaterial*>(material3_ptr->getCopy());
			this->AIDMTag_3_ = AIDMTag_3;
		}
	}
	else this->AIDM_3_ptr = new AIDMMaterial();
	if (material2_ptr != 0)
	{
		if (material2_ptr->getClassTag() != MAT_TAG_AIDMMaterial)
		{
			opserr << "WARNING Material " << material2_ptr->getTag() <<
				"is not AIDMmaterial - section PMMSection\n";
			this->AIDM_2_ptr = new AIDMMaterial();
		}
		else
		{
			this->AIDM_2_ptr = dynamic_cast<AIDMMaterial*>(material2_ptr->getCopy());
			this->AIDMTag_2_ = AIDMTag_2;
		}
	}
	else this->AIDM_2_ptr = new AIDMMaterial();
	//设定参数
	this->AIDM_3_ptr->setlammdaT_pos(lammdaT_2_pos);
	this->AIDM_3_ptr->setCapcacity(this->capacity_3pos_ini_, this->capacity_3neg_ini_);
	this->AIDM_2_ptr->setlammdaT_pos(lammdaT_3_pos);
	this->AIDM_2_ptr->setCapcacity(this->capacity_2pos_ini_, this->capacity_2neg_ini_);
}

PMMSection::PMMSection(const int& tag, const int& section_type,
	const std::vector<double> dimension_vec, AutoMesh::SectionType matType)
	: SectionForceDeformation(tag, SEC_TAG_PMMSection),
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	sectionHeight_3_(0), sectionHeight_2_(0),
	capacity_3pos_ini_(0), capacity_3neg_ini_(0),
	capacity_2pos_ini_(0), capacity_2neg_ini_(0),
	AIDMTag_2_(-1), AIDMTag_3_(-1)
{
	//初始化截面
	if (!this->iniSection(section_type, dimension_vec, matType, 30, 345)) return;
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
	AIDMTag_2_(-1), AIDMTag_3_(-1)
{
}

PMMSection::~PMMSection()
{
	if (this->AIDM_3_ptr) delete this->AIDM_3_ptr;
	if (this->AIDM_2_ptr) delete this->AIDM_2_ptr;
}

bool PMMSection::iniSection(const int& section_type, const std::vector<double> dimension_vec,
	AutoMesh::SectionType& matType, const int& fcu, const double& steel_fy)
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
	//剖分截面
	auto mesh_ptr = this->FRAMSection_sp_->get_mesh_sp(50, AutoMesh::MeshMethod::LoopingTri, false);
	if (!mesh_ptr->mesh())
	{
		opserr << "PMMSection::PMMSection Error to mesh\n";
		return false;
	}
	//设定截面参数
	this->section_sp_.reset(new SectionAnalysis::Section(mesh_ptr, fcu, steel_fy));
	return true;
}

void PMMSection::setInitialK(const double& L, const int& factor)
{
	/*设定初始刚度*/
	this->initial_K2_ = this->E() * this->Iz() * factor / L;
	this->initial_K3_ = this->E() * this->Iy() * factor / L;
	/*AIDM刚度调整*/
	if (!this->ensureIniK) return;
	//调整初始刚度
	this->AIDM_2_ptr->setInitialK(this->initial_K2_);
	this->AIDM_3_ptr->setInitialK(this->initial_K3_);
}

void PMMSection::setLammda(const Vector& force_vec, bool is_I)
{
	//防止非法计算
	if (force_vec(2) == 0 || force_vec(8) == 0 || force_vec(1) == 0 || force_vec(7) == 0)
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
	if (this->AIDMTag_2_ == -1 && this->AIDMTag_3_ == -1) return;
	//计算内力
	auto my = (is_I ? force_vec(4) : force_vec(10) * -1) / 1E6;
	auto mz = (is_I ? force_vec(5) : force_vec(11) * -1) / 1E6;
	auto P = (is_I ? force_vec(0) : force_vec(6) * -1) / 1E3;
	//初始化承载力
	double myca_pos, myca_neg, mzca_pos, mzca_neg;
	//计算承载力 单偏压
	if (!this->consider_double_bending)
	{
		double nouse_m1, nouse_m1;
		this->section_sp_->get_moment(my, 0, P, myca_pos, myca_neg, nouse_m1, nouse_m1);
		this->section_sp_->get_moment(0, mz, P, nouse_m1, nouse_m1, mzca_pos, mzca_neg);
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
	//设定承载力
	this->AIDM_3_ptr->setCapcacity(myca_pos, myca_neg);
	this->AIDM_2_ptr->setCapcacity(mzca_pos, mzca_neg);
}

bool PMMSection::checkCapacity(const Vector& force_vec, const int& eleTag, bool is_I)
{
	//计算内力
	auto my = (is_I ? force_vec(4) : force_vec(10) * -1);
	auto mz = (is_I ? force_vec(5) : force_vec(11) * -1);
	//判断内力是否满足需求
	bool is2Bool = this->AIDM_2_ptr->checkCapacity(mz, eleTag);
	bool is3Bool = this->AIDM_3_ptr->checkCapacity(my, eleTag);
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

SectionForceDeformation* PMMSection::getCopy(void)
{
	auto section = new PMMSection(this->getTag());
	section->section_sp_ = this->section_sp_;
	section->FRAMSection_sp_ = this->FRAMSection_sp_;
	//重构AIDM指针
	section->AIDM_2_ptr = dynamic_cast<AIDMMaterial*>(this->AIDM_2_ptr->getCopy());
	section->AIDM_3_ptr = dynamic_cast<AIDMMaterial*>(this->AIDM_3_ptr->getCopy());
	section->AIDMTag_2_ = this->AIDMTag_2_;
	section->AIDMTag_3_ = this->AIDMTag_3_;
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
	return s;
}

const Vector& PMMSection::getStressResultant(void)
{
	return s;
}

const Matrix& PMMSection::getSectionTangent(void)
{
	//初始化
	ks.Zero();
	//获得初始刚度
	ks(0, 0) = this->initial_K2_;
	ks(1, 1) = this->initial_K3_;
	//判断是否为有效AIDM
	if (this->AIDM_2_ptr->isAvailabelAIDM())
		ks(0, 0) = this->AIDM_2_ptr->getTangent();
	if (this->AIDM_3_ptr->isAvailabelAIDM())
		ks(0, 0) = this->AIDM_3_ptr->getTangent();
	return ks;
}

const Matrix& PMMSection::getInitialTangent(void)
{
	//初始化
	ks.Zero();
	//获得初始刚度
	ks(0, 0) = this->initial_K2_;
	ks(1, 1) = this->initial_K3_;
	return ks;
}

int PMMSection::commitState(void)
{
	//AIDM CommitState
	int retVal = this->AIDM_2_ptr->commitState();
	retVal += this->AIDM_3_ptr->commitState();
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
	this->AIDM_3_ptr->setCapcacity(this->capacity_3pos_ini_, this->capacity_3neg_ini_);
	this->AIDM_3_ptr->revertToStart();
	this->AIDM_2_ptr->setCapcacity(this->capacity_2pos_ini_, this->capacity_2neg_ini_);
	return 0;
}

const ID& PMMSection::getType(void)
{
	return NULL;
}

int PMMSection::getOrder(void) const
{
	return 0;
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


Vector PMMSection::s(0);

Matrix PMMSection::ks(2, 2);

#pragma endregion


