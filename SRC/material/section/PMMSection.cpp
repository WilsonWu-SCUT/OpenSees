#include "PMMSection.h"
#include "AutoMesh.h"
#include "SectionAnalysisExport.h"
#include <elementAPI.h>
#include "AIDMMaterial.h"

using namespace AutoMesh;

PMMSection* OPS_PMMSectionRC(const int& plasticSize, const int& elasticSize, const int& sec_type,
	const int& dimention_size, const int& as_size, bool isBeam)
{
	// get tag
	int tag;
	int numData = 1;
	if (OPS_GetIntInput(&numData, &tag) < 0) return 0;
	//参数是不是满足
	if (OPS_GetNumRemainingInputArgs() != (plasticSize - 1) && 
		OPS_GetNumRemainingInputArgs() != (plasticSize) &&
		OPS_GetNumRemainingInputArgs() != (elasticSize - 1)) 
	{
		opserr << "insufficient arguments for AIDMSection with Tag:" << tag << endln;
		return 0;
	}
	//判断是否未弹性
	bool isPlastic = OPS_GetNumRemainingInputArgs() > (elasticSize - 1);
	
	//参数容器
	std::vector<double> dimension_vec = {};
	std::vector<int> As_vec = {};
	//如果是梁
	if (isBeam) As_vec.push_back(0);

	//材料强度
	double fck;
	double fy;
	//面积配箍系数
	double lammdaSV_y = 0;
	double lammdaSV_z = 0;
	//轴向约束刚度
	int ark = 0;

	//尺寸
	numData = dimention_size;
	double* double_data = new double[dimention_size + 1 ];
	//读取尺寸信息
	if (OPS_GetDoubleInput(&numData, &double_data[0]) < 0) return 0;
	//读取尺寸信息
	for (int i = 0; i < dimention_size; i++)
	{
		if(isBeam && i == 1 && sec_type == 29)
		{
			dimension_vec.push_back(-double_data[i]);
		}
		else  dimension_vec.push_back(double_data[i]);
	}
	//是否塑性
	if (isPlastic)
	{
		//配筋
		numData = as_size;
		int* int_data = new int[as_size + 1];
		if (OPS_GetIntInput(&numData, &int_data[0]) < 0) return 0;
		for (int i = 0; i < as_size; i++)
				As_vec.push_back(int_data[i]);
		//材料强度
		numData = 1;
		if (OPS_GetDoubleInput(&numData, &fck) < 0) return 0;
		if (OPS_GetDoubleInput(&numData, &fy) < 0) return 0;
		//梁
		if (!isBeam)
		{
			if (OPS_GetDoubleInput(&numData, &lammdaSV_y) < 0) return 0;
		}
		if (OPS_GetDoubleInput(&numData, &lammdaSV_z) < 0) return 0;
		delete int_data;
		delete double_data;
		//创建截面
		auto section = new PMMSection(tag, sec_type, dimension_vec, As_vec, isBeam,
			fck, fy, 345, AutoMesh::SectionType::RC, lammdaSV_y, lammdaSV_z);
		//如果是梁
		if (isBeam && OPS_GetNumRemainingInputArgs() > 0)
		{
			if (OPS_GetIntInput(&numData, &ark) < 0) return 0;
			//如果为T形梁
			if (sec_type == 29)
			{
				auto Arect = std::abs(dimension_vec[0] * dimension_vec[1]);
				section->setARK(ark, Arect);
			}
			else section->setARK(ark, -1);
		}
		return section;
	}
	else
	{
		numData = 1;
		if (OPS_GetDoubleInput(&numData, &fck) < 0) return 0;
		delete double_data;
		return new PMMSection(tag, sec_type, dimension_vec, fck, AutoMesh::SectionType::RC);
	}
}

void* OPS_PMMSectionRectBeam()
{
	return OPS_PMMSectionRC(8, 4, 1, 2, 2, true);
}

void* OPS_PMMSectionRectColumn()
{
	return OPS_PMMSectionRC(10, 4, 1, 2, 3, false);
}

void* OPS_PMMSectionCirColumn()
{
	return OPS_PMMSectionRC(7, 3, 3, 1, 1, false);
}

void* OPS_PMMSectionTBeam()
{
	return OPS_PMMSectionRC(10, 6, 29, 4, 2, true);
}

PMMSection::PMMSection()
	:SectionForceDeformation(0, SEC_TAG_PMMSection), 
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	sectionHeight_3_(0), sectionHeight_2_(0),
	capacity_3pos_ini_(0), capacity_3neg_ini_(0),
	capacity_2pos_ini_(0), capacity_2neg_ini_(0),
	fck_(0), steel_fy_(0),
	A_(0), Iy_(0), Iz_(0), E_(0), G_(0), Jx_(0)

{
	//初始化指针
	this->AIDM_2_ptr = new AIDMMaterial();
	this->AIDM_3_ptr = new AIDMMaterial();
}

PMMSection::PMMSection(const int& tag, const int& section_type,
	const std::vector<double> dimension_vec, const std::vector<int> As_vec,
	bool is_beam, const double& fck, const double& bar_fy, const double& steel_fy,
	AutoMesh::SectionType matType, const double& lammdaSV_y, const double& lammdaSV_z)
	: SectionForceDeformation(tag, SEC_TAG_PMMSection),
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	sectionHeight_3_(0), sectionHeight_2_(0),
	capacity_3pos_ini_(0), capacity_3neg_ini_(0),
	capacity_2pos_ini_(0), capacity_2neg_ini_(0),
	fck_(fck), steel_fy_(steel_fy), as_vec_(0),
	A_(0), Iy_(0), Iz_(0), E_(0), G_(0), Jx_(0)
{
	//初始化截面
	this->as_vec_ = As_vec;
	if (!this->iniSection(section_type, dimension_vec, matType, bar_fy, is_beam)) return;
	//计算受拉受压配筋比
	auto lammdaT_2_pos = 1;
	auto lammdaT_3_pos = this->FRAMSection_sp_->isBeam() ? (double)this->as_vec_[2] / this->as_vec_[1] : 1;
	//截面分析
	this->section_sp_->analysis(4);
	//计算承载力
	this->capacity_3pos_ini_ = this->section_sp_->get_moment(0, 180) * 1E6;
	this->capacity_3neg_ini_ = this->section_sp_->get_moment(0, 0) * 1E6;
	this->capacity_2pos_ini_ = this->section_sp_->get_moment(0, 270) * 1E6;
	this->capacity_2neg_ini_ = this->section_sp_->get_moment(0, 90) * 1E6;
	//纵筋配筋特征值
	auto lammdaS = this->section_sp_->A(AutoMesh::MatType::ReinforceBar) / this->section_sp_->A()
		* bar_fy / this->fck_;
	//创建材料指针
	this->AIDM_2_ptr = lammdaSV_y <= 0 ? new AIDMMaterial() : new AIDMMaterial(lammdaSV_y, lammdaS, lammdaT_2_pos);
	this->AIDM_3_ptr = lammdaSV_z <= 0 ? new AIDMMaterial() : new AIDMMaterial(lammdaSV_z, lammdaS, lammdaT_3_pos);
	//一者不存在则不考虑双偏压
	if (lammdaSV_y <= 0 || lammdaSV_z <= 0 || this->FRAMSection_sp_->isBeam()) this->consider_double_bending = false;
	//设定参数
	this->AIDM_3_ptr->setCapcacity(this->capacity_3pos_ini_, this->capacity_3neg_ini_, 0);
	this->AIDM_3_ptr->updateSkeletonParams();
	this->AIDM_2_ptr->setCapcacity(this->capacity_2pos_ini_, this->capacity_2neg_ini_, 0);
	this->AIDM_2_ptr->updateSkeletonParams();
	//输过是柱 不杀死
	//if (!this->FRAMSection_sp_->isBeam())
	//{
	//	this->AIDM_2_ptr->NotKill();
	//	this->AIDM_3_ptr->NotKill();
	//}
}

PMMSection::PMMSection(const int& tag, const int& section_type,
	const std::vector<double> dimension_vec, const double& main_strength, AutoMesh::SectionType matType)
	: SectionForceDeformation(tag, SEC_TAG_PMMSection),
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	sectionHeight_3_(0), sectionHeight_2_(0),
	capacity_3pos_ini_(0), capacity_3neg_ini_(0),
	capacity_2pos_ini_(0), capacity_2neg_ini_(0),
	fck_(30), steel_fy_(345), as_vec_(0),
	A_(0), Iy_(0), Iz_(0), E_(0), G_(0), Jx_(0)
{
	//初始化
	this->as_vec_ = {};
	//不考虑双偏压
	this->consider_double_bending = false;
	//判断强度
	switch (matType)
	{
	case AutoMesh::SectionType::RC:
		this->fck_ = main_strength; break;
	case AutoMesh::SectionType::Steel:
		this->steel_fy_ = main_strength; break;
	default: break;
	}
	//初始化截面
	if (!this->iniSection(section_type, dimension_vec, matType, 400, true)) return;
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
	fck_(0), steel_fy_(0), as_vec_(0),
	A_(0), Iy_(0), Iz_(0), E_(0), G_(0), Jx_(0)
{
	//初始化
	this->as_vec_ = {};
	//初始化指针
	this->AIDM_2_ptr = new AIDMMaterial();
	this->AIDM_3_ptr = new AIDMMaterial();
}

PMMSection::~PMMSection()
{
	if (this->AIDM_3_ptr) delete this->AIDM_3_ptr;
	if (this->AIDM_2_ptr) delete this->AIDM_2_ptr;
}

bool PMMSection::iniSection(const int& section_type, const std::vector<double> dimension_vec,
	AutoMesh::SectionType& matType, const double& bar_fy,bool isbeam)
{
	//构造截面
	this->FRAMSection_sp_.reset(new AutoMesh::FRAMSection(section_type, matType));
	//设定截面基本参数（保护层默认20）
	if (!this->FRAMSection_sp_->set_dimension(dimension_vec, 0))
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
	this->section_sp_.reset(new SectionAnalysis::Section(mesh_ptr, this->fck_, 0.1 * this->fck_, this->steel_fy_));
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
	//设定截面参数
	this->SetSectionBaicsInfo();
	return true;
}

void PMMSection::setARK(const int& ark, const double& ARect)
{
	if (ARect > 0)
	{
		auto Ac = this->section_sp_->A(AutoMesh::MatType::Concrete) +
			this->section_sp_->A(AutoMesh::MatType::CoverConcrete);
		auto factor = Ac == 0 ? 1.0 : ARect / Ac;
		this->AIDM_3_ptr->setARK(ark, factor);
	}
	else this->AIDM_3_ptr->setARK(ark, 1.0);
}

void PMMSection::setInitialK(const double& L, const int& factor, bool ensureIniK)
{
	/*设定初始刚度*/
	auto initial_K2 = this->E_ * this->Iz_ * factor / L;
	auto initial_K3 = this->E_ * this->Iy_ * factor / L;
	//调整初始刚度
	this->AIDM_2_ptr->setInitialK(initial_K2, ensureIniK);
	this->AIDM_3_ptr->setInitialK(initial_K3, ensureIniK);
}

void PMMSection::SetSectionBaicsInfo()
{
	this->A_ = this->section_sp_->A();
	this->Iy_ = this->section_sp_->Iy();
	this->Iz_ = this->section_sp_->Iz();
	this->E_ = this->section_sp_->E();
	this->G_ = this->E_ * 0.4;
	this->Jx_ = this->Iy_ + this->Iz_;
}

void PMMSection::setLammda(const double& shearSpan3, const double& shearSpan2)
{
	//防止非法计算
	if (this->sectionHeight_3_ == 0 || this->sectionHeight_2_ == 0)
		return;
	//计算剪跨比
	auto shearspanRatio_3 = shearSpan3 / this->sectionHeight_3_;
	auto shearspanRatio_2 = shearSpan2 / this->sectionHeight_2_;
	//设定剪跨比
	if(shearspanRatio_2 > 0)
		this->AIDM_2_ptr->setLammda(shearspanRatio_2);
	if (shearspanRatio_3 > 0)
	this->AIDM_3_ptr->setLammda(shearspanRatio_3);
}

void PMMSection::setCapacity(const Vector& force_vec, const Vector& deformation_vec, bool is_I)
{
	//AIDM不存在
	if (this->as_vec_.size() == 0) return;
	//计算内力
	auto deform3 = is_I ? deformation_vec(3) : -deformation_vec(4);
	auto deform2 = is_I ? deformation_vec(1) : -deformation_vec(2);

	auto my = (is_I ? force_vec(4) : force_vec(10) * -1) / 1E6;
	auto mz = (is_I ? force_vec(5) : force_vec(11) * -1) / 1E6;
	auto P = (is_I ? -1 * force_vec(0) : force_vec(6)) / 1E3;

	//计算轴压系数
	auto Ac = this->section_sp_->A(AutoMesh::MatType::Concrete) +
		this->section_sp_->A(AutoMesh::MatType::CoverConcrete);
	auto fA = this->steel_fy_ * this->section_sp_->A(AutoMesh::MatType::Steel) + Ac *
		(P <= 0 ? this->fck_ : 0.1 * this->fck_);
	auto axialRatio = fA == 0 ? 0 : P / fA * 1E3;

	//计算约束轴压系数(仅RC梁考虑)
	if (this->FRAMSection_sp_->isBeam())
	{
		auto CALR = this->AIDM_3_ptr->GetALR();
		P += CALR * fA / 1E3;
		axialRatio += CALR;
	}

	//初始化承载力
	double myca_pos, myca_neg, mzca_pos, mzca_neg;
	//不考虑双偏压的强度
	double myca_pos_sc, myca_neg_sc, mzca_pos_sc, mzca_neg_sc;
	myca_pos_sc = this->section_sp_->get_moment(P, 180);
	myca_neg_sc = this->section_sp_->get_moment(P, 0);
	mzca_pos_sc = this->section_sp_->get_moment(P, 270);
	mzca_neg_sc = this->section_sp_->get_moment(P, 90);
	//计算承载力 单偏压
	if (!this->consider_double_bending)
	{
		myca_pos = myca_pos_sc;
		myca_neg = myca_neg_sc;
		mzca_pos = mzca_pos_sc;
		mzca_neg = mzca_neg_sc;
	}
	else
	{
		//卸载阶段不更新双偏压承载力
		if (this->AIDM_2_ptr->getLoadingType() == 3 || this->AIDM_3_ptr->getLoadingType() == 3)
			return;
		//二四象限不考虑
		if (deform3 * my <= 0 || deform2 * mz <= 0) return;
		//需更新双偏压承载力的内力边界
		auto my_limit = (my > 0 ? myca_pos_sc : myca_neg_sc) * this->double_bending_factor;
		auto mz_limit = (mz > 0 ? mzca_pos_sc : mzca_neg_sc) * this->double_bending_factor;
		//计算峰值承载力对应的变形
		auto deform3_c = this->AIDM_3_ptr->getCstrain(deform3 > 0);
		auto deform2_c = this->AIDM_2_ptr->getCstrain(deform2 > 0);
		auto deform3_ini = this->AIDM_3_ptr->getInitialStrain(deform3 > 0);
		auto deform2_ini = this->AIDM_2_ptr->getInitialStrain(deform2 > 0);
		//承载力系数：基于变形 基于力
		double myPosDeFactor = 0.0; double myNegDeFactor = 0.0;
		double mzPosDeFactor = 0.0; double mzNegDeFactor = 0.0;
		double myPosFoFactor = 0.0; double myNegFoFactor = 0.0;
		double mzPosFoFactor = 0.0; double mzNegFoFactor = 0.0;
		//不满足力的更新条件 也 不满足变形更新条件
		if (std::abs(my) <= my_limit || std::abs(mz) <= mz_limit)
		{
			//第二 第四 象限
			if (deform3 * my <= 0 || deform2 * mz <= 0)
				return;
			//变形太小
			if (std::abs(deform3) <= deform3_ini || std::abs(deform2) <= deform2_ini)
				return;
		}
		//根据内力判别双偏压
		this->section_sp_->get_moment(my, mz, P, myca_pos, myca_neg, mzca_pos, mzca_neg);
		//如果处于一 三象限
		if (deform3 * my > 0 && deform2 * mz > 0)
		{
			//变形大于峰值
			if (std::abs(deform3) >= deform3_c || std::abs(deform2) >= deform2_c)
			{
				//根据变形判别双偏压
				this->section_sp_->get_moment(deform3, deform2, P, myca_pos, myca_neg, mzca_pos, mzca_neg);
			}
			else
			{
				//获得基于力的承载力系数
				myPosFoFactor = myca_pos / myca_pos_sc;
				myNegFoFactor = myca_neg / myca_neg_sc;
				mzPosFoFactor = mzca_pos / mzca_pos_sc;
				mzNegFoFactor = mzca_neg / mzca_neg_sc;
				//根据变形判别双偏压
				this->section_sp_->get_moment(deform3, deform2, P, myca_pos, myca_neg, mzca_pos, mzca_neg);
				//获得基于变形的承载力系数
				myPosDeFactor = myca_pos / myca_pos_sc;
				myNegDeFactor = myca_neg / myca_neg_sc;
				mzPosDeFactor = mzca_pos / mzca_pos_sc;
				mzNegDeFactor = mzca_neg / mzca_neg_sc;
				//获得系数
				auto factor = std::max(std::abs(deform3) / deform3_c, std::abs(deform2) / deform2_c);
				//重构承载力
				myca_pos = myca_pos_sc * (myPosDeFactor * factor + myPosFoFactor * (1 - factor));
				myca_neg = myca_neg_sc * (myNegDeFactor * factor + myNegFoFactor * (1 - factor));
				mzca_pos = mzca_pos_sc * (mzPosDeFactor * factor + mzPosFoFactor * (1 - factor));
				mzca_neg = mzca_neg_sc * (mzNegDeFactor * factor + mzNegFoFactor * (1 - factor));
			}
		}
	}
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
	auto deform3 = is_I ? deformation_vec(3) : -deformation_vec(4);
	auto deform2 = is_I ? deformation_vec(1) : -deformation_vec(2);
	//更新刚度
	int retVal = this->AIDM_3_ptr->setTrialStrain(deform3);
	retVal += this->AIDM_2_ptr->setTrialStrain(deform2);
	//获得加载类型
	auto loadingType3 = this->AIDM_3_ptr->getLoadingType();
	auto loadingType2 = this->AIDM_2_ptr->getLoadingType();
	//如果加载类型不同
	if (loadingType3 != loadingType2)
	{
		//如果其中一者为重加载一者为骨架：走骨架
		if (loadingType3 == 1 && loadingType2 == 2)
			this->AIDM_2_ptr->ResetTrialStrain(deform2, loadingType3);
		else if (loadingType3 == 2 && loadingType2 == 1)
			this->AIDM_3_ptr->ResetTrialStrain(deform3, loadingType2);
	}
	//判断刚度是否相反
	else if (this->AIDM_2_ptr->getTangent() * this->AIDM_3_ptr->getTangent() < 0)
	{
		//如果均是重加载
		if (loadingType3 == 2 && loadingType2 == 2)
		{
			this->AIDM_3_ptr->ResetTrialStrain(deform3, 1);
			this->AIDM_2_ptr->ResetTrialStrain(deform2, 1);
		}
	}
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
		"AxialRatio" + descp,
		"ConstraintAR" + descp,
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
		this->AIDM_3_ptr->getAxialRatio(),
		this->AIDM_3_ptr->getCALR(),
	};
}

bool PMMSection::isKill(void) const
{
	return this->AIDM_2_ptr->isKill() || this->AIDM_3_ptr->isKill();
}

SectionForceDeformation* PMMSection::getCopy(void)
{
	//初始化
	auto section = new PMMSection(this->getTag());
	section->section_sp_ = this->section_sp_;
	section->FRAMSection_sp_ = this->FRAMSection_sp_;
	//基本变量
	section->sectionHeight_3_ = this->sectionHeight_3_;
	section->sectionHeight_2_ = this->sectionHeight_2_;
	//初始承载力
	section->capacity_3pos_ini_ = this->capacity_3pos_ini_;
	section->capacity_3neg_ini_ = this->capacity_3neg_ini_;
	section->capacity_2pos_ini_ = this->capacity_2pos_ini_;
	section->capacity_2neg_ini_ = this->capacity_2neg_ini_;
	//配筋
	section->as_vec_ = this->as_vec_;
	section->fck_ = this->fck_;
	section->steel_fy_ = this->steel_fy_;
	section->consider_double_bending = this->consider_double_bending;
	//截面参数
	section->SetSectionBaicsInfo();

	//重构AIDM指针
	section->AIDM_2_ptr = (AIDMMaterial*)this->AIDM_2_ptr->getCopy();
	section->AIDM_3_ptr = (AIDMMaterial*)this->AIDM_3_ptr->getCopy();
	section->AIDM_3_ptr->setCapcacity(this->capacity_3pos_ini_, this->capacity_3neg_ini_, 0);
	section->AIDM_3_ptr->updateSkeletonParams();
	section->AIDM_2_ptr->setCapcacity(this->capacity_2pos_ini_, this->capacity_2neg_ini_, 0);
	section->AIDM_2_ptr->updateSkeletonParams();

	//输过是柱 不杀死
	//if (!this->FRAMSection_sp_->isBeam())
	//{
	//	section->AIDM_2_ptr->NotKill();
	//	section->AIDM_3_ptr->NotKill();
	//}

	return section;
}

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