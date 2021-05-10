#ifndef PMMSection_h
#define PMMSection_h

#include <memory>
#include <vector>
#include <algorithm>
#include <map>

#include <SectionForceDeformation.h>
#include "SectionAnalysisExport.h"

class AIDMMaterial;

class PMMSection: public SectionForceDeformation
{
public:
	PMMSection();
	PMMSection(const int& tag, const int& AIDMTag_3, const int& AIDMTag_2,
		const int& section_type, 
		const std::vector<double> dimension_vec,
		const std::vector<int> As_vec, bool is_beam,
		const int& fcu, const double& bar_fy, const double& steel_fy);
	PMMSection(const int& tag, std::shared_ptr<SectionAnalysis::Section> section_sp);
	~PMMSection();

public:
	int get_moment(const double& My, const double& Mz, const double& axial_load, bool isI);
	double A(void) const {
		return this->section_sp_->A();
	}
	double Iy(void) const {
		return this->section_sp_->Iy();
	}
	double Iz(void) const {
		return this->section_sp_->Iz();
	}
	double E(void) const {
		return this->section_sp_->E();
	}
	double Jx(void) const {
		return this->Iy() + this->Iz();
	}
	double G(void) const {
		return 0.4 * this->E();
	}

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
	//Section Analysis
	std::shared_ptr<SectionAnalysis::Section> section_sp_;
	//AIDM≤ƒ¡œ÷∏’Î »∆3÷·
	AIDMMaterial* AIDM_3_ptr;
	//AIDM≤ƒ¡œ÷∏’Î »∆2÷·
	AIDMMaterial* AIDM_2_ptr;
	//≤ƒ¡œ±‡∫≈
	int AIDM_3_Tag_;
	int AIDM_2_Tag_;

private:
	static Vector s;
	static Matrix ks;
};

#endif
