#pragma once

namespace AutoMesh
{
	//型钢类型
	enum class SectionBasicShapeType
	{
		//工字钢
		IShaped = 0,

		//箱型
		BoxShaped = 1,

		//环型
		RingShaped = 2,

		//十字钢
		XShaped = 3,

		//矩形
		RectangularShaped = 4,

		//圆形
		CircleShaped = 5,
	};

	//材料类型
	enum class MatType
	{
		//混凝土
		Concrete = 0,

		//纵筋
		ReinforceBar = 1,

		//钢材
		Steel = 2,

		//边界混凝土
		CoverConcrete = 3,
	};

	//界面类型
	enum class SectionType
	{
		RC = 6,

		Steel = 5,
	};

	//剖分方法
	enum class MeshMethod
	{
		//四边形剖分
		LoopingQuad = 0,
		SubmappingQuad = 1,
		TriToQuad = 2,
		PavingQuad = 3,
		MappingQuad = 4,
		GridQuadTri = 5,
		DelaunayTri = 6,
		//三角形剖分
		LoopingTri = 7,
		AFMTri = 8,
		CADTri = 9,
	};
}
