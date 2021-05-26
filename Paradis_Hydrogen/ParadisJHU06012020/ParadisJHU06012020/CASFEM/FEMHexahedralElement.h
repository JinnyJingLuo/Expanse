// Ahmed M. Hussein

#ifndef FEMHEXAHEDRALELEMENT_H_
#define FEMHEXAHEDRALELEMENT_H_

#include "FEMElementGeometry.h"
#include "TriPatch.h"

namespace FEMSystem
{
	class FEMHexahedralElement : public FEMElementGeometry
	{
	public:
		FEMHexahedralElement();
		FEMHexahedralElement(const FEMHexahedralElement& oElement);
		~FEMHexahedralElement();
		FEMHexahedralElement& operator=(const FEMHexahedralElement& oElement);
		void Reset();
		FEMElementGeometryType GetType() const;
		unsigned int GetFacesCount() const;
		
		virtual double GetVolume() const;
		virtual Point GetPoint(const Matrix& oNaturalCoordinates) const;
		virtual vector<Point> GetPoints(const Matrix& oNaturalCoordinates) const;
		virtual void GenerateSurfacePatches(vector< vector<unsigned int> >& vviNodesIndices,vector<GenericNode*>& vpoMidPoints) const;
		virtual bool IsPointInside(Point* poPoint,vector<double>& vdNaturalCoordinates) const;
		virtual double GetDistanceToElementCenterPoint(Point* poPoint) const;
		virtual vector<double> GetNearestNodeNaturalCoordinates(Point* poPoint) const;
		virtual bool IsFaceOnSurface(const unsigned int& iFaceID) const;
		virtual Point GetFaceCenter(const unsigned int& iFaceID) const;
		virtual vector<double> GetNaturalCoordinates(Point* poPoint,const unsigned int& iMaxIterationsCount = 750) const;
		virtual bool IsSurfaceElement() const;
		virtual bool IsValid() const;
		virtual Point GetCenterPoint() const;
		virtual vector<unsigned int> GetSurfaceFacesIndices() const;
		virtual Vector GetFaceOuterNormal(const unsigned int& iFaceID,const Matrix& oNaturalCoordinates) const;
		virtual Point GetPointOnFace(const unsigned int& iFaceID,const Matrix& oNaturalCoordinates) const;
		virtual Matrix GetFaceNaturalCoordinates(const unsigned int& iFaceID,const double& dXi,const double& dEta) const;
		GenericNode* GetFaceQuadPatch(const unsigned int& iFaceIndex,vector<unsigned int>& viOriginalNodesIndices) const;
		unsigned int GetNodesCount() const;
		virtual void Read(FILE* fpFile,const vector<FEMNode*>* pvpoNodes);
		virtual void Write(FILE* fpFile) const;
		Matrix GetShapeFunctions(const double& dXi,const double& dEta,const double& dZeta) const;
		Matrix GetShapeFunctionsXYZDerivatives(const double& dXi,const double& dEta,const double& dZeta,double& dJacobian) const;
		double GetJacobianMatrixDeterminant(const double& dXi,const double& dEta,const double& dZeta) const;
		Matrix GetShapeFunctionsExtrapolationsMatrix() const;

		vector< vector<double> > GetBodyGaussPointsCoordinates() const;
		vector<double> GetBodyGaussPointsWeights() const;
		vector< vector<double> > GetFaceGaussPointsCoordinates(const unsigned int& iFaceIndex) const;
		vector<double> GetFaceGaussPointsWeights() const;
		bool IsMidSideNode(const unsigned int& iIndex) const;
		double GetMassLumpingFactor() const;
		virtual FEMElementGeometry* Clone() const;
				
	private:
	
	protected:
		void Initialize();
		Matrix GetShapeFunctionsDerivatives(const double& dXi,const double& dEta,const double& dZeta) const;
		Matrix GetJacobianMatrix(const double& dXi,const double& dEta,const double& dZeta) const;
		
		static Matrix EvaluateShapeFunctions(const double& dXi,const double& dEta,const double& dZeta);
		static Matrix GetShapeFunctionsExtrapolations();
		static Matrix GetShapeFunctionsExtrapolationsOld();
		static Matrix EvaluateExtrapolatedShapeFunctions(const double& dXi,const double& dEta,const double& dZeta);
		static vector< vector<double> > GenerateBodyGaussPointsCoordinates();
		static vector<double> GenerateBodyGaussPointsWeights();
		static vector< vector< vector<double> > > GenerateFaceGaussPointsCoordinates();
		static vector<double> GenerateFaceGaussPointsWeights();
		
		static Matrix m_oShapeFunctionsExtrapolations;
		static vector< vector<double> > m_vvdGaussPointsCoordinates;
		static vector<double> m_vdGaussPointsWeights;
		static vector< vector< vector<double> > > m_vvdFaceGaussPointsCoordinates;
		static vector<double> m_vdFaceGaussPointsWeights;
		static unsigned int NodesPerElement;
		static unsigned int FacesCount;
		static unsigned int GaussPointsCount;
		static unsigned int FaceGaussPointsCount;
	};
}

#endif



