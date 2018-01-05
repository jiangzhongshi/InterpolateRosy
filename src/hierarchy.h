#pragma once

#include "common.h"
#include "aabb.h"
#include "adjacency.h"
#include <unordered_map>
#include <algorithm> 
#include <queue>
#include "global_types.h"
#include <set>

using namespace std;

struct MeshStats {
	AABB mAABB;
	Vector3f mWeightedCenter;
	double mAverageEdgeLength;
	double mMaximumEdgeLength;
	double mSurfaceArea;

	MeshStats() :
		mWeightedCenter(Vector3f::Zero()),
		mAverageEdgeLength(0.0f),
		mMaximumEdgeLength(0.0f),
		mSurfaceArea(0.0f) { }
};

class MultiResolutionHierarchy {
public:
    MultiResolutionHierarchy();

    bool load(const std::string &filename);
	bool load(const MatrixXf& V, const MatrixXu& F);

  //protected:
	void build(const Eigen::VectorXi &b,
                      const MatrixXf &bc,
                      bool use_boundary = true);
	void construct_tEs_tFEs(MatrixXu & F, std::vector<std::vector<uint32_t>> &mtFes, std::vector<tuple_E> &mtEs);

  MatrixXf &V(uint32_t i = 0) { return mV[i]; }
    const MatrixXf &V(uint32_t i = 0) const { return mV[i]; }

    MatrixXf &N(uint32_t i = 0) { return mN[i]; }
    const MatrixXf &N(uint32_t i = 0) const { return mN[i]; }

    MatrixXf &Q(uint32_t i = 0) { return mQ[i]; }
    const MatrixXf &Q(uint32_t i = 0) const { return mQ[i]; }

    MatrixXf &O(uint32_t i = 0) { return mO[i]; }
    const MatrixXf &O(uint32_t i = 0) const { return mO[i]; }

    MatrixXf &C(uint32_t i = 0) { return mC[i]; }
    const MatrixXf &C(uint32_t i = 0) const { return mC[i]; }

    MatrixXu &F() { return mF; }
    const MatrixXu &F() const { return mF; }

    MatrixXu &T() { return mT; }
    const MatrixXu &T() const { return mT; }

    SMatrix &L(uint32_t i = 0) { return mL[i]; }
    const SMatrix &L(uint32_t i) const { return mL[i]; }

    void smoothOrientationsTri(uint32_t l, bool alignment, bool randomization, bool extrinsic);

	void prolongOrientations(int level);

  AABB aabb() const { return mAABB; };;

    Float scale() const { return ratio_scale; }
	void setScale(Float scale) { 
		ratio_scale = scale; 
		mScale = diagonalLen * scale; 
    }


    int levels() const { return mL.size(); }
public:
	//for both 2D & 3D 
    std::vector<MatrixXf> mV;
    std::vector<MatrixXf> mN;
    std::vector<MatrixXf> mQ;
    std::vector<MatrixXf> mO;
    std::vector<MatrixXf> mC;
    std::vector<SMatrix> mL;
    std::vector<SMatrix> mP;
    MatrixXu mF;
    MatrixXu mT;

	std::vector<std::vector<uint32_t>> nFes;
	std::vector<tuple_E> nEs;
	vector<vector<bool>> nV_boundary_flag;
	std::vector<std::vector<uint32_t>> nV_nes;

	vector<vector<uint32_t>> vnfs;

    MatrixXf mNF, mCF;
    uint32_t mOrientationIterations;
    AABB mAABB;
    Float mScale, mInvScale;
	Float diagonalLen;
	Float ratio_scale;

public:
    statistics sta;
	bool triangles, splitting;

	//for t-mesh vertex tag
	vector<int> V_flag;
};
