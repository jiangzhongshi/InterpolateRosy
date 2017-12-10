#include <numeric>
#include "hierarchy.h"
#include "quat.h"
#include "timer.h"
#include "positions.h"

void MultiResolutionHierarchy::smoothOrientationsTri(uint32_t l, bool alignment, bool randomization, bool extrinsic) {
    const MatrixXf &N = mN[l];
    const SMatrix &L = mL[l];
    MatrixXf &Q = mQ[l];

    Timer<> timer;
    double error = 0;
    int nLinks = 0;
    MatrixXf Q_new(Q.rows(), Q.cols());
    tbb::spin_mutex mutex;

#if PARALLEL
    tbb::parallel_for(
        tbb::blocked_range<uint32_t>(0u, (uint32_t) L.outerSize(), GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
            std::vector<std::pair<uint32_t, Float>> neighbors;
            double errorLocal = 0;
            int nLinksLocal = 0;
            for (uint32_t k = range.begin(); k != range.end(); ++k) {
#else
            std::vector<std::pair<uint32_t, Float>> neighbors;
            double errorLocal = 0;
            int nLinksLocal = 0;
            for (uint32_t k = 0; k < (uint32_t) L.outerSize(); k++) {
#endif
                SMatrix::InnerIterator it(L, k);

                uint32_t i = it.row();
                Vector3f q_i = Vector3f::Zero();
                Vector3f n_i = N.col(i);

				if (nV_boundary_flag[l][i]) {
					Q_new.col(i) = Q.col(i);
					continue;
				}

                neighbors.clear();
                for (; it; ++it) {
                    uint32_t j = it.col();
                    if (i == j)
                        continue;
                    neighbors.push_back(std::make_pair(j, it.value()));
                }

                if (randomization && neighbors.size() > 0)
                    pcg32(mOrientationIterations, k)
                        .shuffle(neighbors.begin(), neighbors.end());

                for (auto n : neighbors) {
                    uint32_t j = n.first;
                    Float value = n.second;
                    Float dp;

                    Vector3f q_j = Q.col(j), n_j = N.col(j);
                    if (extrinsic) {
                        q_j = applyRotationExtrinsic((q_i == Vector3f::Zero()) ? Q.col(i) : q_i, n_i, q_j, n_j);
                        dp = Q.col(i).dot(applyRotation(Q.col(i), N.col(i), Q.col(j), N.col(j)));
                    } else {
                        q_j = applyRotation((q_i == Vector3f::Zero()) ? Q.col(i) : q_i, n_i, q_j, n_j);
                        dp = Q.col(i).dot(applyRotation(Q.col(i), N.col(i), Q.col(j), N.col(j)));
                    }

                    errorLocal += std::abs(std::acos(std::min(dp, (Float) 1)));
                    ++nLinksLocal;

                    q_i += q_j * value;
                }

                if (q_i != Vector3f::Zero())
                    Q_new.col(i) = (q_i - n_i.dot(q_i) * n_i).normalized();
            }
            tbb::spin_mutex::scoped_lock guard(mutex);
            error += errorLocal;
            nLinks += nLinksLocal;
#if PARALLEL
        }
    );
#else
#endif
    mOrientationIterations++;
    Q = std::move(Q_new);
}

void MultiResolutionHierarchy::prolongOrientations(int level) {
	const SMatrix &P = mP[level];
	const MatrixXf &N = mN[level];
	for (int k = 0; k < P.outerSize(); ++k) {
		SMatrix::InnerIterator it(P, k);
		for (; it; ++it) {
			if (!tetMesh()) {
				if (nV_boundary_flag[level][it.row()]) continue;
				Vector3f q_j = mQ[level + 1].col(it.col());
				Vector3f n_i = N.col(it.row());
				mQ[level].col(it.row()) = q_j - n_i * n_i.dot(q_j);
			}
			else {
				Quaternion q_j = mQ[level + 1].col(it.col());
				Vector3f n_i = N.col(it.row());
				if (n_i != Vector3f::Zero()) {
					Float magnitude = q_j.norm();
					q_j = Quaternion(q_j / magnitude).align(n_i) * magnitude;
				}
				mQ[level].col(it.row()) = q_j;
			}
		}
	}
}

