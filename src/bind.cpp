
#if PYBIND

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "common.h"
#include "hierarchy.h"
#include "timer.h"
using MatrixXdR = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

bool decimate_to_fn(const Eigen::MatrixXd& V,
                    const Eigen::MatrixXi& F,
                    int max_m,
                    Eigen::MatrixXd& U,
                    Eigen::MatrixXi& G,
                    Eigen::VectorXi& MG);

MatrixXf rosy_process(const MatrixXf &V,
                      const MatrixXu &F,
                      const Eigen::VectorXi &b,
                      const MatrixXf &bc,
                      Float scale,
                      int smooth_iter) {

    MultiResolutionHierarchy mRes;
    Timer<> timer;
    timer.beginStage("data pre-processing");
    mRes.load(V,F);

    mRes.build(b,bc);
    assert(!mRes.mQ[0].hasNaN());

    mRes.setScale(scale);
    timer.endStage();

    timer.beginStage("rosy optimization");

    int mLevelIterations = 0;
    int mMaxIterations = smooth_iter;
    int mLevel = mRes.levels()-2;

    while (true) {
        mRes.smoothOrientationsTri(mLevel, true, true, true);
        mLevelIterations++;

        if (mLevelIterations >= mMaxIterations) {
            mLevelIterations = 0;
            if (mLevel == 0) {
                break;
            }

            mLevel--;
            mRes.prolongOrientations(mLevel);
        }
    }
    timer.endStage();
    return mRes.mQ[0];
}


bool qslim_to_fn(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const size_t max_m,
  Eigen::MatrixXd & U,
  Eigen::MatrixXi & G,
  Eigen::VectorXi & MG);

namespace py = pybind11;
PYBIND11_MODULE(rosy, m) {
    m.doc() = R"(as a test)";

    m.def("add_any", [](py::EigenDRef<Eigen::MatrixXd> x, int r, int c, double v) { x(r,c) += v; });
    m.def("smooth_field", [](const MatrixXf &V, const MatrixXu &F, const Eigen::VectorXi&b, 
                                const MatrixXf &bc, Float s, int iter){
        return rosy_process(V,F,b,bc, s,iter);
    });
    m.def("decimate_to_maxf", [](const Eigen::MatrixXd& V, const Eigen::MatrixXi&F, int max_m) {
       Eigen::MatrixXd U;
       Eigen::MatrixXi G;
       Eigen::VectorXi MG;
       decimate_to_fn(V,F,max_m, U,G,MG);
       return std::make_tuple(U,G,MG);
    });
        m.def("qslim_to_maxf", [](const Eigen::MatrixXd& V, const Eigen::MatrixXi&F, int max_m) {
       Eigen::MatrixXd U;
       Eigen::MatrixXi G;
       Eigen::VectorXi MG;
       qslim_to_fn(V,F,max_m, U,G,MG);
       return std::make_tuple(U,G,MG);
    });
}
#endif