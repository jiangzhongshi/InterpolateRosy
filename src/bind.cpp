
#if PYBIND

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <boost/function_output_iterator.hpp>
#include <fstream>
#include <vector>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <igl/readOFF.h>

namespace MyCGAL{
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;
struct halfedge2edge
{
  halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges)
    : m_mesh(m), m_edges(edges)
  {}
  void operator()(const halfedge_descriptor& h) const
  {
    m_edges.push_back(edge(h, m_mesh));
  }
  const Mesh& m_mesh;
  std::vector<edge_descriptor>& m_edges;
};



template <typename P>
bool read_off(const Eigen::MatrixXd& V,
              const Eigen::MatrixXi& F,
              CGAL::Surface_mesh<P>& sm)
{
  using namespace CGAL;
  typedef Surface_mesh<P> Mesh;
  typedef typename Kernel_traits<P>::Kernel K;
  typedef typename K::Vector_3 Vector_3;
  typedef typename Mesh::Face_index Face_index;
  typedef typename Mesh::Vertex_index Vertex_index;
  typedef typename Mesh::size_type size_type;
  int n, f, e = 0;
  std::string off;

  n = V.rows();
  f = F.rows();

  sm.reserve(V.rows(), F.rows()*2, F.rows());
  std::vector<Vertex_index> vertexmap(n);
  P p;
  Vector_3 v;
  typename Mesh::template Property_map<Vertex_index,CGAL::Color> vcolor;
  typename Mesh::template Property_map<Vertex_index,Vector_3> vnormal;
  bool vcolored = false, v_has_normals = false;

  char ci;

  for(int i=0; i < n; i++){
    Vertex_index vi = sm.add_vertex(P(V(i,0), V(i,1), V(i,2)));
    vertexmap[i] = vi;
  }
  std::vector<Vertex_index> vr;
  size_type d, vi;
  bool fcolored = false;
  typename Mesh::template Property_map<Face_index,CGAL::Color> fcolor;

  for(int i=0; i < f; i++){
    d = 3;
    vr.resize(3);
    for(std::size_t j=0; j<d; j++){
      vi = F(i,j);
      vr[j] = vertexmap[vi];
    }
    Face_index fi = sm.add_face(vr);
    if(fi == sm.null_face())
    {
      sm.clear();
      return false;
    }
   
  }
  return true;
}

  template <typename P>
  bool write_off(const CGAL::Surface_mesh<P>& sm, Eigen::MatrixXd& V,
              Eigen::MatrixXi& F) {
    using namespace CGAL;
    typedef Surface_mesh<P> Mesh;
    typedef typename Mesh::Vertex_index Vertex_index;
    typedef typename Mesh::Face_index Face_index;

    V.resize(sm.number_of_vertices(), 3);
    F.resize(sm.number_of_faces(), 3);
    std::vector<int> reindex;
    reindex.resize(sm.num_vertices());
    int n = 0;
    BOOST_FOREACH(Vertex_index v, sm.vertices()){
      auto p = sm.point(v);
      V.row(n) << p.x(), p.y(), p.z();
      reindex[v]=n++;
    }

    int n_f = 0;
    BOOST_FOREACH(Face_index f, sm.faces()){
      // os << sm.degree(f);
      int fj = 0;
      BOOST_FOREACH(Vertex_index v, CGAL::vertices_around_face(sm.halfedge(f),sm)){
        F(n_f, fj) = reindex[v];
        fj++;
      }
      n_f ++;
    }
    return true;
  }
}

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

    m.def("isotropic_remeshing", [](const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double target_edge_length){
        using namespace MyCGAL;
        MyCGAL::Mesh mesh;
        read_off(V,F, mesh);
        // for(int i=0; i<V.rows(); i++)
        unsigned int nb_iter = 3;
          std::vector<edge_descriptor> border;
          PMP::border_halfedges(faces(mesh),
            mesh,
            boost::make_function_output_iterator(halfedge2edge(mesh, border)));
          PMP::split_long_edges(border, target_edge_length, mesh);
        PMP::isotropic_remeshing(
            faces(mesh),
            target_edge_length,
            mesh,
            PMP::parameters::number_of_iterations(nb_iter)
            .protect_constraints(true)//i.e. protect border, here
            );
        Eigen::MatrixXd V2;
        Eigen::MatrixXi F2;
        write_off(mesh, V2, F2);
        return std::make_tuple(V2,F2);
    });
}
#endif