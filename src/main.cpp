
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <boost/function_output_iterator.hpp>
#include <fstream>
#include <vector>
#include <igl/readOFF.h>

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

int main(int argc, char **argv) {
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
  const char* filename = (argc > 1) ? argv[1] : "data/pig.off";
  Mesh mesh;
  // if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
  //   std::cerr << "Not a valid input file." << std::endl;
  //   return 1;
  // }
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::readOBJ(filename, V,F);
  read_off(V,F, mesh);
  double target_edge_length = 0.008296;
  unsigned int nb_iter = 3;
  std::cout << "Split border...";
    std::vector<edge_descriptor> border;
    PMP::border_halfedges(faces(mesh),
      mesh,
      boost::make_function_output_iterator(halfedge2edge(mesh, border)));
    PMP::split_long_edges(border, target_edge_length, mesh);
  std::cout << "done." << std::endl;
  std::cout << "Start remeshing of " << filename
    << " (" << num_faces(mesh) << " faces)..." << std::endl;
  PMP::isotropic_remeshing(
      faces(mesh),
      target_edge_length,
      mesh,
      PMP::parameters::number_of_iterations(nb_iter)
      .protect_constraints(false)//i.e. protect border, here
      );
  std::cout << "Remeshing done." << std::endl;
  return 0;
}
