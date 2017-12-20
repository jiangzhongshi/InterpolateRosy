
#include <igl/decimate.h>
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/is_edge_manifold.h>
#include <igl/remove_unreferenced.h>
#include <igl/slice_mask.h>
#include <igl/slice.h>
#include <igl/connect_boundary_to_infinity.h>
#include <igl/max_faces_stopping_condition.h>
#include <igl/shortest_edge_and_midpoint.h>


bool decimate_by_face(
  const Eigen::MatrixXd & OV,
  const Eigen::MatrixXi & OF,
  int max_face,
  Eigen::MatrixXd & U,
  Eigen::MatrixXi & G,
  Eigen::VectorXi & J,
  Eigen::VectorXi & I
  )
{
  using namespace igl;

  // Decimate 1
  using namespace Eigen;
  using namespace std;
  // Working copies
  Eigen::MatrixXd V = OV;
  Eigen::MatrixXi F = OF;
  Eigen::MatrixXi E = OE;
  Eigen::VectorXi EMAP = OEMAP;
  Eigen::MatrixXi EF = OEF;
  Eigen::MatrixXi EI = OEI;
  typedef std::set<std::pair<double,int> > PriorityQueue;
  PriorityQueue Q;
  std::vector<PriorityQueue::iterator > Qit;
  Qit.resize(E.rows());
  // If an edge were collapsed, we'd collapse it to these points:
  MatrixXd C(E.rows(),V.cols());
  for(int e = 0;e<E.rows();e++)
  {
    double cost = e;
    RowVectorXd p(1,3);
    cost_and_placement(e,V,F,E,EMAP,EF,EI,cost,p);
    C.row(e) = p;
    Qit[e] = Q.insert(std::pair<double,int>(cost,e)).first;
  }
  int prev_e = -1;
  bool clean_finish = false;

  while(true)
  {
    if(Q.empty())
    {
      break;
    }
    if(Q.begin()->first == std::numeric_limits<double>::infinity())
    {
      // min cost edge is infinite cost
      break;
    }
    int e,e1,e2,f1,f2;
    if(collapse_edge(
       cost_and_placement, pre_collapse, post_collapse,
       V,F,E,EMAP,EF,EI,Q,Qit,C,e,e1,e2,f1,f2))
    {
      if(stopping_condition(V,F,E,EMAP,EF,EI,Q,Qit,C,e,e1,e2,f1,f2))
      {
        clean_finish = true;
        break;
      }
    }else
    {
      if(prev_e == e)
      {
        assert(false && "Edge collapse no progress... bad stopping condition?");
        break;
      }
      // Edge was not collapsed... must have been invalid. collapse_edge should
      // have updated its cost to inf... continue
    }
    prev_e = e;
  }
  // remove all IGL_COLLAPSE_EDGE_NULL faces
  MatrixXi F2(F.rows(),3);
  J.resize(F.rows());
  int m = 0;
  for(int f = 0;f<F.rows();f++)
  {
    if(
      F(f,0) != IGL_COLLAPSE_EDGE_NULL || 
      F(f,1) != IGL_COLLAPSE_EDGE_NULL || 
      F(f,2) != IGL_COLLAPSE_EDGE_NULL)
    {
      F2.row(m) = F.row(f);
      J(m) = f;
      m++;
    }
  }
  F2.conservativeResize(m,F2.cols());
  J.conservativeResize(m);
  VectorXi _1;
  remove_unreferenced(V,F2,U,G,_1,I);
  return clean_finish;
}
