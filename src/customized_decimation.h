#ifndef CUSTOMIZED_DECIMATION_H
#define CUSTOMIZED_DECIMATION_H

#include "timer.h"
#include <igl/readOBJ.h>
#include <igl/decimate.h>
#include <igl/max_faces_stopping_condition.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/collapse_edge.h>
#include <igl/connect_boundary_to_infinity.h>
#include <igl/decimate.h>
#include <igl/edge_flaps.h>
#include <igl/max_faces_stopping_condition.h>
#include <igl/per_vertex_point_to_plane_quadrics.h>
#include <igl/qslim_optimal_collapse_edge_callbacks.h>
#include <igl/quadric_binary_plus_operator.h>
#include <igl/remove_unreferenced.h>
#include <igl/slice.h>
#include <igl/slice_mask.h>

inline void composite_combine_mapping(const std::vector<std::pair<int,int>>& combined_pairs,
                                int nv,
                               Eigen::VectorXi& MG) 
{
    // verify with these two lines
    //std::vector<std::pair<int,int>> pairs {{3,5},{2,3},{1,2},{0,4}};
    //composite_combine_mapping(pairs, 6, I) ;
    MG = Eigen::VectorXi::LinSpaced(nv, 0,nv-1);
    auto valid =std::vector<bool>(nv,true);// Eigen::Matrix<bool, -1, 1>::Constant(nv,true);
    
    // get grouping and valid from collapsed_verts
    for (auto v0v1: combined_pairs) {
      auto v_min = std::min(v0v1.first,v0v1.second);
      auto v_max = std::max(v0v1.first,v0v1.second);

      valid[v_max] = false;
      assert(MG(v_min) == v_min && "v_min shouldn't be touched");
      MG(v_max) = v_min;
    }

    // recover from valid
    for(int i=0; i<nv; i++) {
        if (!valid[MG(i)]) {
            MG(i) = MG(MG(i));
        }
    }
}


inline bool qslim_to_fn(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const size_t max_m,
  Eigen::MatrixXd & U,
  Eigen::MatrixXi & G,
  Eigen::VectorXi & MG)
{
  using namespace igl;
   Eigen::VectorXi I,J;
  // Original number of faces
  const int orig_m = F.rows();
  // Tracking number of faces
  int m = F.rows();
  typedef Eigen::MatrixXd DerivedV;
  typedef Eigen::MatrixXi DerivedF;
  DerivedV VO;
  DerivedF FO;
  igl::connect_boundary_to_infinity(V,F,VO,FO);
  // decimate will not work correctly on non-edge-manifold meshes. By extension
  // this includes meshes with non-manifold vertices on the boundary since these
  // will create a non-manifold edge when connected to infinity.
  if(!is_edge_manifold(FO))
  {
    return false;
  }
  Eigen::VectorXi EMAP;
  Eigen::MatrixXi E,EF,EI;
  edge_flaps(FO,E,EMAP,EF,EI);
  // Quadrics per vertex
  typedef std::tuple<Eigen::MatrixXd,Eigen::RowVectorXd,double> Quadric;
  std::vector<Quadric> quadrics;
  per_vertex_point_to_plane_quadrics(VO,FO,EMAP,EF,EI,quadrics);
  // State variables keeping track of edge we just collapsed
  std::vector<std::pair<int,int>> collapsed_verts;
  int v1 = -1;
  int v2 = -1;
  // Callbacks for computing and updating metric
  std::function<void(
    const int e,
    const Eigen::MatrixXd &,
    const Eigen::MatrixXi &,
    const Eigen::MatrixXi &,
    const Eigen::VectorXi &,
    const Eigen::MatrixXi &,
    const Eigen::MatrixXi &,
    double &,
    Eigen::RowVectorXd &)> cost_and_placement;
  std::function<bool(
    const Eigen::MatrixXd &                                         ,/*V*/
    const Eigen::MatrixXi &                                         ,/*F*/
    const Eigen::MatrixXi &                                         ,/*E*/
    const Eigen::VectorXi &                                         ,/*EMAP*/
    const Eigen::MatrixXi &                                         ,/*EF*/
    const Eigen::MatrixXi &                                         ,/*EI*/
    const std::set<std::pair<double,int> > &                        ,/*Q*/
    const std::vector<std::set<std::pair<double,int> >::iterator > &,/*Qit*/
    const Eigen::MatrixXd &                                         ,/*C*/
    const int                                                        /*e*/
    )> pre_collapse;
  std::function<void(
    const Eigen::MatrixXd &                                         ,   /*V*/
    const Eigen::MatrixXi &                                         ,   /*F*/
    const Eigen::MatrixXi &                                         ,   /*E*/
    const Eigen::VectorXi &                                         ,/*EMAP*/
    const Eigen::MatrixXi &                                         ,  /*EF*/
    const Eigen::MatrixXi &                                         ,  /*EI*/
    const std::set<std::pair<double,int> > &                        ,   /*Q*/
    const std::vector<std::set<std::pair<double,int> >::iterator > &, /*Qit*/
    const Eigen::MatrixXd &                                         ,   /*C*/
    const int                                                       ,   /*e*/
    const int                                                       ,  /*e1*/
    const int                                                       ,  /*e2*/
    const int                                                       ,  /*f1*/
    const int                                                       ,  /*f2*/
    const bool                                                  /*collapsed*/
    )> post_collapse;
  qslim_optimal_collapse_edge_callbacks(
    E,quadrics,v1,v2, cost_and_placement, pre_collapse,post_collapse);


     post_collapse = [&v1,&v2,&quadrics,& collapsed_verts](
      const Eigen::MatrixXd &                                         ,   /*V*/
      const Eigen::MatrixXi &                                         ,   /*F*/
      const Eigen::MatrixXi &                                         ,   /*E*/
      const Eigen::VectorXi &                                         ,/*EMAP*/
      const Eigen::MatrixXi &                                         ,  /*EF*/
      const Eigen::MatrixXi &                                         ,  /*EI*/
      const std::set<std::pair<double,int> > &                        ,   /*Q*/
      const std::vector<std::set<std::pair<double,int> >::iterator > &, /*Qit*/
      const Eigen::MatrixXd &                                         ,   /*C*/
      const int                                                       ,   /*e*/
      const int                                                       ,  /*e1*/
      const int                                                       ,  /*e2*/
      const int                                                       ,  /*f1*/
      const int                                                       ,  /*f2*/
      const bool                                                  collapsed
      )->void
  {
    if(collapsed)
    {
      quadrics[v1<v2?v1:v2] = quadrics[v1] + quadrics[v2];
       collapsed_verts.push_back(std::make_pair(v1,v2));
    }
  };
  // Call to greedy decimator
  bool ret = decimate(
    VO, FO,
    cost_and_placement,
    max_faces_stopping_condition(m,orig_m,max_m),
    pre_collapse,
    post_collapse,
    E, EMAP, EF, EI,
    U, G, J, I);
  // Remove phony boundary faces and clean up
  const Eigen::Array<bool,Eigen::Dynamic,1> keep = (J.array()<orig_m);
  igl::slice_mask(Eigen::MatrixXi(G),keep,1,G);
  igl::slice_mask(Eigen::VectorXi(J),keep,1,J);
  Eigen::VectorXi _1,I2;
  igl::remove_unreferenced(Eigen::MatrixXd(U),Eigen::MatrixXi(G),U,G,_1,I2);
  igl::slice(Eigen::VectorXi(I),I2,1,I);
    composite_combine_mapping(collapsed_verts, V.rows(), MG);


  return ret;
}


inline bool decimate_to_fn(const Eigen::MatrixXd& V,
                    const Eigen::MatrixXi& F,
                    int max_m,
                    Eigen::MatrixXd& U,
                    Eigen::MatrixXi& G,
                    Eigen::VectorXi& MG)
{
    // modified pre/post collapse to add a MG based on igl::decimate
    using namespace Eigen;
    using namespace std;
    using namespace igl;

    Eigen::VectorXi I,J;

     // Original number of faces
  const int orig_m = F.rows();
  // Tracking number of faces
  int m = F.rows();
  typedef Eigen::MatrixXd DerivedV;
  typedef Eigen::MatrixXi DerivedF;
  DerivedV VO;
  DerivedF FO;
  igl::connect_boundary_to_infinity(V,F,VO,FO);
  // decimate will not work correctly on non-edge-manifold meshes. By extension
  // this includes meshes with non-manifold vertices on the boundary since these
  // will create a non-manifold edge when connected to infinity.
  if(!is_edge_manifold(FO))
  {
    return false;
  }

  std::vector<std::pair<int,int>> collapsed_verts;
  int v0, v1;
  const auto always_try = [&collapsed_verts, &v0,&v1](
    const Eigen::MatrixXd &                                         ,/*V*/
    const Eigen::MatrixXi &                                         ,/*F*/
    const Eigen::MatrixXi &E                                         ,/*E*/
    const Eigen::VectorXi &                                         ,/*EMAP*/
    const Eigen::MatrixXi &                                         ,/*EF*/
    const Eigen::MatrixXi &                                         ,/*EI*/
    const std::set<std::pair<double,int> > &                        ,/*Q*/
    const std::vector<std::set<std::pair<double,int> >::iterator > &,/*Qit*/
    const Eigen::MatrixXd &                                         ,/*C*/
    const int e                                                       /*e*/
    ) -> bool {
        // std::cout<<"pre"<<E(e,0)<<std::endl;
        v0 = E(e,0);
        v1 = E(e,1);
        return true;
     };
  const auto never_care = [&collapsed_verts,&v0,&v1](
    const Eigen::MatrixXd &                                         ,   /*V*/
    const Eigen::MatrixXi &                                         ,   /*F*/
    const Eigen::MatrixXi &                                         ,   /*E*/
    const Eigen::VectorXi &                                         ,/*EMAP*/
    const Eigen::MatrixXi &                                         ,  /*EF*/
    const Eigen::MatrixXi &                                         ,  /*EI*/
    const std::set<std::pair<double,int> > &                        ,   /*Q*/
    const std::vector<std::set<std::pair<double,int> >::iterator > &, /*Qit*/
    const Eigen::MatrixXd &                                         ,   /*C*/
    const int                                                       ,   /*e*/
    const int                                                       ,  /*e1*/
    const int                                                       ,  /*e2*/
    const int                                                       ,  /*f1*/
    const int                                                       ,  /*f2*/
    const bool collapsed                                                 /*collapsed*/
    )-> void { 
        if (collapsed) {
            collapsed_verts.push_back(std::make_pair(v0,v1));
        }
    };

  bool ret = decimate(
    VO,
    FO,
    shortest_edge_and_midpoint,
    max_faces_stopping_condition(m,orig_m,max_m),
    always_try,
    never_care,
    U,
    G,
    J,
    I);
  const Eigen::Array<bool,Eigen::Dynamic,1> keep = (J.array()<orig_m);
  igl::slice_mask(Eigen::MatrixXi(G),keep,1,G);
  igl::slice_mask(Eigen::VectorXi(J),keep,1,J);
  Eigen::VectorXi _1,I2;
  igl::remove_unreferenced(Eigen::MatrixXd(U),Eigen::MatrixXi(G),U,G,_1,I2);
  igl::slice(Eigen::VectorXi(I),I2,1,I);

  composite_combine_mapping(collapsed_verts, V.rows(), MG);
  return ret;
}

#endif