// Distorted-grid adapter for the verified official LibShell MidedgeTan code.
//
// This file changes no shell formula.  It only replaces the procedural
// uniform reference vertices in the historical benchmark wrappers by an
// externally supplied vertex table and, for distributed loads, supplied
// continuum nodal integration weights.

#include "ElasticShell.h"
#include "MeshConnectivity.h"
#include "MidedgeAngleTanFormulation.h"
#include "RestState.h"
#include "StVKMaterial.h"
#include "types.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using SFF = LibShell::MidedgeAngleTanFormulation;
struct Grid { Eigen::MatrixXd vertices; Eigen::MatrixXi faces; };

Grid read_grid(int n, const std::string& path) {
    Grid grid;
    const int nv = (n + 1) * (n + 1);
    grid.vertices.resize(nv, 3);
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open vertex table");
    for (int v = 0; v < nv; ++v) {
        if (!(input >> grid.vertices(v,0) >> grid.vertices(v,1)
                    >> grid.vertices(v,2))) {
            throw std::runtime_error("incomplete vertex table");
        }
    }
    auto vertex = [n](int i, int j) { return j * (n + 1) + i; };
    grid.faces.resize(2*n*n, 3);
    int face = 0;
    for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) {
        const int q00=vertex(i,j), q10=vertex(i+1,j);
        const int q11=vertex(i+1,j+1), q01=vertex(i,j+1);
        if ((i+j)%2 == 0) {
            grid.faces.row(face++) << q00,q10,q11;
            grid.faces.row(face++) << q00,q11,q01;
        } else {
            grid.faces.row(face++) << q00,q10,q01;
            grid.faces.row(face++) << q10,q11,q01;
        }
    }
    return grid;
}

Eigen::VectorXd read_weights(int n, const std::string& path) {
    const int nv = (n + 1) * (n + 1);
    Eigen::VectorXd weights(nv);
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open weight table");
    for (int v = 0; v < nv; ++v) {
        if (!(input >> weights[v])) throw std::runtime_error("incomplete weight table");
    }
    return weights;
}

Eigen::SparseMatrix<double> from_triplets(
        int dimension, const std::vector<Eigen::Triplet<double>>& trips) {
    Eigen::SparseMatrix<double> result(dimension, dimension);
    result.setFromTriplets(trips.begin(), trips.end());
    result.makeCompressed();
    return result;
}

LibShell::MonolayerRestState rest_state(
        const LibShell::MeshConnectivity& mesh, const Eigen::MatrixXd& vertices,
        const Eigen::VectorXd& edge, double young, double poisson,
        double thickness) {
    LibShell::MonolayerRestState rest;
    rest.thicknesses.resize(mesh.nFaces(), thickness);
    rest.lameAlpha.resize(mesh.nFaces(), young*poisson/(1.0-poisson*poisson));
    rest.lameBeta.resize(mesh.nFaces(), young/(2.0*(1.0+poisson)));
    LibShell::ElasticShell<SFF>::firstFundamentalForms(mesh,vertices,rest.abars);
    LibShell::ElasticShell<SFF>::secondFundamentalForms(mesh,vertices,edge,rest.bbars);
    return rest;
}

Eigen::SparseMatrix<double> stiffness(
        const LibShell::MeshConnectivity& mesh, const Eigen::MatrixXd& vertices,
        const Eigen::VectorXd& edge, const LibShell::MonolayerRestState& rest,
        bool bending_only = false) {
    LibShell::StVKMaterial<SFF> material;
    std::vector<Eigen::Triplet<double>> trips;
    if (!bending_only) {
        LibShell::ElasticShell<SFF>::elasticEnergy(
            mesh,vertices,edge,material,rest,nullptr,&trips,
            LibShell::HessianProjectType::kNone);
    } else {
        LibShell::ElasticShell<SFF>::elasticEnergy(
            mesh,vertices,edge,material,rest,
            LibShell::ElasticShell<SFF>::ET_BENDING,nullptr,&trips,
            LibShell::HessianProjectType::kNone);
    }
    return from_triplets(3*vertices.rows()+SFF::numExtraDOFs*mesh.nEdges(), trips);
}

Eigen::SparseMatrix<double> principal_submatrix(
        const Eigen::SparseMatrix<double>& full, const std::vector<int>& active,
        std::vector<int>& inverse) {
    inverse.assign(full.rows(), -1);
    for (int k=0;k<int(active.size());++k) inverse[active[k]]=k;
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(full.nonZeros());
    for(int k=0;k<full.outerSize();++k)
        for(Eigen::SparseMatrix<double>::InnerIterator it(full,k);it;++it) {
            const int row=inverse[it.row()], col=inverse[it.col()];
            if(row>=0 && col>=0) trips.emplace_back(row,col,it.value());
        }
    return from_triplets(active.size(), trips);
}

Eigen::VectorXd solve_spd(const Eigen::SparseMatrix<double>& matrix,
                          const Eigen::VectorXd& force) {
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(matrix);
    if(solver.info()!=Eigen::Success) throw std::runtime_error("factorization failed");
    Eigen::VectorXd result=solver.solve(force);
    if(solver.info()!=Eigen::Success || !result.allFinite())
        throw std::runtime_error("linear solve failed");
    return result;
}

int plate(int n, const std::string& vertices_path, const std::string& weights_path, double thickness=0.01) {
    Grid grid=read_grid(n,vertices_path);
    const Eigen::VectorXd weights=read_weights(n,weights_path);
    LibShell::MeshConnectivity mesh(grid.faces);
    Eigen::VectorXd edge; SFF::initializeExtraDOFs(edge,mesh,grid.vertices);
    auto rest=rest_state(mesh,grid.vertices,edge,1.0,0.3,thickness);
    const auto full=stiffness(mesh,grid.vertices,edge,rest,true);
    const int nv=grid.vertices.rows();
    auto vertex=[n](int i,int j){return j*(n+1)+i;};
    std::vector<int> active;
    for(int j=1;j<n;++j)for(int i=1;i<n;++i)active.push_back(3*vertex(i,j)+2);
    for(int e=0;e<SFF::numExtraDOFs*mesh.nEdges();++e)active.push_back(3*nv+e);
    std::vector<int> inverse;
    const auto matrix=principal_submatrix(full,active,inverse);
    Eigen::VectorXd force=Eigen::VectorXd::Zero(active.size());
    for(int j=1;j<n;++j)for(int i=1;i<n;++i)
        force[inverse[3*vertex(i,j)+2]]=1.0e-8*weights[vertex(i,j)];
    const Eigen::VectorXd u=solve_spd(matrix,force);
    const double observable=u[inverse[3*vertex(n/2,n/2)+2]];
    std::cout<<std::setprecision(17)
      <<"{\"benchmark\":\"plate\",\"method\":\"MidedgeTan\",\"n\":"<<n
      <<",\"raw_scalar_states\":"<<full.rows()<<",\"active_states\":"<<active.size()
      <<",\"tangent_nnz\":"<<matrix.nonZeros()<<",\"observable\":"<<observable<<"}\n";
    return 0;
}

int roof(int n, const std::string& vertices_path, const std::string& weights_path, double thickness=0.25) {
    Grid grid=read_grid(n,vertices_path);
    const Eigen::VectorXd weights=read_weights(n,weights_path);
    LibShell::MeshConnectivity mesh(grid.faces);
    Eigen::VectorXd edge; SFF::initializeExtraDOFs(edge,mesh,grid.vertices);
    auto rest=rest_state(mesh,grid.vertices,edge,4.32e8,0.0,thickness);
    const auto full=stiffness(mesh,grid.vertices,edge,rest);
    const int nv=grid.vertices.rows(), ndof=full.rows();
    auto vertex=[n](int i,int j){return j*(n+1)+i;};
    std::vector<char> fixed(ndof,0);
    for(int j=0;j<=n;++j){
        fixed[3*vertex(0,j)]=1;
        fixed[3*vertex(n,j)+1]=1; fixed[3*vertex(n,j)+2]=1;
    }
    for(int e=0;e<SFF::numExtraDOFs*mesh.nEdges();++e){
        const int me=e/SFF::numExtraDOFs;
        const int v0=mesh.edgeVertex(me,0),v1=mesh.edgeVertex(me,1);
        if(v0%(n+1)==0&&v1%(n+1)==0)fixed[3*nv+e]=1;
    }
    std::vector<int> active;
    for(int d=0;d<ndof;++d)if(!fixed[d])active.push_back(d);
    std::vector<int> inverse;
    const auto matrix=principal_submatrix(full,active,inverse);
    Eigen::VectorXd force=Eigen::VectorXd::Zero(active.size());
    for(int j=0;j<=n;++j)for(int i=0;i<=n;++i){
        const int a=inverse[3*vertex(i,j)+2];
        if(a>=0)force[a]=-90.0*weights[vertex(i,j)];
    }
    const Eigen::VectorXd u=solve_spd(matrix,force);
    const double observable=-u[inverse[3*vertex(0,n)+2]];
    std::cout<<std::setprecision(17)
      <<"{\"benchmark\":\"roof\",\"method\":\"MidedgeTan\",\"n\":"<<n
      <<",\"raw_scalar_states\":"<<full.rows()<<",\"active_states\":"<<active.size()
      <<",\"tangent_nnz\":"<<matrix.nonZeros()<<",\"observable\":"<<observable<<"}\n";
    return 0;
}

int m3(int n, const std::string& vertices_path, double thickness=0.04) {
    Grid grid=read_grid(n,vertices_path);
    LibShell::MeshConnectivity mesh(grid.faces);
    Eigen::VectorXd edge; SFF::initializeExtraDOFs(edge,mesh,grid.vertices);
    auto rest=rest_state(mesh,grid.vertices,edge,6.825e7,0.3,thickness);
    const auto full=stiffness(mesh,grid.vertices,edge,rest);
    const int nv=grid.vertices.rows(), ndof=full.rows();
    auto vertex=[n](int i,int j){return j*(n+1)+i;};
    auto pdof=[&](int i,int j,int c){return 3*vertex(i,j)+c;};
    std::vector<char> fixed(ndof,0),dependent(ndof,0);
    for(int j=0;j<=n;++j){
        fixed[pdof(0,j,1)]=1; fixed[pdof(n,j,0)]=1;
        for(int c:{0,2})dependent[pdof(0,j,c)]=1;
        for(int c:{1,2})dependent[pdof(n,j,c)]=1;
    }
    fixed[pdof(n/2,0,2)]=1;
    for(int e=0;e<SFF::numExtraDOFs*mesh.nEdges();++e){
        const int me=e/SFF::numExtraDOFs;
        const int v0=mesh.edgeVertex(me,0),v1=mesh.edgeVertex(me,1);
        const int i0=v0%(n+1),i1=v1%(n+1);
        if((i0==0&&i1==0)||(i0==n&&i1==n))fixed[3*nv+e]=1;
    }
    std::vector<int> column(ndof,-1); int reduced_dimension=0;
    for(int d=0;d<ndof;++d)if(!fixed[d]&&!dependent[d])column[d]=reduced_dimension++;
    std::vector<Eigen::Triplet<double>> projection_trips;
    for(int d=0;d<ndof;++d)if(column[d]>=0)projection_trips.emplace_back(d,column[d],1.0);
    auto add=[&](int row,int source,double value){if(column[source]>=0)
        projection_trips.emplace_back(row,column[source],value);};
    for(int j=0;j<=n;++j){
        for(int c:{0,2})if(!fixed[pdof(0,j,c)]){
            add(pdof(0,j,c),pdof(1,j,c),4.0/3.0);
            add(pdof(0,j,c),pdof(2,j,c),-1.0/3.0);
        }
        for(int c:{1,2})if(!fixed[pdof(n,j,c)]){
            add(pdof(n,j,c),pdof(n-1,j,c),4.0/3.0);
            add(pdof(n,j,c),pdof(n-2,j,c),-1.0/3.0);
        }
    }
    Eigen::SparseMatrix<double> projection(ndof,reduced_dimension);
    projection.setFromTriplets(projection_trips.begin(),projection_trips.end());
    projection.makeCompressed();
    const Eigen::SparseMatrix<double> matrix=projection.transpose()*full*projection;
    Eigen::VectorXd external=Eigen::VectorXd::Zero(ndof);
    external[pdof(0,n,0)]=1.0; external[pdof(n,n,1)]=-1.0;
    const Eigen::VectorXd force=projection.transpose()*external;
    const Eigen::VectorXd reduced=solve_spd(matrix,force);
    const Eigen::VectorXd displacement=projection*reduced;
    const double observable=-displacement[pdof(n,n,1)];
    std::cout<<std::setprecision(17)
      <<"{\"benchmark\":\"M3\",\"method\":\"MidedgeTan\",\"n\":"<<n
      <<",\"raw_scalar_states\":"<<ndof<<",\"active_states\":"<<reduced_dimension
      <<",\"tangent_nnz\":"<<matrix.nonZeros()<<",\"observable\":"<<observable<<"}\n";
    return 0;
}
} // namespace

int main(int argc,char** argv){
    if(argc!=5 && argc!=6){std::cerr<<"usage: DRIVER benchmark n vertices weights [thickness]\n";return 2;}
    try{
        const std::string benchmark(argv[1]); const int n=std::atoi(argv[2]);
        if(benchmark=="plate")return plate(n,argv[3],argv[4],argc==6?std::atof(argv[5]):0.01);
        if(benchmark=="roof")return roof(n,argv[3],argv[4],argc==6?std::atof(argv[5]):0.25);
        if(benchmark=="M3")return m3(n,argv[3],argc==6?std::atof(argv[5]):0.04);
        return 2;
    }catch(const std::exception& error){std::cerr<<error.what()<<"\n";return 3;}
}
