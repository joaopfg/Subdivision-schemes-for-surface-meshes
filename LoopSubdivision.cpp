/* ========================================================================= *
 *                                                                           *
 *                       Luca Castelli Aleardi                       		 *
 *           Copyright (c) 2019, Ecole Polytechnique                		 *
 *           Department of Computer Science                  				 *
 *                          All rights reserved.                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of the course material developed for		             *
 *   INF574 Digital Representation and Analysis of Shapes (2019/20)			 *
 * ========================================================================= */
#include <bits/stdc++.h>
#include <igl/opengl/glfw/Viewer.h>

#ifndef HALFEDGE_DS_HEADER
  #define HALFEDGE_DS_HEADER
  #include "HalfedgeDS.cpp"
#endif

using namespace Eigen;
using namespace std;

typedef long long int ll;

const int mod2 = 1e9 + 9;
/**
 * @author Luca Castelli Aleardi (2019)
 */

struct Triangle2
{
  int a, b, c;

  Triangle2(const int &v1, const int &v2, const int &v3)
  {
    a = v1;
    b = v2;
    c = v3;
  }

  bool operator==(const Triangle2 &t) const
  {
    if (a == t.a && b == t.b && c == t.c)
      return true;
    else if (a == t.a && b == t.c && c == t.b)
      return true;
    else if (a == t.b && b == t.a && c == t.c)
      return true;
    else if (a == t.b && b == t.c && c == t.a)
      return true;
    else if (a == t.c && b == t.a && c == t.b)
      return true;
    else if (a == t.c && b == t.b && c == t.a)
      return true;
    else
      return false;
  }
};

class Triangle2HashFunction
{
public:
  size_t operator()(const Triangle2 &t) const
  {
    return (int)(((ll)t.a * (ll)t.b * (ll)t.c)%mod2);
  }
};

class LoopSubdivision
{

public:
    /** 
	 * Initialize the data structures
	 **/
    LoopSubdivision(MatrixXd &V_original, MatrixXi &F_original, HalfedgeDS &mesh)
    {
        he = &mesh;
        V = &V_original;
		F = &F_original; // NOT NEEDED if using the half-edge data structure
        int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
        int n = V_original.rows();         // number of vertices in the original mesh
        int f = he->sizeOfFaces();

        nVertices = n+e;
        nFaces = 4*f;
        
        // COMPLETED (initialize arrays V1 and F1)
        V1 = MatrixXd::Zero(n+e, 3);
        F1 = MatrixXi::Zero(4*f, 3);
    }

    /** 
	 * Perform the subdivision of the mesh (just perform one round subdivision). <b>
	 * As result, the resulting subdivided mesh is stored in arrays 'V1' and 'F1' 
	 **/
    void subdivide()
    {
        std::cout << "Performing one round subdivision" << endl;
        int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
        int n = he->sizeOfVertices();      // number of vertices in the original mesh
        int f = he->sizeOfFaces();         // number of vertices in the original mesh

        // COMPLETED
        // first step: perform a copy of the original points

        for(int i=0;i<n;i++){
            MatrixXd vnew = updateOriginalPoint(i);

            V1(i,0) = vnew(0,0);
            V1(i,1) = vnew(0,1);
            V1(i,2) = vnew(0,2);
        }

        // second step: compute new midpoint vertices and assign a number, between 0..e-1, to all halfedges
        int dirToInd[e], indToDir[2*e];
        bool visitEdge[2*e];

        memset(visitEdge, false, sizeof(visitEdge));

        int cont = 0;

        for(int i=0;i<2*e;i++){
            if(!visitEdge[i]){
                int iOpp = he->getOpposite(i);
                visitEdge[i] = true;
                visitEdge[iOpp] = true;
                dirToInd[cont] = i;
                indToDir[i] = cont;
                indToDir[iOpp] = cont;
                cont++;
            }
        }

        for(int i=n;i<n+e;i++){
            MatrixXd Pmid = computeEdgePoint(dirToInd[i-n]);
            V1(i,0) = Pmid(0,0);
            V1(i,1) = Pmid(0,1);
            V1(i,2) = Pmid(0,2);
        }

        // third step: set the face/vertex incidence relations
        cont = 0;
        unordered_map<Triangle2, bool, Triangle2HashFunction> visit;
        unordered_map<Triangle2, bool, Triangle2HashFunction>::iterator it;

        for(int i=0;i<e;i++){
            int P4P1 = dirToInd[i];
            int Pc = n + indToDir[P4P1];
            int P1 = he->getTarget(P4P1);
            int P5P4 = he->getPrev(P4P1);
            int P5 = n + indToDir[P5P4];
            int P4 = he->getTarget(P5P4);
            int P1P6 = he->getPrev(P5P4);
            int P6 = n + indToDir[P1P6];
            int P1P4 = he->getOpposite(P4P1);
            int P2P1 = he->getPrev(P1P4);
            int P2 = n + indToDir[P2P1];
            int P4P3 = he->getNext(P1P4);
            int P3 = n + indToDir[P4P3];

            Triangle2 P2P1Pc = Triangle2(P2, P1, Pc);
            Triangle2 P1P6Pc = Triangle2(P1, P6, Pc);
            Triangle2 P4P3Pc = Triangle2(P4, P3, Pc);
            Triangle2 P5P4Pc = Triangle2(P5, P4, Pc);
            Triangle2 P3P2Pc = Triangle2(P3, P2, Pc);
            Triangle2 P6P5Pc = Triangle2(P6, P5, Pc);

            if(visit.find(P2P1Pc) == visit.end()){
                visit[P2P1Pc] = true;
                F1(cont,0) = P2;
                F1(cont,1) = P1;
                F1(cont,2) = Pc;
                cont++;

                if(cont == 4*f) break;
            }

            if(visit.find(P1P6Pc) == visit.end()){
                visit[P1P6Pc] = true;
                F1(cont,0) = P1;
                F1(cont,1) = P6;
                F1(cont,2) = Pc;
                cont++;

                if(cont == 4*f) break;
            }

            if(visit.find(P4P3Pc) == visit.end()){
                visit[P4P3Pc] = true;
                F1(cont,0) = P4;
                F1(cont,1) = P3;
                F1(cont,2) = Pc;
                cont++;

                if(cont == 4*f) break;
            }

            if(visit.find(P5P4Pc) == visit.end()){
                visit[P5P4Pc] = true;
                F1(cont,0) = P5;
                F1(cont,1) = P4;
                F1(cont,2) = Pc;
                cont++;

                if(cont == 4*f) break;
            }

            if(visit.find(P3P2Pc) == visit.end()){
                visit[P3P2Pc] = true;
                F1(cont,0) = P3;
                F1(cont,1) = P2;
                F1(cont,2) = Pc;
                cont++;

                if(cont == 4*f) break;
            }

            if(visit.find(P6P5Pc) == visit.end()){
                visit[P6P5Pc] = true;
                F1(cont,0) = P6;
                F1(cont,1) = P5;
                F1(cont,2) = Pc;
                cont++;

                if(cont == 4*f) break;
            }
        }
    }

    /** 
	 * Return the number of half-edges
	 **/
    MatrixXd getVertexCoordinates()
    {
        return V1;
    }

    /** 
	 * Return the number of faces
	 **/
    MatrixXi getFaces()
    {
        return F1;
    }

    /** 
	 * Print the combinatorial information of the subdivided mesh <b>
	 * verbosity=0: print only the number of vertices and faces <b>
	 * verbosity=1: print all incidence relations
	 **/
    void print(int verbosity)
    {
        cout << "\tn=" << nVertices << ", f=" << nFaces << endl;

        if (verbosity > 0) // print all vertex coordinates and face/vertex incidence relations
        {
            for (int i = 0; i < nVertices; i++)
            {
                cout << "v" << i << ": " << V1.row(i) << endl;
            }

            std::cout << "new faces: " << nFaces << endl;
            for (int i = 0; i < nFaces; i++)
            {
                cout << "f" << i << ": " << F1.row(i) << endl;
            }
        }
    }

private:
    /**
	 * Compute the midpoint of the given half-edge 'h=(u,v)'
	 */
    MatrixXd computeEdgePoint(int h)
    {
        // COMPLETED
        int e2Ind = he->getTarget(h);
        int t2e1 = he->getPrev(h);
        int e1Ind = he->getTarget(t2e1);
        int e2t2 = he->getPrev(t2e1);
        int t2Ind = he->getTarget(e2t2);
        int hOpp = he->getOpposite(h);
        int e1t1 = he->getNext(hOpp);
        int t1Ind = he->getTarget(e1t1);

        MatrixXd e1(1,3), e2(1,3), t1(1,3), t2(1,3);

        e1(0,0) = (*V)(e1Ind, 0);
        e1(0,1) = (*V)(e1Ind, 1);
        e1(0,2) = (*V)(e1Ind, 2);

        e2(0,0) = (*V)(e2Ind, 0);
        e2(0,1) = (*V)(e2Ind, 1);
        e2(0,2) = (*V)(e2Ind, 2);

        t1(0,0) = (*V)(t1Ind, 0);
        t1(0,1) = (*V)(t1Ind, 1);
        t1(0,2) = (*V)(t1Ind, 2);

        t2(0,0) = (*V)(t2Ind, 0);
        t2(0,1) = (*V)(t2Ind, 1);
        t2(0,2) = (*V)(t2Ind, 2);

        MatrixXd ans(1,3);

        ans(0,0) = 3*e1(0,0)/8 + 3*e2(0,0)/8 + t1(0,0)/8 + t2(0,0)/8;  
        ans(0,1) = 3*e1(0,1)/8 + 3*e2(0,1)/8 + t1(0,1)/8 + t2(0,1)/8;
        ans(0,2) = 3*e1(0,2)/8 + 3*e2(0,2)/8 + t1(0,2)/8 + t2(0,2)/8;

        return ans;
    }

    /**
	 * Given a vertex 'v' of the original mesh, compute and return its new coordinates
	 */
    MatrixXd updateOriginalPoint(int v)
    {
        // COMPLETED
        int d = vertexDegree(v);
        float alpha;

        if(d == 3) alpha = 3.0/16;
        else alpha = 3.0/(8*d);

        MatrixXd vnew(1,3);

        vnew(0,0) = (1.0 - alpha*d)*(*V)(v,0);
        vnew(0,1) = (1.0 - alpha*d)*(*V)(v,1);
        vnew(0,2) = (1.0 - alpha*d)*(*V)(v,2);

        int eCrawl = he->getEdge(v), vCrawl;

        while(true){
            eCrawl = he->getOpposite(eCrawl);
            vCrawl = he->getTarget(eCrawl);

            vnew(0,0) += alpha*(*V)(vCrawl, 0);
            vnew(0,1) += alpha*(*V)(vCrawl, 1);
            vnew(0,2) += alpha*(*V)(vCrawl, 2);

            eCrawl = he->getPrev(eCrawl);

            if(eCrawl == he->getEdge(v)) break;
        }

        return vnew;
    }

    int vertexDegree(int v)
    {
        // COMPLETED
        int deg = 0, eCrawl = he->getEdge(v);

        while(true){
            deg++;
            eCrawl = he->getOpposite(eCrawl);
            eCrawl = he->getPrev(eCrawl);

            if(eCrawl == he->getEdge(v)) break;
        }

        return deg;
    }

    /** Half-edge representation of the original input mesh */
    HalfedgeDS *he;
    MatrixXd *V; // vertex coordinates of the original input mesh

	/** faces/vertex incidence relations in the original mesh */
	MatrixXi *F; // REMARK: not needed if using the half-edge data structure

    int nVertices, nFaces; // number of vertices, faces in the new subdivided mesh
    MatrixXd V1;          // vertex coordinates of the new subdivided mesh
    MatrixXi F1;          // faces of the new subdivided mesh
};
