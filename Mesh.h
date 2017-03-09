/* CS3242 3D Modeling and Animation
 * Programming Assignment I
 * School of Computing
 * National University of Singapore
 */
 
#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

typedef unsigned int uint;

//data structure of indexed face set
typedef struct {
	//3d coordinates
	float x;
	float y;
	float z;
} Vertex;

typedef struct{
	//three vertex ids
	uint a,b,c;
} Face;

//data structure of Winged Edge
typedef struct {
	uint a, b;
} Edge;

typedef struct {
	// write your code here
	Vertex v;
	Edge e;
	std::vector<int> nearbyVertices;
	std::vector<Face> nearbyFaces;
} WEVertex;

typedef struct {
	// write your code here
	Edge e, pcw, pccw, ncw, nccw;
	Face pf, nf;
	int cost;
} WEEdge;

typedef struct {
	// write your code here
	uint a, b, c;
	Edge e;
} WEFace;

class Mesh{
private:
	std::vector<Vertex> V;
	std::vector<Face> F;

	std::vector<WEVertex> WEV;
	std::vector<WEEdge> WEE;
	std::vector<WEFace> WEF;
public:
	Mesh() {};
	Mesh(const char*);
	//load a Mesh from .mesh file
	void loadMF(const char*);
	//write a Mesh to .mesh file (no header)
	void writeMF(const char*);
	//simplify a mesh
	void simplifyMesh(const char* input, const char* output, int faceCnt);
	//turn indexed face set to winged edge
	void convertMesh();
	//turn winged edge to indexed face set
	void revertMesh();

	//helper methods
	std::vector<WEVertex> neighborVertices(WEVertex v);
	std::vector<WEFace> neighborFaces(WEVertex v);
	std::vector<WEVertex> adjacentVertices(WEFace f);
	std::vector<WEEdge> adjacentEdges(WEFace f);
	
};
#endif