/* CS3242 3D Modeling and Animation
 * Programming Assignment I
 * School of Computing
 * National University of Singapore
 */

#include "Mesh.h"
#include <map>
#include <queue>
#include <cmath>
#include <Eigen\Dense>
#include <iostream>

Mesh::Mesh(const char* filename){	
	loadMF(filename);
}

struct Comp
{
	bool operator()(const WEEdge& a, const WEEdge& b)
	{
		return a.cost > b.cost;
	}
};

std::map<int, WEVertex> wevMap;
std::map<std::pair<int, int>, WEEdge> weeMap;
std::priority_queue<WEEdge, std::vector<WEEdge>, Comp> pq;

void Mesh::loadMF(const char* filename){
	if(V.size()>0) V.clear();
	if(F.size()>0) F.clear();
	std::ifstream infile;
	infile.open(filename, std::ios::in);
	std::string strbuff;
	while(std::getline(infile,strbuff)){
		std::stringstream ss;
		ss<<strbuff;
		char type;
		ss>>type;
		if(type=='v'){
			Vertex v;
			ss>>v.x>>v.y>>v.z;
			V.push_back(v);
		}
		else if(type=='f'){
			Face f;
			ss>>f.a>>f.b>>f.c;
			F.push_back(f);
		}
	}
	infile.close();
}

void Mesh::writeMF(const char* filename){
	std::ofstream outfile;
	outfile.open(filename, std::ios::out);
	std::string strbuff;
	for(uint i=0;i<V.size();i++){
		outfile<<"v "<<V[i].x<<" "<<V[i].y<<" "<<V[i].z<<std::endl;
	}
	for(uint i=0;i<F.size();i++){
		outfile<<"f "<<F[i].a<<" "<<F[i].b<<" "<<F[i].c<<std::endl;
	}
	outfile.close();
}

void Mesh::simplifyMesh(const char* input, const char* output, int faceCnt){
	// you may assume inputs are always valid
	loadMF(input);
	convertMesh();
	std::cout<<"Original face count: "<<F.size()<<std::endl;
	std::cout << "Original vertices count: " << V.size() << std::endl;
	std::cout << "Vertex map count: " << wevMap.size() << std::endl;
	std::cout << "Edge map count: " << weeMap.size() << std::endl;
	// do mesh simplification
	// write your code here
	// sort winged edge list into ascending cost
	for (auto const &ent : weeMap)
	{
		WEEdge focus = ent.second;
		uint v1id = ent.first.first;
		uint v2id = ent.first.second;
		WEVertex v1, v2;
		v1 = wevMap[v1id];
		v2 = wevMap[v2id];
		// melax algo calculation
		float diff = sqrt(pow((v1.v.x - v2.v.x), 2) + pow((v1.v.y - v2.v.y), 2) + pow((v1.v.z - v2.v.z), 2));
		// faces shareed by v1, v2 = edge next and prev faces
		Face shared1, shared2;
		shared1 = focus.pf;
		shared2 = focus.nf;
		float max = -99999;
		for (size_t i = 0; i < v1.nearbyFaces.size(); i++)
		{
			Face f = v1.nearbyFaces.at(i);
			Eigen::Vector3f a(V.at(f.a-1).x, V.at(f.a-1).y, V.at(f.a - 1).z);
			Eigen::Vector3f b(V.at(f.b - 1).x, V.at(f.b - 1).y, V.at(f.b - 1).z);
			Eigen::Vector3f c(V.at(f.c - 1).x, V.at(f.c - 1).y, V.at(f.c - 1).z);
			float minF;
			Eigen::Vector3f fNormal = ((b - a).cross(c - a));
			// calculate normal for shared face 1
			Eigen::Vector3f s1a(V.at(shared1.a-1).x, V.at(shared1.a-1).y, V.at(shared1.a-1).z);
			Eigen::Vector3f s1b(V.at(shared1.b-1).x, V.at(shared1.b-1).y, V.at(shared1.b-1).z);
			Eigen::Vector3f s1c(V.at(shared1.c-1).x, V.at(shared1.c-1).y, V.at(shared1.c-1).z);
			Eigen::Vector3f s1Normal = ((s1b - s1a).cross(s1c - s1a));
			// calculate normal for shared face 2
			Eigen::Vector3f s2a(V.at(shared2.a - 1).x, V.at(shared2.a - 1).y, V.at(shared2.a - 1).z);
			Eigen::Vector3f s2b(V.at(shared2.b - 1).x, V.at(shared2.b - 1).y, V.at(shared2.b - 1).z);
			Eigen::Vector3f s2c(V.at(shared2.c - 1).x, V.at(shared2.c - 1).y, V.at(shared2.c - 1).z);
			Eigen::Vector3f s2Normal = ((s2b - s2a).cross(s2c - s2a));
			minF = std::min((1 - fNormal.dot(s1Normal)) / 2, (1 - fNormal.dot(s2Normal)) / 2);
			if (max < minF * diff)
			{
				max = minF * diff;
			}
		}
		focus.cost = max;
		pq.push(focus);
	}
	// collapse edges till faceCnt reached
	for (int i = faceCnt; i < (int)F.size(); i++)
	{
		WEEdge topEdge = pq.top();
		Edge focusEdge = topEdge.e;
		Vertex finalV;
		finalV.x = 0.5*(V.at(focusEdge.a).x + V.at(focusEdge.b).x);
		finalV.y = 0.5*(V.at(focusEdge.a).y + V.at(focusEdge.b).y);
		finalV.z = 0.5*(V.at(focusEdge.a).z + V.at(focusEdge.b).z);
		V.at(focusEdge.a) = finalV;
		V.at(focusEdge.b) = finalV;
		pq.pop();
	}

	revertMesh();
	writeMF(output);
}

void Mesh::convertMesh()
{
	// write your code here
	uint a, b, c;;
	for (size_t i = 0; i < F.size(); i++)
	{
		Face current = F.at(i);
		WEVertex v1, v2, v3;
		a = current.a-1;
		b = current.b-1;
		c = current.c-1;
		// make sure a>b>c for easier search
		if (c > b)
		{
			int temp = b;
			b = c;
			c = temp;
		}
		if (b > a) {
			int temp = a;
			a = b;
			b = temp;
			if (c > b)
			{
				int temp = b;
				b = c;
				c = temp;
			}
		}
		v1.v = V.at(a);
		v2.v = V.at(b);
		v3.v = V.at(c);

		// check if the vertex is stored and subsequently check if the edge is stored
		WEEdge matchAB, matchBC, matchAC;
		bool AB, BC, AC;
		AB = false;
		BC = false;
		AC = false;
		Edge edgeAB, edgeBC, edgeAC;
		// init edges
		edgeAB.a = 0;
		edgeAB.b = 0;
		edgeBC.a = 0;
		edgeBC.b = 0;
		edgeAC.a = 0;
		edgeAC.b = 0;
		std::map<int, WEVertex>::iterator itV;
		if ((itV = wevMap.find(a)) != wevMap.end())
		{
			v1 = itV->second;
			// check if a, b connects
			if (std::find(v1.nearbyVertices.begin(), v1.nearbyVertices.end(), b) != v1.nearbyVertices.end())
			{
				// found edge connecting a, b
				matchAB = weeMap[std::make_pair(a, b)];
				edgeAB = matchAB.e;
				AB = true;
			}
			else
			{
				// b is not in a's list so update
				v1.nearbyVertices.push_back(b);
				// check if b is already stored
				if ((itV = wevMap.find(b)) != wevMap.end())
				{
					v2 = itV->second;
					// update b's list too
					v2.nearbyVertices.push_back(a);
					v2.nearbyFaces.push_back(current);
					wevMap[b] = v2;
				}
			}
			// check if a, c connects
			if (std::find(v1.nearbyVertices.begin(), v1.nearbyVertices.end(), c) != v1.nearbyVertices.end())
			{
				// found edge connecting a, c
				matchAC = weeMap[std::make_pair(a, c)];
				edgeAC = matchAC.e;
				AC = true;
			}
			else
			{
				// c is not in a's list so update
				v1.nearbyVertices.push_back(c);
				// check if b is already stored
				if ((itV = wevMap.find(c)) != wevMap.end())
				{
					v3 = itV->second;
					// update b's list too
					v3.nearbyVertices.push_back(a);
					v3.nearbyFaces.push_back(current);
					wevMap[c] = v3;
				}
			}
			// update v1
			v1.nearbyFaces.push_back(current);
			wevMap[a] = v1;
		}
		else
		{
			// a is not found in wevMap so add a in, edge ab and ac are not stored
			v1.nearbyVertices.push_back(b);
			v1.nearbyVertices.push_back(c);
			v1.nearbyFaces.push_back(current);
			wevMap[a] = v1;
			// check if b is stored
			if ((itV = wevMap.find(b)) != wevMap.end())
			{
				v2 = itV->second;
				// since a is new, add a to b's list
				v2.nearbyVertices.push_back(a);
				// b is found, check if b, c connects
				if (std::find(v2.nearbyVertices.begin(), v2.nearbyVertices.end(), c) != v2.nearbyVertices.end())
				{
					// found edge connecting b, c
					matchBC = weeMap[std::make_pair(b, c)];
					edgeBC = matchBC.e;
					BC = true;
				}
				else
				{
					// c is not in b's list so update
					v2.nearbyVertices.push_back(c);
					// check if c is already stored
					if ((itV = wevMap.find(c)) != wevMap.end())
					{
						v3 = itV->second;
						// update c's list too
						v3.nearbyVertices.push_back(b);
						v3.nearbyVertices.push_back(a);
					}
					else
					{
						// c is new
						v3.nearbyVertices.push_back(b);
						v3.nearbyVertices.push_back(a);
						v3.nearbyFaces.push_back(current);
						wevMap[c] = v3;
					}
				}
				v2.nearbyFaces.push_back(current);
				wevMap[b] = v2;
			}
			else
			{
				// b is not found in wevMap so add b in
				v2.nearbyVertices.push_back(a);
				v2.nearbyVertices.push_back(c);
				v2.nearbyFaces.push_back(current);
				wevMap[b] = v2;
				// check if c is in wevMap
				if ((itV = wevMap.find(c)) != wevMap.end())
				{
					v3 = itV->second;
					v3.nearbyVertices.push_back(b);
					v3.nearbyVertices.push_back(a);
					v3.nearbyFaces.push_back(current);
					wevMap[c] = v3;
				}
				else
				{
					// c is not found in mevMap so add c in
					v3.nearbyVertices.push_back(b);
					v3.nearbyVertices.push_back(a);
					v3.nearbyFaces.push_back(current);
					wevMap[c] = v3;
				}
			}
		}

		std::vector<int> testing;
		testing.push_back(a);
		// empty edge
		Edge placeholder;
		placeholder.a = 1;
		placeholder.b = 1;
		Face placeholderF;
		placeholderF.a = 1;
		placeholderF.b = 1;
		placeholderF.c = 1;
		// update edges
		if (!AB) {
			edgeAB.a = a;
			edgeAB.b = b;
		}
		if (!BC) {
			edgeBC.a = b;
			edgeBC.b = c;
		}
		if (!AC) {
			edgeAC.a = a;
			edgeAC.b = c;
		}
		// update WEEdge
		if (!AB) {
			matchAB.e = edgeAB;
			matchAB.ncw = edgeBC;
			matchAB.nccw = edgeAC;
			matchAB.nf = current;
			matchAB.pcw = placeholder;
			matchAB.pccw = placeholder;
			matchAB.pf = placeholderF;
		} 
		else
		{
			matchAB.pcw = edgeBC;
			matchAB.pccw = edgeAC;
			matchAB.pf = current;
		}
		if (!BC) {
			matchBC.e = edgeBC;
			matchBC.ncw = edgeAC;
			matchBC.nccw = edgeAB;
			matchBC.nf = current;
			matchBC.pcw = placeholder;
			matchBC.pccw = placeholder;
			matchBC.pf = placeholderF;
		}
		else
		{
			matchBC.pcw = edgeAC;
			matchBC.pccw = edgeAB;
			matchBC.pf = current;
		}
		if (!AC) {
			matchAC.e = edgeAC;
			matchAC.ncw = edgeAB;
			matchAC.nccw = edgeBC;
			matchAC.nf = current;
			matchAC.pcw = placeholder;
			matchAC.pccw = placeholder;
			matchAC.pf = placeholderF;
		}
		else
		{
			matchAC.pcw = edgeAB;
			matchAC.pccw = edgeBC;
			matchAC.pf = current;
		}
		weeMap[std::make_pair(a, b)] = matchAB;
		weeMap[std::make_pair(b, c)] = matchBC;
		weeMap[std::make_pair(a, c)] = matchAC;
	}
}

void Mesh::revertMesh()
{
	// write your code here
	F.clear();
	for (auto const &ent : weeMap) {
		Face tempF = ent.second.nf;
		F.push_back(tempF);
	}
}

std::vector<WEVertex> Mesh::neighborVertices(WEVertex v)
{
	// write your code here

	return std::vector<WEVertex>();
}

std::vector<WEFace> Mesh::neighborFaces(WEVertex v)
{
	// write your code here
	
	return std::vector<WEFace>();
}

std::vector<WEVertex> Mesh::adjacentVertices(WEFace f)
{
	// write your code here
	
	return std::vector<WEVertex>();
}

std::vector<WEEdge> adjacentEdges(WEFace f)
{
	// write your code here
	
	return std::vector<WEEdge>();
}