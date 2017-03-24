// This file is part of ugmat, a program for analysing and comparing matrices
//
// Copyright (C) 2017 Sebastian Reiter, G-CSC Frankfurt <sreiter@gcsc.uni-frankfurt.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>
#include <cstring>
#include <map>

using namespace std;

typedef double	number;
typedef unsigned int uint;


namespace params
{
	bool	makeAbs			= false;
	bool	verbose			= false;
}


class CommonError{};
#define CHECK(expr, msg) 	{if(!(expr)) {cout << "ERROR: " << msg << endl; throw CommonError();}}

struct Node{
	Node() : x(0), y(0), z(0), ci(0)	{}
	
	bool operator == (const Node& p) const
	{
		return x == p.x && y == p.y && z == p.z && ci == p.ci;
	}

	bool operator < (const Node& p) const
	{
		if(ci < p.ci)
			return true;
		else if(ci > p.ci)
			return false;
		
		if(x < p.x)
			return true;
		else if(x > p.x)
			return false;

		if(y < p.y)
			return true;
		else if(y > p.y)
			return false;

		if(z < p.z)
			return true;
		else if(z > p.z)
			return false;

		return false;
	}
	
	number x, y, z;
	int ci;			///< component index
};

ostream& operator << (ostream& out, const Node& p)
{
	out << "(" << p.x << ", " << p.y << ", " << p.z << ")[" << p.ci << "]";
	return out;
}


typedef pair<size_t, number>	Link;

struct Matrix{
	Matrix() : worldDim(0)	{}
	
	int	worldDim;
	vector<Node> 			nodes;
	vector<vector<Link > >	links;
};


struct SimpleVector{
	SimpleVector() : worldDim(0)	{}

	int worldDim;
	vector<pair<Node, number> >	entries;
};


void MakeAbs(Matrix& m)
{
	for(size_t inode = 0; inode < m.nodes.size(); ++inode){
		vector<Link>& links = m.links[inode];
		for(size_t ilink = 0; ilink < links.size(); ++ilink){
			links[ilink].second = fabs(links[ilink].second);
		}
	}
}


bool LoadMatrix(Matrix& m, const char* filename)
{
	cout << "INFO -- loading matrix from " << filename << endl;
	
	string line;
	ifstream in(filename);
	if(!in){
		cout << "ERROR -- File not found: " << filename << endl;
		return false;
	}
	
	int blockSize;
	in >> blockSize;
	in >> m.worldDim;
	
	int numEntries;
	in >> numEntries;
	m.nodes.clear();
	m.nodes.reserve(numEntries);
	
	std::map<Node, size_t> nodeMap;


	for(int i = 0; i < numEntries; ++i){
		Node p;
		switch(m.worldDim){
			case 1:	in >> p.x; break;
			case 2:	in >> p.x >> p.y; break;
			case 3:	in >> p.x >> p.y >> p.z; break;
			default:
				cout << "ERROR -- Unsupported world-dimension (" << m.worldDim
					 << ") during write: " << filename << endl;
				return false;
		}

	//	during the lookup p.ci should always be 0
		size_t& ci = nodeMap[p];
		p.ci = ci;
		++ci;
		m.nodes.push_back(p);
	}
	
//	get the rest of the last line
	getline(in, line);

//	read some arbitrary value which separates nodes and links
	getline(in, line);

//	read data values
	m.links.clear();
	m.links.resize(numEntries);

	size_t numNANs = 0;
	vector<double> values;

	while(!in.eof()){
		getline(in, line);
		if(line.empty())
			continue;
		size_t start = 0;
		values.clear();
		while(start < line.size()) {
			const char cur = line[start];
			if(cur == ' '){
				++start;
				continue;
			}
			else if (cur == '[' || cur == ']'){
				++start;
				continue;
			}
			else{
				size_t end = line.find(' ', start);
				if(end == string::npos)
					end = line.size();
				const size_t num = end - start;
				if((line.find("n", start, num) != string::npos)
					|| (line.find("N", start, num) != string::npos))
				{
					++numNANs;
					continue;
				}

				char* endPtr;
				const char* startPtr = line.c_str() + start;
				values.push_back(strtod(startPtr, &endPtr));
				start += endPtr - startPtr;
			}
		}

		if(values.size() < 3){
			cout << "ERROR -- Not enough values specified in connection. In File: "
				<< filename << endl;
			cout << "line read: " << line << endl;
			continue;
		}
		else if(values.size() > 3){
			cout << "ERROR -- Too many values specified in connection. Block matrices currently not supported. In File: "
				<< filename << endl;
			cout << "line read: " << line << endl;
			continue;
		}

		const size_t ind1 = static_cast<size_t>(values[0]);
		const size_t ind2 = static_cast<size_t>(values[1]);

		if((ind1 < 0) || (ind1 >= numEntries)){
			cout << "ERROR -- Bad source index: " << ind1 << ". In File: " << filename << endl;
		}

		if((ind2 < 0) || (ind2 >= numEntries)){
			cout << "ERROR -- Bad target index: " << ind2 << ". In File: " << filename << endl;
		}

		m.links[ind1].push_back(Link(ind2, values[2]));
	}

	if(numNANs > 0)
		cout << "  -> WARNING: vector contains " << numNANs << " 'nan' entries!" << endl;
	
	if(params::makeAbs)
		MakeAbs(m);

	return true;
}

bool SaveSimpleVector(const SimpleVector& v, const char* filename)
{
	cout << "INFO -- saving vector to " << filename << endl;
	
	ofstream out(filename);
	if(!out){
		cout << "ERROR -- File can not be opened for write: " << filename << endl;
		return false;
	}
	
	out << int(1) << endl;
	out << v.worldDim << endl;
	out << v.entries.size() << endl;
	
	for(size_t i = 0; i < v.entries.size(); ++i){
		const Node& n = v.entries[i].first;
		switch(v.worldDim){
			case 1:	out << n.x << endl; break;
			case 2:	out << n.x << " " << n.y << endl; break;
			case 3:	out << n.x << " " << n.y << " " << n.z << endl; break;
			default:
				cout << "ERROR -- Unsupported world-dimension (" << v.worldDim
					 << ") during write: " << filename << endl;
				return false;
		}
	}
	
	out << int(1) << endl;
	
	for(size_t i = 0; i < v.entries.size(); ++i){
		out << int(i) << " " << int(i) << " "
			<< setprecision(numeric_limits<number>::digits10 + 1)
			<< v.entries[i].second << endl;
	}
	
	return true;
}

void Diagonal(SimpleVector& vOut, const Matrix& m)
{
	vOut.worldDim = m.worldDim;
	vOut.entries.resize(m.nodes.size());
	for(size_t inode = 0; inode < m.nodes.size(); ++inode){
		vOut.entries[inode].first = m.nodes[inode];
		vOut.entries[inode].second = 0;
		
		const vector<Link>& links = m.links[inode];
		for(size_t ilink = 0; ilink < links.size(); ++ilink){
			if(links[ilink].first == inode){
				vOut.entries[inode].second = links[ilink].second;
				break;
			}
		}
	}
}

void MaxOffDiag(SimpleVector& vOut, const Matrix& m)
{
	vOut.worldDim = m.worldDim;
	vOut.entries.resize(m.nodes.size());
	for(size_t inode = 0; inode < m.nodes.size(); ++inode){
		vOut.entries[inode].first = m.nodes[inode];
		bool gotOffDiag = false;

		const vector<Link>& links = m.links[inode];
		for(size_t ilink = 0; ilink < links.size(); ++ilink){
			vOut.entries[inode].second = 0;
			if(links[ilink].first != inode){
				vOut.entries[inode].second = links[ilink].second;
				break;
			}
		}

		for(size_t ilink = 0; ilink < links.size(); ++ilink){
			if(links[ilink].first != inode){
				vOut.entries[inode].second
					= max(vOut.entries[inode].second, links[ilink].second);
			}
		}
	}
}

void MinOffDiag(SimpleVector& vOut, const Matrix& m)
{
	vOut.worldDim = m.worldDim;
	vOut.entries.resize(m.nodes.size());
	for(size_t inode = 0; inode < m.nodes.size(); ++inode){
		vOut.entries[inode].first = m.nodes[inode];
		bool gotOffDiag = false;

		const vector<Link>& links = m.links[inode];
		for(size_t ilink = 0; ilink < links.size(); ++ilink){
			vOut.entries[inode].second = 0;
			if(links[ilink].first != inode){
				vOut.entries[inode].second = links[ilink].second;
				break;
			}
		}

		for(size_t ilink = 0; ilink < links.size(); ++ilink){
			if(links[ilink].first != inode){
				vOut.entries[inode].second
					= min(vOut.entries[inode].second, links[ilink].second);
			}
		}
	}
}


int main(int argc, char** argv)
{
	int		defHistoSecs	= 5;

	static const int maxNumFiles = 3;
	const char* file[maxNumFiles];
	int numFiles = 0;

	for(int i = 2; i < argc; ++i){
		if(argv[i][0] == '-'){

			if(strcmp(argv[i], "-abs") == 0){
				params::makeAbs = true;
			}

			else if(strcmp(argv[i], "-verbose") == 0){
				params::verbose = true;
			}

			else{
				cout << "Invalid option supplied: " << argv[i] << endl;
				return 1;
			}
		}

		else if(numFiles < maxNumFiles){
			file[numFiles] = argv[i];
			++numFiles;
		}

		else{
			cout << "Can't interpret parameter " << argv[i] << ": Too many parameters specified." << endl;
			return 1;
		}
	}
	
	string command;
	if(argc > 1)
		command = argv[1];


	try{
		if(command.find("diag") == 0){
			CHECK(numFiles == 2, "An in-file and an out-file have to be specified");
			Matrix m;
			LoadMatrix(m, file[0]);
			SimpleVector v;
			Diagonal(v, m);
			SaveSimpleVector(v, file[1]);
		}
		
		else if(command.find("maxOffDiag") == 0){
			CHECK(numFiles == 2, "An in-file and an out-file have to be specified");
			Matrix m;
			LoadMatrix(m, file[0]);
			SimpleVector v;
			MaxOffDiag(v, m);
			SaveSimpleVector(v, file[1]);
		}
		
		else if(command.find("minOffDiag") == 0){
			CHECK(numFiles == 2, "An in-file and an out-file have to be specified");
			Matrix m;
			LoadMatrix(m, file[0]);
			SimpleVector v;
			MinOffDiag(v, m);
			SaveSimpleVector(v, file[1]);
		}

		else{
			cout << endl;
			cout << "ugmat - (c) 2017 Sebastian Reiter, G-CSC Frankfurt" << endl;
			cout << endl;
			cout << "USAGE: ugmat command [options] [files]" << endl;
			cout << "OR:    ugmat command [files] [options]" << endl << endl;

			cout << "SAMPLE: ugmat diag -consistent A.mat diag.vec" << endl << endl;

			cout << "COMMANDS:" << endl;

			cout << "  diag:       Extracts the diagonal of the specified matrix and writes it to the specified" << endl;
			cout << "              vector file." << endl;
  			cout << "              2 Files required - 1: in-file ('.mat'), 2: out-file ('.ugx')" << endl << endl;

  			cout << "  maxOffDiag: Extracts the maximal off-diagonal values of the specified matrix" << endl;
			cout << "              and writes them to the specified vector file." << endl;
  			cout << "              2 Files required - 1: in-file ('.mat'), 2: out-file ('.ugx')" << endl << endl;

  			cout << "  minOffDiag: Extracts the minimal off-diagonal values of the specified matrix" << endl;
			cout << "              and writes them to the specified vector file." << endl;
  			cout << "              2 Files required - 1: in-file ('.mat'), 2: out-file ('.ugx')" << endl << endl;

			cout << "OPTIONS:" << endl;
			cout << "  -abs:             Only absolute values will be regarded." << endl << endl;

			cout << "  -verbose:         If specified, additional information is printed for each processed vector." << endl << endl;
		}
	}
	catch(...){
		return 1;
	}
	return 0;
}
