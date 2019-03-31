//using burning algorithm to get all details of structure

#include <bits/stdc++.h> 
#include <algorithm>
#include <assert.h>



using namespace std; 


struct Atominfo {
	int natoms;
	vector <int> atomnumber;
	vector <int> atomtype;
	vector <double> xcoords;
	vector <double> ycoords;
	vector <double> zcoords;
	vector <int> xi;
	vector <int> yi;
	vector <int> zi;
}; 

struct Bondinfo {
	int nbonds;
	vector <int> bondnumber;
	vector <int> bondtype;
	vector <int> bondfirst;
	vector <int> bondsecond;
};

struct Lmpsabinfo {
	Atominfo ai;
	Bondinfo bi;
	int atomtypes;
	int bondtypes; 
	vector <double> masses;
	
	vector <int>  molecule;
	double xboxlen;
	double yboxlen;
	double zboxlen;
	double xlo;
	double xhi;
	double ylo;
	double yhi;
	double zlo;
	double zhi;
};

struct P1P2list {
    vector <int> P1;
    vector <int> P2;
};

struct ClusterSize {
    int P1;
    int P2;
    bool bigenough = false;
    bool bigx = false;
    bool bigy = false;
    bool bigz = false;
};

Lmpsabinfo read_finalstructwithbondinfo(string filename, int nlines)
/* filename of lammps data file, nlines doesn't need to be exact but needs to be big enough to cover bonds, okay if larger than actual number of lines */
{
	ifstream infile;
	infile.open(filename.c_str());
	string line;
	char * pch;
	int linecomment = nlines;
	int linetypeinfo = nlines;
	int linebox = nlines;
	int linemasses = nlines;
	int linebondcoeff = nlines;
	int lineatom = nlines;
	int linebondsend = nlines;
	double  boxcoords [3][2];
	int natomlines = nlines;
	int nindex;
	int bindex;
	int natoms;
	int atomtypes;
	int nbonds;
	int bondtypes;
	int nbtype;
	Atominfo ai;
	Bondinfo bi;
	int nmass;
	Lmpsabinfo lmpi;
	for (int i=0; i < linebondsend; i++){
	//cout << "i " << i << endl;
	vector<string> vec_line;
	getline(infile,line);
	if (line.length() > 0){
	//cout << "line " << line << endl;
	//cout << arrayTokens << endl;
	char lineArr[line.length()];
	//cout << "lineArr initialized " << lineArr << endl;
	strcpy(lineArr, line.c_str());
 	//cout << "lineArr copied " << lineArr << endl;
		
	pch = strtok(lineArr, " ");
	while( pch != NULL ) {	vec_line.push_back(pch); pch = strtok( NULL, " " );}
	//cout << "words in line "  << vec_line.size() << endl;
	//for (int ii =0; ii < vec_line.size(); ii++){
	//cout << "word " << ii << " " << vec_line[ii] << endl; 
	//}
	if (i == 0) { 
	linecomment = i;
	//cout << "linecomment " << linecomment << endl;
	}
	else if (i > linecomment && i < linetypeinfo){
		if (vec_line[1] == "atoms"){
			natoms = stoi (vec_line[0]);	
			ai.natoms = natoms;
			cout << "natoms: " << natoms << endl;
			ai.atomnumber.resize(natoms);
        		ai.atomtype.resize(natoms);
       			ai.xcoords.resize(natoms);
        		ai.ycoords.resize(natoms);
        		ai.zcoords.resize(natoms);
			ai.xi.resize(natoms);
			ai.yi.resize(natoms);
			ai.zi.resize(natoms);
 			//make molecule have natoms length
			lmpi.molecule.resize(natoms,0);
		}		
		else if (vec_line[1] == "atom"){
			atomtypes = stoi(vec_line[0]);
			lmpi.masses.assign(atomtypes,0);
		}
		else if (vec_line[1] == "bonds"){
			nbonds = stoi(vec_line[0]);
			cout << "nbonds: " << nbonds << endl;
			bi.nbonds = nbonds;
			bi.bondnumber.resize(nbonds);
			bi.bondtype.resize(nbonds);
			bi.bondfirst.resize(nbonds);
			bi.bondsecond.resize(nbonds);
		}
		else if (vec_line[1] == "bond"){
			bondtypes = stoi(vec_line[0]);
			linetypeinfo = i+1;
			//cout << "linetypeinfo " << linetypeinfo << endl;
			//lmpi.bondconst =  new double[bondtypes];
			//lmpi.bondcoeff =  new double[bondtypes];
		}
	}
	else if (i > linetypeinfo && i < linebox){
		//cout << "next segment " << endl;
		//cout << "vec_line[2] " << vec_line[2] << endl;
		if (vec_line[2] == "xlo"){
			boxcoords[0][0] = atof (vec_line[0].c_str());
			boxcoords[0][1] = atof (vec_line[1].c_str());
			lmpi.xlo = boxcoords[0][0];
			lmpi.xhi = boxcoords[0][1];
			lmpi.xboxlen = boxcoords[0][1] - boxcoords[0][0];
			//cout << "xlo " << boxcoords[0][0] <<  " xhi " << boxcoords[0][1] << endl;
			cout <<"xboxlen " << lmpi.xboxlen << endl;
		}
		else if (vec_line[2] == "ylo"){
			boxcoords[1][0] = atof (vec_line[0].c_str());
			boxcoords[1][1] = atof (vec_line[1].c_str());
			lmpi.ylo = boxcoords[1][0];
			lmpi.yhi = boxcoords[1][1];
			lmpi.yboxlen = boxcoords[1][1] - boxcoords[1][0];
			cout << "yboxlen " << lmpi.yboxlen << endl;
		}
		else if (vec_line[2] == "zlo"){
			boxcoords[2][0] = atof (vec_line[0].c_str());
			boxcoords[2][1] = atof (vec_line[1].c_str());
			lmpi.zlo = boxcoords[2][0];
			lmpi.zhi = boxcoords[2][1];
			linebox = i+1;
			lmpi.zboxlen = boxcoords[2][1] - boxcoords[2][0];
			cout << "zboxlen " << lmpi.zboxlen << endl;
		}
		//else if (vec_line[3] == "xy"){
		//	linebox = i;
		//	cout << "linebox " << linebox << endl;
		//}
	}
	else if (i > linebox && i < linemasses){
		if (vec_line.size() == 1){
		}
		else {
			nmass = stoi(vec_line[0]);
			//cout << "nmass " << nmass;
			lmpi.masses[(nmass-1)] = atof(vec_line[1].c_str());
			//cout << " mass " << lmpi.masses[(nmass-1)] << endl;
			if (nmass == atomtypes){
			linemasses = i;
			}
		}
	}
	else if (i > linemasses && i < linebondcoeff){
		if (vec_line[0] == "Bond"){
		}	
		else{
			nbtype = stoi(vec_line[0]);
			//lmpi.bondconst[(nbtype-1)] = atof(vec_line[1].c_str());
			//lmpi.bondcoeff[(nbtype-1)]= atof(vec_line[2].c_str());
			//cout << "Bondtype " << nbtype << " " << lmpi.bondconst[(nbtype-1)] << " " << lmpi.bondcoeff[(nbtype-1)] << endl;
			if (nbtype == bondtypes){
			linebondcoeff = i;
			}
		}
	}
	else if (i > linebondcoeff and i < lineatom){
		//cout << "lineatom: " << lineatom << endl;
		if (vec_line[0] == "Atoms"){
		lineatom = i + natoms +4;
		}
		else if (vec_line[0] == "Velocities"){
		lineatom = i + natoms + 2;
		//cout << "skipping velocities " << lineatom << endl;
		}
		else{
			if (vec_line.size() >5 ){
				nindex = stoi(vec_line[0]) - 1;
				ai.atomnumber[nindex] = nindex+1;
				//cout << "Atomnumber: " << ai.atomnumber[nindex] << endl;
        			ai.atomtype[nindex] = stoi(vec_line[2]);
				//cout << "Atomtype: " << ai.atomtype[nindex] << endl;	
				lmpi.molecule[nindex] = stoi(vec_line[1]);
        			ai.xcoords[nindex] = atof (vec_line[3].c_str());
        			ai.ycoords[nindex] = atof (vec_line[4].c_str());
        			ai.zcoords[nindex] = atof (vec_line[5].c_str());
				if (vec_line.size() == 9){
				ai.xi[nindex] = stoi(vec_line[6]);
				ai.yi[nindex] = stoi(vec_line[7]);
				ai.zi[nindex] = stoi(vec_line[8]);	}
				else if (vec_line.size() < 9) {
				ai.xi[nindex] = 0;
				ai.yi[nindex] = 0;
				ai.zi[nindex] = 0;
				
				}
				 //cout << "Atom " <<  ai.atomnumber[nindex] << " " << ai.atomtype[nindex] << " " << ai.xcoords[nindex] << " " << ai.ycoords[nindex] << " " << ai.zcoords[nindex] << endl;
				
			}
			else if (vec_line[0] == "Bonds") {
			//cout << "Starting Bonds Section" << endl;
			linebondsend = i+2 + nbonds;
			//cout << "linebondsend: " << linebondsend << " lineatom "<< lineatom << endl;
			}
		    }
		}
	
	else if (i > lineatom and i < linebondsend){
		//cout << "Bonds section" << endl;
		if (vec_line[0] == "Bonds") {
		linebondsend = i + 2 + nbonds;
		//cout << "Bonds section" << " linebondsend " << linebondsend << endl;
		}
		else if (vec_line.size() == 4){
			bindex = stoi(vec_line[0]) -1;
			bi.bondnumber[bindex] = bindex +1;
			bi.bondtype[bindex] = stoi(vec_line[1]);
			bi.bondfirst[bindex]= stoi(vec_line[2]);
			bi.bondsecond[bindex] = stoi(vec_line[3]);
			//cout << "bindex" << bindex << "Bond " << bi.bondnumber[bindex] << " bondtype " << bi.bondtype[bindex] << " bondfirst " << bi.bondfirst[bindex] << " bondsecond " << bi.bondsecond[bindex] << " lineatom " << lineatom << "linebondsend" << linebondsend << " line " << line << endl;
		} 
	
	}
	}
	}
	//cout << "out of lammpsdatafile loop" << endl;
	lmpi.ai = ai;
	lmpi.bi = bi;
	lmpi.atomtypes = atomtypes;
	lmpi.bondtypes = bondtypes;
    
			
	//cout << "Returning data" << endl;
	
	return lmpi;
}

struct OneCluster{
    vector <int> cluster;
    vector <int> edge1List;
    vector <int> edge2List;
    double mass = 0;
};

struct LoopInfo {
    vector <int> loop;
    vector <int> edge1list;
    vector <int> edge2list; 
};


class Burning {
    public:
    vector <vector <int>> cluster;
    vector <double> xcoords;
    vector <double> ycoords;
    vector <double> zcoords;
    vector <int> burnVal;
    vector <int> burn2Val; 
    vector <int> burn3Val;
    vector <int> elasticBackbone;
    vector <int> Eedgelist1;
    vector <int> Eedgelist2;
    vector <int> loopsClosed;
    vector <int> backbone;
    vector <int> Bedgelist1;
    vector <int> Bedgelist2;
    int clusterSites;
    int loops;
    int edges;
    int shortestPath = 0;
    int sitesE;
    int loopsE;
    int cuttingsE;
    int maxE;
    int sitesB;
    int burningTime;
    int burnedSites; 
    vector <double> realDists;
    vector <double> xdists;
    vector <double> ydists;
    vector <double> zdists;
    double shortestRealPath;
    double shortestX;
    double shortestY;
    double shortestZ;
    int firstBurn(int P1, int P2, double maxDist){
        int t = 1;
        realDists.resize(burnVal.size());
        xdists.resize(burnVal.size());
        ydists.resize(burnVal.size());
        zdists.resize(burnVal.size());
        double xdist, ydist, zdist;
        burnedSites = 0;
        bool neighP2 = false;
    
        burnVal[P1] = t;
        realDists[P1] = 0.0;
        xdists[P1] = 0;
        ydists[P1] = 0;
        zdists[P1] = 0;
        
    
        burnedSites++;
        loops = 0;
        //cout << "burnVal[P1] " << burnVal[P1] << endl;
        //cout << " P1 " << P1 << " P2 " << P2 << endl;
        t++;
        vector <int> neighbors;
        
        //cout << P1 << " ";
        //second burn
        for (int n : cluster[P1]){
            
                
                if (burnVal[n] == 0 && n != P2){
                    //cout << n << " n ";
                    burnVal[n] = t;
                    realDists[n] = realDists[P1] + pow(pow((xcoords[n] - xcoords[P1]),2) + pow((ycoords[n] - ycoords[P1]),2) + pow((zcoords[n] - zcoords[P1]),2) ,0.5);   
                    
                    neighbors.push_back(n);
                    burnedSites++;
                    
                }
                else if (n == P2){
                    neighP2 = true;
                }
                else if (burnVal[n] ==t) {
                        loops++;
                        loopsClosed.push_back(n);
                }
            
        }
        
        //subsequent burns
        vector <int> newNeigh;
        
        t++;
        while (burnedSites < clusterSites && neighbors.size() > 0){
            //cout << "t " << t << endl;
            for (int i = 0; i< neighbors.size(); i++){
                int n = neighbors[i];
                
                for (int j = 0;  j < cluster[n].size(); j++){
                    int v = cluster[n][j];
                   
                    
                    if (burnVal[v] == 0){
                        //cout <<  v << " v ";
                        burnVal[v] = t;
                        burnedSites++;
                        newNeigh.push_back(v);
                        realDists[v] = realDists[n] + pow(pow((xcoords[n] - xcoords[v]),2) + pow((ycoords[n] - ycoords[v]),2) + pow((zcoords[n] - zcoords[v]),2) ,0.5);   
                        xdist = xcoords[v] - xcoords[n];
                        ydist = ycoords[v] - ycoords[n];
                        zdist = zcoords[v] - zcoords[n];
                        
                        while (xdist > maxDist){
                        xdist = xdist - maxDist*2;
                        }
                        while (xdist < -maxDist){
                            xdist = xdist + maxDist*2;
                        }
                        while (ydist > maxDist){
                            ydist = ydist - maxDist*2;
                        }
                        while (ydist < -maxDist){
                            ydist = ydist + maxDist*2;
                        }
                        while (zdist > maxDist){
                            zdist = zdist - maxDist*2;
                        }
                        while (zdist < -maxDist){
                            zdist = zdist + maxDist*2;
                        }
                        xdists[v] = xdists[n] + (xdist);
                        ydists[v] = ydists[n] + (ydist);
                        zdists[v] = zdists[n] + (zdist);
                        //cout << "newNeigh.size() " << newNeigh.size() <<  "clusterSites " << clusterSites << " burnedSites " << burnedSites <<  endl;
                        
                        if (v == P2){
                            shortestPath = t;
                            shortestRealPath = realDists[v];
                            shortestX = abs(xdists[v]);
                            shortestY = abs(ydists[v]);
                            shortestZ = abs(zdists[v]);
                            if (neighP2){
                                loops++;
                                loopsClosed.push_back(v);
                            }
                            //cout << "shortestX " << shortestX << " shortestY " << shortestY << " shortestZ " << shortestZ << endl;
                        }
                    }
                    else if (burnVal[v] == t) {
                        loops++;
                        loopsClosed.push_back(v);
                    }

                }   
                
            }
            
            t++;
            neighbors = newNeigh;
            newNeigh = {};
        }
        burningTime = t;
        //cout << "" << endl;
        
        return t;

    }

    int firstBurnComplete(int P1, int P2){
        int t = 1;
        burnedSites = 0;
        
        burnVal[P1] = t;
    
        burnedSites++;
        loops = 0;
        //cout << "burnVal[P1] " << burnVal[P1] << endl;
        t++;
        vector <int> neighbors;
        
        //second burn
        for (int n : cluster[P1]){
            
            
                if (burnVal[n] == 0 && n != P2){
                    burnVal[n] = t;
                    
                    neighbors.push_back(n);
                    burnedSites++;
                    if (n == P2){
                            shortestPath = t;
                        }
                }
                else if (burnVal[n] ==t) {
                        loops++;
                        loopsClosed.push_back(n);
                }
            
        }
        
        //subsequent burns
        vector <int> newNeigh;
        int nbondsP2 = cluster[P2].size();
        if (nbondsP2 == 1){
            shortestPath = 1;
            return t;

        }
        
        t++;
        while (burnedSites < clusterSites && neighbors.size() > 0){
            //cout << "t " << t << endl;
            for (int i = 0; i< neighbors.size(); i++){
                int n = neighbors[i];
                
                for (int j = 0;  j < cluster[n].size(); j++){
                    int v = cluster[n][j];
                    
                    if (burnVal[v] == 0){
                        burnVal[v] = t;
                        burnedSites++;
                        newNeigh.push_back(v);
                        //cout << "newNeigh.size() " << newNeigh.size() <<  "clusterSites " << clusterSites << " burnedSites " << burnedSites <<  endl;
                        
                        if (v == P2){
                            shortestPath = t;
                        }
                    }
                    else if (burnVal[v] == t) {
                        loops++;
                        loopsClosed.push_back(v);
                    }

                }   
                
            }
            
            t++;
            neighbors = newNeigh;
            newNeigh = {};
        }
        burningTime = t;
        
        return t;

    }

    int secondBurn(int P1, int P2){
        vector <int> vals = {burnVal[P2]};
        elasticBackbone.push_back(P2);
        int valT = burnVal[P2];
        
        int t = 1;
        sitesE = 0;
        int burn2sites = 0;
        burn2Val[P2] = t;
        t++;
    
        sitesE++;
        loopsE  = 0;
        vector <int> neighbors;
        int newVal;
        int burningT = 10;
        maxE = 1;
        cuttingsE = 0;
        vector <int> newVals;
        Eedgelist1 = Bedgelist1;
        Eedgelist2 = Bedgelist2;
        Eedgelist1.push_back(P1);
        Eedgelist2.push_back(P2);
        //second burn
        for (int n : cluster[P2]){
            
                if (burnVal[n] < valT && burn2Val[n] ==0){
                    burn2Val[n] = t;
                    elasticBackbone.push_back(n);
                    Eedgelist1.push_back(P2);
                    Eedgelist2.push_back(n);
                    newVals.push_back(burnVal[n]);
                    //cout << " n " << n << " burnVal[n] " << burnVal[n] << endl;
                    neighbors.push_back(n);
                    sitesE++;
                    burn2sites++;
                    burningT++;
                }
                else if (burn2Val[n] == 0){
                    burn2sites++;
                }
                else if (burn2Val[n] == t){
                    loopsE++;
                }
                
            
        }
        if (burningT  > maxE){
            maxE = burningT;
        }
        if (burningT == 1){
            cuttingsE++;
        }
        burningT = 0;

        vals = newVals;
        newVals = {}; 
        //subsequent burns
        vector <int> newNeigh;
        int nbondsP1 = cluster[P1].size();
        if (nbondsP1 == 1){
            return t;
        }
        t++;
        valT--;
        while (neighbors.size() >= 1 and burn2sites < clusterSites){
            //cout << "t " << t << endl;
            for (int i = 0; i< neighbors.size(); i++){
                int n = neighbors[i];
                //cout << " n " << n << "neighbors.size " << neighbors.size() << " i " << i << endl;
                for (int v: cluster[n]){
                    if (burnVal[v] < valT  && burn2Val[v] == 0){
                        burn2Val[v] = t;
                        elasticBackbone.push_back(v);
                        Eedgelist1.push_back(n);
                        Eedgelist2.push_back(v);
                        newVals.push_back(burnVal[v]);
                        //cout << "i " << i << " newVal " << newVal << endl;
                        sitesE++;
                        newNeigh.push_back(v);
                        burn2sites++;
                        burningT++;
                        
                    }
                    else if (burn2Val[v] == 0){
                        burn2sites++;
                    }
                    else if (burn2Val[v] == t){
                        loopsE++;
                    }
                }            
            }
            t++;
            valT--;
            neighbors = newNeigh;
            newNeigh = {};
            vals = newVals;
            newVals = {};
            if (burningT  > maxE){
            maxE  = burningT;
            }
            if (burningT == 1){
            cuttingsE++;
            }
            burningT = 0;


        }

        return t;

    }
    
    int thirdBurn(int P1, int P2){
        vector <int> neighbors = loopsClosed; 
        vector <int> newNeigh;
        backbone = elasticBackbone;
        Bedgelist1 = Eedgelist1;
        Bedgelist2 = Eedgelist2;
        int burn3sites = 0;
        sitesB = sitesE;
        int t = 1;
        vector <int> vals;  
        vector <int> newVals;
        int reachedEsites;

        //burn from Pi, loopsClosed

        for (int i = 0; i < loopsClosed.size(); i++){
            t = 1;
            
            vector <int> siteList;
            vector <int> siteE1;
            vector <int> siteE2;
            reachedEsites = 0;
            int n = loopsClosed[i];
            int valT = burnVal[n];
        
            if (find(elasticBackbone.begin(), elasticBackbone.end(),n) == elasticBackbone.end()){
                burn3Val[n] = t;
                burn3sites++;
                vals.push_back(burnVal[n]);
                siteList.push_back(n);
                t++;

                for (int v: cluster[n]){
                    if (burn3Val[v] == 0 && burnVal[v] < valT  && find(elasticBackbone.begin(), elasticBackbone.end(), v) == elasticBackbone.end()){
                        burn3Val[v] = t;
                        burn3sites++;
                        neighbors.push_back(v);
                        newVals.push_back(burnVal[v]);
                        siteList.push_back(v);
                        siteE1.push_back(n);
                        siteE2.push_back(v);
                    }
                    else if (find(elasticBackbone.begin(), elasticBackbone.end(),v) != elasticBackbone.end()) {
                        reachedEsites++;
                        siteE1.push_back(n);
                        siteE2.push_back(v);

                    }
                    
                }
                t++;
                valT--;
                vals = newVals;
                newVals = {};
                while (neighbors.size() > 0 && burn3sites < clusterSites){
                    for (int j = 0; j < neighbors.size(); j++){
                        int u = neighbors[j];
                        for (int w: cluster[u]){
                            if (burn3Val[w] == 0 && burnVal[w] < valT && find(elasticBackbone.begin(), elasticBackbone.end(), w) == elasticBackbone.end()){
                                burn3Val[w] = t;
                                burn3sites++;
                                newNeigh.push_back(w);
                                newVals.push_back(burn3Val[w]);
                                siteList.push_back(w);
                                siteE1.push_back(u);
                                siteE2.push_back(w);
                                //cout << "added site " << endl;
                            }
                            else if (find(elasticBackbone.begin(), elasticBackbone.end(), w) != elasticBackbone.end()){
                                reachedEsites++;
                                siteE1.push_back(u);
                                siteE2.push_back(w);
                            }
                            
                        }

                    }
                    neighbors = newNeigh;
                    newNeigh = {};
                    vals = newVals;
                    newVals= {};
                    t++;
                    valT--;
                }
            }
            //cout << " reachedEsites " << reachedEsites << endl;
            if (reachedEsites > 1){
                //add all sites to growing backbone
                backbone.insert(backbone.end(), siteList.begin(), siteList.end());
                Bedgelist1.insert(Bedgelist1.end(), siteE1.begin(), siteE1.end());
                Bedgelist2.insert(Bedgelist2.end(), siteE2.begin(), siteE2.end());
                sitesB += siteList.size();

            }

        }
            
       
        
        
        
        return t;

    }
    
};

struct DanglingInfo{
    int danglingEnds;
    int xlinkEnds;
    int primaryEnds;
    int secondaryEnds;
}; 

struct ActiveInfo{
    int nActiveLinks;
    int nActiveStrands;
    vector <int> lensActive;
    double avgActiveLen;
    double avgActiveMass;
    int xLinks;
};

class GraphC {
public:
    
    vector <vector<int>> graph;
    vector <int> atominds;
    int edges = 0;
    vector <double> xcoords;
    vector <double> ycoords;
    vector <double> zcoords;
    vector <double> masses;
    vector <bool> marked;
    vector <int> atomtypes;

    vector <bool> xLinker;
    vector <vector <int>> actualAinds;
    bool foundLoop = false;

    GraphC(int n){
        graph.resize(n);
        xcoords.resize(n);
        ycoords.resize(n);
        zcoords.resize(n);
        masses.resize(n);
        marked.resize(n);
        xLinker.resize(n);
        atomtypes.resize(n);


    }

    class ClustersInfo{
    public: 
        vector <vector <int>> allClusters;
        int clusters; 
        vector <int> clusterSizes;
        vector <vector <int>> edge1Lists;
        vector <vector <int>> edge2Lists; 
        vector <double> Cmasses; 
        vector <GraphC> clusterGraphs;
    };

    ClusterSize checkSize(double xlo, double xhi, double ylo, double yhi, double zlo, double zhi){
        
        double minX;
        double minY;
        double minZ;
        double maxX;
        double maxY;
        double maxZ;
        int P1x;
        int P2x; 
        int P1y;
        int P2y;
        int P1z;
        int P2z; 
        double boxlenx = xhi - xlo;
        double boxleny = yhi - ylo;
        double boxlenz = zhi - zlo;
        bool x1s = false;
        bool x2s = false;
        bool x3s = false;
        bool x4s = false;
        bool x5s = false;
        bool x6s = false;
        bool x7s = false;
        bool x8s = false;
        bool x9s = false;
        bool x10s = false;
        bool y1s = false;
        bool y2s = false;
        bool y3s = false;
        bool y4s = false;
        bool y5s = false;
        bool y6s = false;
        bool y7s = false;
        bool y8s = false;
        bool y9s = false;
        bool y10s = false;
        bool z1s = false;
        bool z2s = false;
        bool z3s = false;
        bool z4s = false;
        bool z5s = false;
        bool z6s = false;
        bool z7s = false;
        bool z8s = false;
        bool z9s = false;
        bool z10s = false;
        
        for (int i = 0; i < atominds.size(); i++){
            int atomind = atominds[i];
            double x = xcoords[atomind];
            double xic = x - xlo;
            double y = ycoords[atomind];
            double yic = y - ylo;
            double z = zcoords[atomind];
            double zic = z - zlo;
            if (x > maxX){
                maxX = x;
                P2x = atomind;
            }
            if (x < minX){
                minX = x;
                P1x = atomind;
            }
            if (y > maxY){
                maxY = y;
                P2y = atomind;
            }
            if (y < minY){
                minY = y;
                P1y = atomind;
            }
            if (z > maxZ){
                maxZ = z;
                P2z = atomind;
            }
            if (z < minZ){
                minZ =  z;
                P1z = atomind;
            }
            if (xic >= boxlenx*0  && xic < boxlenx*0.1){
                x1s = true;
            }
            else if (xic >= boxlenx*0.1 && xic < boxlenx*0.2){
                x2s = true;
            }
            else if (xic >= boxlenx*0.2 && xic < boxlenx*0.3){
                x3s = true;
            }
            else if (xic >= boxlenx*0.3 && xic < boxlenx*0.4){
                x4s = true;
            }
            else if (xic >= boxlenx*0.4 && xic < boxlenx*0.5){
                x5s = true;
            }
            else if (xic >= boxlenx*0.5 && xic < boxlenx*0.6){
                x6s = true;
            }
            else if (xic >= boxlenx*0.6 && xic < boxlenx*0.7){
                x7s = true;
            }
            else if (xic >= boxlenx*0.7 && xic < boxlenx*0.8){
                x8s = true;
            }
            else if (xic >= boxlenx*0.8 && xic < boxlenx*0.9){
                x9s = true;
            }
            else if (xic >= boxlenx*0.9 && xic < boxlenx){
                x10s = true;
            }

            if (yic >= boxleny*0  && yic < boxleny*0.1){
                y1s = true;
            }
            else if (yic >= boxleny*0.1 && yic < boxleny*0.2){
                y2s = true;
            }
            else if (yic >= boxleny*0.2 && yic < boxleny*0.3){
                y3s = true;
            }
            else if (yic >= boxleny*0.3 && yic < boxleny*0.4){
                y4s = true;
            }
            else if (yic >= boxleny*0.4 && yic < boxleny*0.5){
                y5s = true;
            }
            else if (yic >= boxleny*0.5 && yic < boxleny*0.6){
                y6s = true;
            }
            else if (yic >= boxleny*0.6 && yic < boxleny*0.7){
                y7s = true;
            }
            else if (yic >= boxleny*0.7 && yic < boxleny*0.8){
                y8s = true;
            }
            else if (yic >= boxleny*0.8 && yic < boxleny*0.9){
                y9s = true;
            }
            else if (yic >= boxleny*0.9 && yic < boxleny){
                y10s = true;
            }

            if (zic >= boxlenz*0  && zic < boxlenz*0.1){
                z1s = true;
            }
            else if (zic >= boxlenz*0.1 && zic < boxlenz*0.2){
                z2s = true;
            }
            else if (zic >= boxlenz*0.2 && zic < boxlenz*0.3){
                z3s = true;
            }
            else if (zic >= boxlenz*0.3 && zic < boxlenz*0.4){
                z4s = true;
            }
            else if (zic >= boxlenz*0.4 && zic < boxlenz*0.5){
                z5s = true;
            }
            else if (zic >= boxlenz*0.5 && zic < boxlenz*0.6){
                z6s = true;
            }
            else if (zic >= boxlenz*0.6 && zic < boxlenz*0.7){
                z7s = true;
            }
            else if (zic >= boxlenz*0.7 && zic < boxlenz*0.8){
                z8s = true;
            }
            else if (zic >= boxlenz*0.8 && zic < boxlenz*0.9){
                z9s = true;
            }
            else if (zic >= boxlenz*0.9 && zic < boxlenz){
                z10s = true;
            }
        
        }

        ClusterSize cSize;
        if (x1s && x2s && x3s && x4s && x5s && x6s && x7s && x8s && x9s && x10s){
            cSize.bigx = true;
            cSize.P1 = P1x;
            cSize.P2 = P2x;
        }
        if (y1s && y2s && y3s && y4s && y5s && y6s && y7s && y8s && y9s && y10s){
            cSize.bigy = true;
            cSize.P1 = P1y;
            cSize.P2 = P2y;
        }
        if (z1s && z2s && z3s && z4s && z5s && z6s && z7s && z8s && z9s && z10s){
            cSize.bigz = true;
            cSize.P1 = P1z;
            cSize.P2 = P2z;
        }
        if (cSize.bigx || cSize.bigy || cSize.bigz){
            cSize.bigenough = true;
        }
        return cSize;
    }
    bool isCyclicUtil(int v, vector <bool> &visited, int parent) { 
    // Mark the current node as visited 
    visited[v] = true; 
    //cout << " v " << v << " parent " << parent << endl;
  
    // Recur for all the vertices adjacent to this vertex 
    vector<int>::iterator i; 
    for (i = graph[v].begin(); i != graph[v].end(); ++i) 
    { 
        // If an adjacent is not visited, then recur for that adjacent 
        if (!visited[*i]) 
        { 
            
           if (isCyclicUtil(*i, visited, v)) {
              return true; 
           }
        } 
  
        // If an adjacent is visited and not parent of current vertex, 
        // then there is a cycle. 
        else if (*i != parent) {
           return true; 
        }
    } 
    return false; 
    } 
    bool isCyclic() { 
    // Mark all the vertices as not visited and not part of recursion 
    // stack 
    vector <bool> visited (graph.size()); 
    for (int i = 0; i < graph.size(); i++) {
        visited[i] = false; 
    }
  
    // Call the recursive helper function to detect cycle in different 
    // DFS trees 
    for (int u = 0; u < graph.size(); u++) {
        
        if (!visited[u]) {// Don't recur for u if it is already visited 
          if (isCyclicUtil(u, visited, -1)){ 
             cout << "cycle found" << endl;
             return true; 

          }
        }
    }
  
    return false; 
    } 
    LoopInfo findCycle(int v, int u, LoopInfo &path) {

        marked[v] = true;
        path.loop.push_back(v);

        for (int w : graph[v] ){
            if(!marked[w]) {
                marked[w] = true;
                path.edge1list.push_back(v);
                path.edge2list.push_back(w); 
                findCycle(w,v, path);
            } else if (v != u && path.loop.size() >= 3) {
                foundLoop = true;
                return path;
            }
        }
        return path;
    }

    bool visited (int node, vector <int> path){
        if (find(path.begin(), path.end(), node) != path.end()){
            return true;
        }
        else {
            return false;
        }
    }

    vector <int> invert(vector <int> path){
        reverse(path.begin(), path.end());
        return path;
    }
    vector <int> rotate_to_smallest(vector <int> path){
        int n = find(path.begin(), path.end(), *min_element(path.begin(), path.end())) - path.begin();
        vector <int> newPath;
        newPath.insert(newPath.end(), path.begin(), path.begin()+n);
       
        newPath.insert(newPath.end(), path.begin()+n, path.end());
        return newPath; 

    }

    bool isNew(vector <int> path, vector <LoopInfo> allLoops){
        for (int i = 0; i < allLoops.size(); i++){
            if (allLoops[i].loop == path){
                return false;
            }
        }
        return true;

    }
    vector <LoopInfo> findNewCycles(vector <int> &path, vector <LoopInfo> & allLoops){
        int start_node  = path[0];
        int next_node;
        vector <int> sub;
    

        for (int u = 0; u < graph.size(); u++){
            for (int v: graph[u]){
                int node1 = u;
                int node2 = v;
                if (start_node == node1 || start_node == node2){
                    if (node1 == start_node){
                        next_node = node2;
                    }
                    else {
                        next_node = node1;
                    }
                    if (!visited(next_node, path)){
                        sub = {next_node};
                        sub.insert(sub.end(), path.begin(), path.end());
                        allLoops = findNewCycles(sub, allLoops);

                    }
                    else if (path.size() > 2 and next_node == path[path.size()-1]){
                        //cycle found
                        vector <int> p = rotate_to_smallest(path);
                        vector <int> inv = invert(p);
                        if (isNew(p, allLoops) && isNew(inv, allLoops)){
                            LoopInfo loop;
                            loop.loop = p;
                            allLoops.push_back(loop);
                        }
                    }
                }

            }
        }
        return allLoops;


    }
    vector <LoopInfo> findAllCycles(){
        LoopInfo loop;

    }

    DanglingInfo removeDanglingEnds(){
        //removes dangling ends from cluster graph and changes graph
        DanglingInfo dInfo;
        dInfo.danglingEnds = 0;
        dInfo.primaryEnds = 0;
        dInfo.secondaryEnds = 0;
        dInfo.xlinkEnds = 0;

        for (int u = 0; u < graph.size(); u++){
            if (xLinker[u]){
                int nxlinkB = graph[u].size();
                for (int v: graph[u]){
                    if (nxlinkB == 1){
                        deleteEdge(u,v);
                        dInfo.danglingEnds++;
                        int atype = atomtypes[u];
                        if (atype == 3 || atype == 1 || atype == 10 || atype == 7 || atype == 14 || atype ==11){
                            dInfo.primaryEnds++;
                        }
                        else if(atype == 6 || atype == 4) {
                            dInfo.secondaryEnds++;

                        }
                        else if (atype == 25){
                            dInfo.xlinkEnds++;
                        }
                    }
                }
            }
        }
        for (int u = 0; u < graph.size(); u++){
            int nBonds = graph[u].size();
            if (nBonds == 1){
                deleteEdge(u, graph[u][0]);
            }
        }
        return dInfo;

    }
    
   vector <LoopInfo> findP1P2(double Lx, double Ly, double Lz){
        //check if loop exists, if not return 
        bool loopExists = isCyclic();
        cout << " loopExists " << loopExists << endl;
        vector <LoopInfo> loops;
        if (loopExists){
            cout << "loopExists" << endl;
            DanglingInfo danglingEnds = removeDanglingEnds();
            /*
            for (int u = 0; u < graph.size(); u++){
                for (int v : graph[u]){
                    vector <int> path = {u};
                    loops = findNewCycles(path,loops);
                    path = {v};
                    loops= findNewCycles(path, loops);
                }
            } 
            */  
        }
        else {
            return loops;
        }
        return loops;
    }
    
    void addEdge(int u, int v) {
        if (find(graph[u].begin(), graph[u].end(), v) == graph[u].end()){
            graph[u].push_back(v);
            if (find(atominds.begin(), atominds.end(), u) == atominds.end()){
                atominds.push_back(u);
            }
        
            if (find(graph[v].begin(), graph[v].end(), u) == graph[v].end()){
                 graph[v].push_back(u);
                 if (find(atominds.begin(), atominds.end(),v) == atominds.end()){
                     atominds.push_back(v);
                 }
                 edges++;
            }
        }   
    }
    void deleteEdge(int u, int v) {
        auto it1 = find(graph[u].begin(), graph[u].end(), v);
        graph[u].erase(it1);
        auto it2 = find(graph[v].begin(), graph[v].end(), u);
        graph[v].erase(it2);
    }
    P1P2list boundaryBonds(double len){
        vector <int> P1;
        vector <int> P2;
        P1P2list p1p2;
        for (int u = 0; u < graph.size(); u++){
            for (int v: graph[u]){
                double xdist = abs(xcoords[u] - xcoords[v]);
                double ydist = abs(ycoords[u] - ycoords[v]);
                double zdist = abs(zcoords[u] - zcoords[v]);
                if (xdist > len || ydist > len || zdist > len){
                    P1.push_back(u);
                    P2.push_back(v);
                } 
            }
        }
        p1p2.P1 = P1;
        p1p2.P2 = P2;
        return p1p2;
            
    }

    ClustersInfo connectedComponents() { 
    // Mark all the vertices as not visited 
    int V = graph.size();
    //cout << " V " << V << endl;
    vector <bool> visited (V,true); 
    //cout << " atominds.size() " << atominds.size() << endl;
    
    for(int v = 0; v < atominds.size(); v++) {
        
        visited[atominds[v]] = false;
    }
    
    int clusters = 0; 
    vector <vector <int>> allClusters;
    vector <int> clusterSizes; 
    vector <vector <int>> edge1Lists;
    vector <vector <int>> edge2Lists;
    vector <double> Cmasses;
    vector <GraphC> cGraphs;
    for (int v=0; v<V; v++) 
    { 
        if (visited[v] == false) 
        { 
            // print all reachable vertices 
            // from v 
            //cout << " v " << v << endl;
            OneCluster ocluster;
                   
            ocluster = DFSUtil(v, visited,ocluster, -1); 
        
            //cout << "ocluster.cluster.size() " << ocluster.cluster.size() << endl;
            allClusters.push_back(ocluster.cluster);
            clusters++;
            clusterSizes.push_back(ocluster.cluster.size());
            GraphC gOneCluster = GraphC(graph.size());
            edge1Lists.push_back(ocluster.edge1List);
            edge2Lists.push_back(ocluster.edge2List);
            //cout << "cluster mass " << ocluster.mass << endl;
            Cmasses.push_back(ocluster.mass);
            
            for (int i = 0; i < ocluster.cluster.size(); i++){
                
                int atomind = ocluster.cluster[i];
                //cout << " vertex i " << i << " atomind " << atomind << " graph.size() " << graph.size() << endl;
                gOneCluster.xcoords[atomind] = xcoords[atomind];
                //cout << "xcoord" << endl;
                gOneCluster.ycoords[atomind] = ycoords[atomind];
                //cout << "ycoords" << endl;
                gOneCluster.zcoords[atomind] = zcoords[atomind];
                gOneCluster.masses[atomind] = masses[atomind];
                gOneCluster.xLinker[atomind] = xLinker[atomind];
                gOneCluster.atominds.push_back(atomind);
                gOneCluster.atomtypes[atomind] = atomtypes[atomind];
                if (masses[atomind] < 1){
                    cout << " atomind " << atomind << " masses[atomind] " << masses[atomind] << endl;
                }
                //cout << "xLinker " << endl;
            }
            for (int i = 0; i < ocluster.edge1List.size(); i++){
                //cout << "edge i " << i << " " << ocluster.edge1List[i] << " " << ocluster.edge2List[i] << endl;
                gOneCluster.addEdge(ocluster.edge1List[i], ocluster.edge2List[i]);
            }
            cGraphs.push_back(gOneCluster);
            
            //cout << "\n"; 
        } 
    }
    ClustersInfo ci;
    ci.allClusters = allClusters;
    ci.clusters = clusters;
    ci.clusterSizes = clusterSizes;
    ci.edge1Lists = edge1Lists;
    ci.edge2Lists = edge2Lists; 
    ci.Cmasses = Cmasses;
    ci.clusterGraphs = cGraphs;

    return ci;
    } 

    ClustersInfo connectedActive(vector <int> activeXlinks) { 
    // Mark all the vertices as not visited 
    int V = graph.size();
    //cout << " V " << V << endl;
    vector <bool> visited (V,true); 
    //cout << " atominds.size() " << atominds.size() << endl;
    
    for(int v = 0; v < atominds.size(); v++) {
        
        visited[atominds[v]] = false;
    }
    
    int clusters = 0; 
    vector <vector <int>> allClusters;
    vector <int> clusterSizes; 
    vector <vector <int>> edge1Lists;
    vector <vector <int>> edge2Lists;
    vector <double> Cmasses;
    vector <GraphC> cGraphs;
    for (int i=0; i < activeXlinks.size(); i++) 
    { 
        int v = activeXlinks[i];
        if (visited[v] == false) 
        { 
            // print all reachable vertices 
            // from v 
            //cout << " v " << v << endl;
            OneCluster ocluster;
                   
            ocluster = DFSUtil(v, visited,ocluster, -1); 
        
            //cout << "ocluster.cluster.size() " << ocluster.cluster.size() << endl;
            allClusters.push_back(ocluster.cluster);
            clusters++;
            clusterSizes.push_back(ocluster.cluster.size());
            GraphC gOneCluster = GraphC(graph.size());
            edge1Lists.push_back(ocluster.edge1List);
            edge2Lists.push_back(ocluster.edge2List);
            //cout << "cluster mass " << ocluster.mass << endl;
            Cmasses.push_back(ocluster.mass);
            
            for (int i = 0; i < ocluster.cluster.size(); i++){
                
                int atomind = ocluster.cluster[i];
                //cout << " vertex i " << i << " atomind " << atomind << " graph.size() " << graph.size() << endl;
                gOneCluster.xcoords[atomind] = xcoords[atomind];
                //cout << "xcoord" << endl;
                gOneCluster.ycoords[atomind] = ycoords[atomind];
                //cout << "ycoords" << endl;
                gOneCluster.zcoords[atomind] = zcoords[atomind];
                gOneCluster.masses[atomind] = masses[atomind];
                gOneCluster.xLinker[atomind] = xLinker[atomind];
                gOneCluster.atominds.push_back(atomind);
                gOneCluster.atomtypes[atomind] = atomtypes[atomind];
                if (masses[atomind] < 1){
                    cout << " atomind " << atomind << " masses[atomind] " << masses[atomind] << endl;
                }
                //cout << "xLinker " << endl;
            }
            for (int i = 0; i < ocluster.edge1List.size(); i++){
                //cout << "edge i " << i << " " << ocluster.edge1List[i] << " " << ocluster.edge2List[i] << endl;
                gOneCluster.addEdge(ocluster.edge1List[i], ocluster.edge2List[i]);
            }
            cGraphs.push_back(gOneCluster);
            
            //cout << "\n"; 
        } 
    }
    ClustersInfo ci;
    ci.allClusters = allClusters;
    ci.clusters = clusters;
    ci.clusterSizes = clusterSizes;
    ci.edge1Lists = edge1Lists;
    ci.edge2Lists = edge2Lists; 
    ci.Cmasses = Cmasses;
    ci.clusterGraphs = cGraphs;

    return ci;
    } 

    OneCluster DFSUtil(int v, vector <bool> &visited, OneCluster &ocluster, int parent) { 
    // Mark the current node as visited and print it 
    visited[v] = true; 
    //cout << v << " " ;
    ocluster.cluster.push_back(v);
    //cout << ocluster.cluster.size() << "* ";
    ocluster.mass = ocluster.mass +  masses[v];
    
    
    //cout << "mass " << masses[v] << endl;
  
    // Recur for all the vertices 
    // adjacent to this vertex 
    vector<int>::iterator i; 
    for(i = graph[v].begin(); i != graph[v].end(); ++i) {
        //cout << "graph[v].size() " << graph[v].size() << endl;
        if(!visited[*i]) {
            //cout << " v " << v << " *i " << *i << endl;
            //cout << "connected ocluster.cluster.size() " << ocluster.cluster.size() <<  endl;
           
            
            ocluster.edge1List.push_back(v);
            ocluster.edge2List.push_back(*i);
            DFSUtil(*i, visited, ocluster, v); 
        }
        else if (*i != parent){
            ocluster.edge1List.push_back(v);
            ocluster.edge2List.push_back(*i);
            

        }
    }
    //cout << "ocluster.cluster.size() " << ocluster.cluster.size() << endl;
    return ocluster;

    }

    ActiveInfo activeAnalysis(){
        ActiveInfo actI;
        actI.xLinks = 0;
        actI.nActiveLinks = 0;
        actI.nActiveStrands = 0;
        actI.lensActive; 
        actI.avgActiveLen = 0;
        actI.avgActiveMass = 0;
       
        vector <int> activeXlinks;
        //first find nActiveLinks
        for (int i=0; i < atominds.size(); i++){
            int aind = atominds[i];
            if (xLinker[aind]){
                cout << "FOUND XLINKER" << endl;
                actI.xLinks++;
                int nbonds = graph[aind].size();
                cout << "actI.xLinks " << actI.xLinks << endl;
                cout << " nbonds " << nbonds << endl;
                if (nbonds < 3){
                    activeXlinks.push_back(aind);
                    actI.nActiveLinks++;
                }
            }
        }
        //next find clusters 
        ClustersInfo activeChains = connectedActive(activeXlinks);
        actI.nActiveStrands = activeChains.clusters;

     
        for (int i = 0; i < activeChains.clusterSizes.size(); i++){
            actI.lensActive.push_back(activeChains.clusterSizes[i]);
            actI.avgActiveLen += activeChains.clusterSizes[i];
            actI.avgActiveMass += activeChains.Cmasses[i];
        }
        actI.avgActiveLen = actI.avgActiveLen/actI.nActiveStrands;
        actI.avgActiveMass = actI.avgActiveMass/actI.nActiveStrands;
        return actI;

    }


 
};


GraphC consolidateCrosslinker(Lmpsabinfo lmpi, vector <int> xlinktypes, double xLinkMass, int solventType){
    cout << "consolidate Crosslinker" << endl;
	vector <int> bondnumber;
	vector <int> bondtype;
	vector <int> bondfirst;
	vector <int> bondsecond;
    vector <int> foundxmol;

    
    int notXlinkAtoms = 0;
    int newAtoms = 0;
    vector <double> newXcoords;
    vector <double> newYcoords;
    vector <double> newZcoords;
    
    

    
    for (int i =0; i < lmpi.ai.natoms; i++){
        int atype = lmpi.ai.atomtype[i];
        //cout << " atype " << atype << endl;
        if (find(xlinktypes.begin(), xlinktypes.end(), atype) == xlinktypes.end() && atype != solventType){
            newXcoords.push_back(lmpi.ai.xcoords[i]);
            newYcoords.push_back(lmpi.ai.ycoords[i]);
            newZcoords.push_back(lmpi.ai.zcoords[i]);
            notXlinkAtoms++;
            newAtoms++;
            
          
            //cout << "not crosslinker " << endl;
        }    
        else if (find(xlinktypes.begin(), xlinktypes.end(), atype) != xlinktypes.end()) {
            int moli = lmpi.molecule[i]; 
            //cout << "crosslinker" << endl;
            
                if (find(foundxmol.begin(),foundxmol.end(),moli) == foundxmol.end()){	
	                foundxmol.push_back(moli);
                    newAtoms++;
                }

        }
    }

    vector <double> masses (newAtoms, 0);
    vector <bool> xLinker (newAtoms,false);
    vector <int> newAtominds;
    vector <int> newatomtypes (newAtoms);

    vector <vector <int>> newActualAtominds (newAtoms);
    
    for (int i =0; i < lmpi.ai.natoms; i++){
        int atype = lmpi.ai.atomtype[i];
        int aind = lmpi.ai.atomnumber[i]-1;
        if (aind != i){
            cout << "aind " << aind << " doesn't equal i" <<i << endl;
        }
        assert (aind == i);
    
        if (find(xlinktypes.begin(), xlinktypes.end(), atype) == xlinktypes.end() && atype != solventType){
            masses[aind] = lmpi.masses[atype-1];
            newatomtypes[aind] = atype;
            if (find(newAtominds.begin(), newAtominds.end(), i) == newAtominds.end()){
                newAtominds.push_back(aind);
            }
            newActualAtominds[i].push_back(i);

        }
        else if (find(xlinktypes.begin(), xlinktypes.end(), atype) != xlinktypes.end()){
            int cmol = lmpi.molecule[i];
	        int cind = find(foundxmol.begin(),foundxmol.end(),cmol) - foundxmol.begin();
            
            if (find(newAtominds.begin(), newAtominds.end(), notXlinkAtoms+cind) == newAtominds.end()){
                newAtominds.push_back(notXlinkAtoms + cind);
                masses[notXlinkAtoms +cind] = xLinkMass;
                xLinker[notXlinkAtoms + cind] = true;
                newatomtypes[notXlinkAtoms + cind] = 25;
            }
            
        }
    }
    
    //deal with xlinker and bonds 
    int nxlinkers = foundxmol.size();
    int molc;
    int sind;
    vector <int> avgxnew;
    vector <int> avgynew;
    vector <int> avgznew;
    for (int i=0; i<nxlinkers; i++){
        molc = foundxmol[i];
        sind = 0;
	    vector <int> currentxlinker;

        while (find(lmpi.molecule.begin()+sind,lmpi.molecule.end(),molc) != lmpi.molecule.end()){
	        //cout << " in while loop " << "molc " << molc << "sind " << sind <<endl;
	        sind = find(lmpi.molecule.begin()+sind+1,lmpi.molecule.end(),molc) - lmpi.molecule.begin();
	        if (find(xlinktypes.begin(),xlinktypes.end(),lmpi.ai.atomtype[sind]) != xlinktypes.end()){ //if same mole and crosslinker type
	        currentxlinker.push_back(lmpi.ai.atomnumber[sind]);
	        //cout << "currentxlinker pushed back " << lmpi.ai.atomnumber[sind] << endl;
            newActualAtominds[notXlinkAtoms + i].push_back(lmpi.ai.atomnumber[sind]-1);
	        }
	    }

        int nbeadsxlink = currentxlinker.size();
	    //cout << "nbeadsxlink " << nbeadsxlink << endl;
        double xabsrelcoords [nbeadsxlink];
	    double yabsrelcoords [nbeadsxlink];
	    double zabsrelcoords [nbeadsxlink];
	    double xdist;
	    double ydist;
	    double zdist;
	    //lmpi.xboxlen,lmpi.yboxlen,lmpi.zboxlen
	
	    double avgx = 0;
	    double avgy = 0;
	    double avgz = 0;
	    double xlinkmass=0;

        vector <int> actualanumsxlinker;
	    xabsrelcoords[0] = lmpi.ai.xcoords[currentxlinker[0]-1];
	    yabsrelcoords[0] = lmpi.ai.ycoords[currentxlinker[0]-1];
	    zabsrelcoords[0] = lmpi.ai.zcoords[currentxlinker[0]-1];
	    actualanumsxlinker.push_back(lmpi.ai.atomnumber[currentxlinker[0]-1]);

        for (int j=1; j<nbeadsxlink; j++){
		//cout << "nbeadsxlink " << nbeadsxlink << " j " << j << endl;
		xdist = lmpi.ai.xcoords[currentxlinker[j]-1] - xabsrelcoords[j-1];
		ydist = lmpi.ai.ycoords[currentxlinker[j]-1] - yabsrelcoords[j-1];
		zdist = lmpi.ai.zcoords[currentxlinker[j]-1] - zabsrelcoords[j-1];
		//xlinkmass = xlinkmass + lmpi.masses[lmpi.ai.atomtype[j]-1];
		actualanumsxlinker.push_back(lmpi.ai.atomnumber[currentxlinker[j]-1]);
		while (xdist > lmpi.xboxlen/2){
			xdist = xdist - lmpi.xboxlen;
			//cout << " while xdist " << xdist << endl;
		}
		while (xdist <= -lmpi.xboxlen/2){
			xdist = xdist + lmpi.xboxlen;
			//cout << " while xdist " << xdist << endl;
			
		}
		while (ydist > lmpi.yboxlen/2){
			ydist = ydist - lmpi.yboxlen;
			//cout << " while ydist " << ydist << endl;
		}
		while (ydist <= -lmpi.yboxlen/2){
			ydist = ydist + lmpi.yboxlen;
			//cout << " while ydist " << ydist << endl;
		}
		while (zdist > lmpi.zboxlen/2){
			zdist = zdist - lmpi.zboxlen;
			//cout << " while zdist " << zdist << endl;
		}
		while (zdist <= lmpi.zboxlen/2){
			zdist = zdist + lmpi.zboxlen;
		}	//cout << " while zdist " << zdist << endl;
		xabsrelcoords[j] = xabsrelcoords[j-1] + xdist;
		yabsrelcoords[j] = yabsrelcoords[j-1] + ydist;
		zabsrelcoords[j] = zabsrelcoords[j-1] + zdist;
	    }

        //calculate absolute average
	for (int j=0; j<nbeadsxlink;j++){
	//cout << "calculate abs. avg. j " << j << endl;
	avgx = avgx + xabsrelcoords[j];
	avgy = avgy + yabsrelcoords[j];
	avgz = avgz + zabsrelcoords[j];
	}
	
	avgx = avgx/double(nbeadsxlink);
	avgy = avgy/double(nbeadsxlink);
	avgz = avgz/double(nbeadsxlink);
	
	//put avg coords in pbc box
	
	while(avgx < lmpi.xlo){
		avgx = avgx + lmpi.xboxlen;
	}
	while(avgx >= lmpi.xhi){
		avgx = avgx - lmpi.xboxlen;
	}
	while(avgy < lmpi.ylo){
		avgy = avgy + lmpi.yboxlen;
	}
	while(avgy >= lmpi.yhi){
		avgy = avgy - lmpi.yboxlen;
	}
	while(avgz < lmpi.zlo){
		avgz = avgz + lmpi.zboxlen;
	}
	while(avgz >= lmpi.zhi){
		avgz = avgz - lmpi.zboxlen;
	}
        newXcoords.push_back(avgx);
        newYcoords.push_back(avgy);
        newZcoords.push_back(avgz);
        
        avgxnew.push_back(avgx);
	    avgynew.push_back(avgy);
	    avgznew.push_back(avgz);
    }

    //now deal with all the bonds and graph
    GraphC gspace = GraphC(newAtoms);
    
    gspace.xcoords = newXcoords;
    gspace.ycoords = newYcoords;
    gspace.zcoords = newZcoords;
    gspace.atominds = newAtominds;
    
    vector <int> tookcarexlinkmol;
    vector <int> xlinkfunctionality(foundxmol.size(),0);
    //cout << "foundxmol size" << foundxmol.size() << endl;
    int addxlinkbonds = 0;
    int bondcounter = 0;
    for (int i=1; i<=lmpi.bi.nbonds;i++){
        //cout << "bond i: " << i << endl;
        int firstatomnum = lmpi.bi.bondfirst[i-1];
        int secondatomnum = lmpi.bi.bondsecond[i-1];
        //cout << " firstatomnum " << firstatomnum << " secondatomnum " << secondatomnum << " newAtoms " << newAtoms << endl;
        //find the indices 
        int firstatype = lmpi.ai.atomtype[firstatomnum-1];
        int secondatype = lmpi.ai.atomtype[secondatomnum-1];
        //cout << " firstatype " << firstatype << " secondatype " << secondatype << endl;
   
        //add edge (calculate bondlength, addedge and updatevertexedges only if xlink type is not involved
        if (find(xlinktypes.begin(),xlinktypes.end(),firstatype) == xlinktypes.end() && find(xlinktypes.begin(),xlinktypes.end(),secondatype) == xlinktypes.end() && (firstatype != solventType) && (secondatype != solventType)){
            gspace.addEdge(firstatomnum-1, secondatomnum-1); 
            //cout << "first " << firstatomnum-1 << " second " << secondatomnum-1  << endl;
            bondcounter++;
            ;
        }

        else {
            //cout << "else" << endl;
            //can use the molecule number of crosslinker, the new nodes were added so new atomnumbers are lmpi.ai.natoms + # of crosslinker, where there is a unique molecule number associated with each # of crosslinker in foundxmol
            //first and second are known (original atomnumbers-> at least one is a xlinker but both could also be xlinkers
            //firstatype and secondatype are also know, if both are xlinkers, no need to add  a bond so only add a bond if one is xlinker and one is not
            //cout << " else " << endl;
            if ((find(xlinktypes.begin(),xlinktypes.end(),firstatype) != xlinktypes.end() && find(xlinktypes.begin(),xlinktypes.end(),secondatype) == xlinktypes.end())  && (firstatype != solventType) && (secondatype != solventType)){
	            //cout << "First is xlinker and second is not" << endl;
                //first is xlinker and second is not
                //check mol of first
                int cmol = lmpi.molecule[firstatomnum-1];
	            int cind = find(foundxmol.begin(),foundxmol.end(),cmol) - foundxmol.begin();
	            xlinkfunctionality[cind]++; 
                //cout << "first " << notXlinkAtoms + cind << " second " << secondatomnum-1 << endl;
	            gspace.addEdge(notXlinkAtoms + cind, secondatomnum-1);
                
                addxlinkbonds++;        	
	            tookcarexlinkmol.push_back(cmol);
                bondcounter++;
            }
            else if (find(xlinktypes.begin(),xlinktypes.end(),firstatype) == xlinktypes.end() && find(xlinktypes.begin(),xlinktypes.end(),secondatype) != xlinktypes.end()  && (firstatype != solventType) && (secondatype != solventType)){
	            //cout << "Second is xlinker and first is not" << endl;
                //second is xlinker and first is not
                int cmol = lmpi.molecule[secondatomnum-1];
                
                //check if cmol is in tookcarexlinkmol
	            int cind = find(foundxmol.begin(),foundxmol.end(),cmol) - foundxmol.begin();
	            xlinkfunctionality[cind]++;
	            addxlinkbonds++;
                //cout << "first " << firstatomnum -1  << " second " << notXlinkAtoms+cind << endl;
	            gspace.addEdge(firstatomnum-1, notXlinkAtoms + cind);
               
                tookcarexlinkmol.push_back(cmol);
                bondcounter++;

            } 

  

        }

    }

    gspace.edges = bondcounter;
    gspace.masses = masses;
    gspace.xLinker = xLinker;
    gspace.actualAinds = newActualAtominds;
    gspace.atomtypes = newatomtypes;
    cout << "total bonds " << bondcounter << endl;
       
    return gspace;
}
   
struct Results{
    int sitesB; 
    int sitesE; 
    int maxE;
    int danglingEnds;
    bool bigEnough;
    int nClusters;
    int nCluster1;
    int nCluster2;
    double massCluster1;
    double massCluster2; 
    int loopsC;
    double avgloopsC;
    int loopsE;
    int loopsB;
    int CburningTime;
    int CshortestPath;
    int BburningTime;
    int BshortestPath;
    int cuttingsE;
    int gAtoms;
    vector <int> provenP1;
    vector <int> provenP2;
    vector <int> eBackAtomInds;
    vector <int> bBackAtomInds;
    vector <int> LclusterAInds;
    vector <double> AclusterMasses;
    double totalMass;
    bool infinite;
    int Cedges;
    int dPrimaryEnds;
    int dSecondaryEnds;
    int dXlinkEnds;
    int nActiveLinks;
    int nActiveStrands;
    double avgActiveLen;
    double avgActiveMass;
    int aXlinks;
    int eActiveChains;
    int eActiveBeads;
    vector <int> elasticBeads;




};



Results analyzeOne(string filename, vector <int> provenP1, vector <int> provenP2, int pastLoops) 
{ 
    // add edges 

    Lmpsabinfo lmpi = read_finalstructwithbondinfo(filename, 100000);
    Results results;

    //get new natoms and nbonds based on crosslinker links
    vector <int> xlinktypes = {22,23,24,25,26};
    GraphC gspace = consolidateCrosslinker(lmpi, xlinktypes, 504.6, 21);
    results.gAtoms = gspace.graph.size();

    GraphC::ClustersInfo clusters = gspace.connectedComponents();
    cout << "finished connectedComponents" << endl;
    cout << "clusters.clusterSizes.size() " << clusters.clusterSizes.size() << " clusters.Cmasses.size() " << clusters.Cmasses.size() << endl;
    //find largest cluster 
    int largestCluster = *max_element(clusters.clusterSizes.begin(), clusters.clusterSizes.end());
    cout << "largestCluster Size " << largestCluster  << endl;
    int Lind = find(clusters.clusterSizes.begin(), clusters.clusterSizes.end(), largestCluster) - clusters.clusterSizes.begin();

    cout << "Lind" << Lind << endl;
    int secondLargest = 0;
    int sInd;
    results.totalMass = 0;
    for (int i =0; i < clusters.clusterSizes.size(); i++){
        cout << " clusterSize " << clusters.clusterSizes[i] << " mass " << clusters.Cmasses[i] << endl;
        if (clusters.clusterSizes[i] > secondLargest && clusters.clusterSizes[i] < largestCluster){
            secondLargest = clusters.clusterSizes[i];
            sInd = i;
        }
        results.totalMass += clusters.Cmasses[i];
	results.AclusterMasses.push_back(clusters.Cmasses[i]);
    }
        
    results.nClusters = clusters.clusters;
    results.nCluster1 = largestCluster;
    results.nCluster2 = secondLargest;
    results.massCluster1 = clusters.Cmasses[Lind];
    results.massCluster2 = clusters.Cmasses[sInd];


    GraphC network = clusters.clusterGraphs[Lind];
    for (int i = 0; i < network.atominds.size(); i++){
        int ai = network.atominds[i];
        for (int j = 0; j < gspace.actualAinds[ai].size(); j++){
            results.LclusterAInds.push_back(gspace.actualAinds[ai][j]); 
        }
    }


    results.Cedges = network.edges;
    Burning networkB;
    networkB.cluster = network.graph;
    networkB.xcoords = network.xcoords;
    networkB.ycoords = network.ycoords;
    networkB.zcoords = network.zcoords;
    
    networkB.edges = network.edges;
    networkB.clusterSites = largestCluster;
    
    //determine if cluster is big enough 
    int maxSitesE = 0;
    GraphC networkMid = network;
    DanglingInfo dInfo = networkMid.removeDanglingEnds();
    results.dPrimaryEnds = dInfo.primaryEnds;
    results.dSecondaryEnds = dInfo.secondaryEnds;
    results.dXlinkEnds = dInfo.xlinkEnds;
    results.danglingEnds = dInfo.danglingEnds;
    
    ClusterSize czp1p2 = networkMid.checkSize(lmpi.xlo, lmpi.xhi, lmpi.ylo, lmpi.yhi, lmpi.zlo, lmpi.zhi);
    results.bigEnough = czp1p2.bigenough;
    cout << "bigEnough " << results.bigEnough << endl;
    results.infinite = false;
    
    



    if (czp1p2.bigenough){
        //find elastically active Xlinks, 3 bonds not counting dangling ends;
        vector <int> eaXlinks; 
        vector <int> eaStrands;
        vector <int> eaBeads; 
        vector <double> eaMass;
        vector <int> allXlinksC; 
        int activeEchains = 0;
        results.aXlinks = 0;

        for (int i = 0; i < networkMid.atominds.size(); i++){
            int aind = networkMid.atominds[i];
            if (networkMid.xLinker[aind]){
                int nbonds = networkMid.graph[aind].size();
                if (nbonds >= 3){
                    eaXlinks.push_back(aind);
                    //cout << "active xlinker " << aind << " " << eaXlinks.size() << endl;
                    results.aXlinks++;
                }
                allXlinksC.push_back(aind);

            }
        }

        //try all elastically active Xlinks
        P1P2list p1p2B = networkMid.boundaryBonds((lmpi.xhi - lmpi.xlo)/2.0);
        P1P2list p1p2;

        for (int i = 0; i < eaXlinks.size(); i++){
            for (int j = i+1; j < eaXlinks.size(); j++){
                p1p2.P1.push_back(eaXlinks[i]);
                p1p2.P2.push_back(eaXlinks[j]);
            }
        }
    
        

        int maxT = 0;
        int maxInd = 0;
        p1p2.P1.insert(p1p2.P1.end(), provenP1.begin(), provenP1.end());
        p1p2.P2.insert(p1p2.P2.end(), provenP2.begin(), provenP2.end());
        double maxDist = (lmpi.xhi - lmpi.xlo)/2.0;
        int maxCloops = 0;
        vector <int> indices;
        vector <int> clusterLoops;
        vector <int> clusterEsites;
        vector <int> clusterShortest;
        vector <int> clusterBsites;

        //find maxCloops first
        for (int i = 0; i < p1p2B.P1.size(); i++){
            vector <int> burnVal (network.graph.size(), 0);
            vector <int> empty;
            networkB.burnVal = burnVal;
            networkB.cluster = network.graph;
            networkB.backbone = empty;
            networkB.elasticBackbone = empty;
            networkB.Eedgelist1 = empty;
            networkB.Eedgelist2 = empty;
            networkB.loopsClosed = empty;
            networkB.Bedgelist1 = empty;
            networkB.Bedgelist2 = empty;
            networkB.clusterSites = largestCluster;

            int t = networkB.firstBurn(p1p2B.P1[i], p1p2B.P2[i], maxDist);
             //cout << "loops" << networkB.loops << endl;
             //cout << "done first burn" << endl;
            double buffer = 0;
           
            //if (networkB.burnedSites < networkB.clusterSites){
                //continue;
                
            //}
            //if (networkB.loops < maxCloops || networkB.loops < 1){
                //continue;
            //}
            if (networkB.loops > maxCloops){

                maxCloops = networkB.loops;
                results.infinite = true;
            }

        }
        for (int i = 0; i < p1p2.P1.size(); i++){
            //cout << " P1 " << p1p2.P1[i] << " P2 " << p1p2.P2[i] << endl;
            
            vector <int> burnVal (network.graph.size(), 0);
            vector <int> empty;
            networkB.burnVal = burnVal;
            networkB.cluster = network.graph;
            networkB.backbone = empty;
            networkB.elasticBackbone = empty;
            networkB.Eedgelist1 = empty;
            networkB.Eedgelist2 = empty;
            networkB.loopsClosed = empty;
            networkB.Bedgelist1 = empty;
            networkB.Bedgelist2 = empty;
            networkB.clusterSites = largestCluster;
            
            //networkB.shortestX = 0;
            //networkB.shortestY = 0;
            //networkB.shortestZ = 0;
            //cout << "first burn" << endl;
             int t = networkB.firstBurn(p1p2.P1[i], p1p2.P2[i], maxDist);
             //cout << "shortestPath " << networkB.shortestPath << endl;
             //cout << "loops" << networkB.loops << endl;
             //cout << "done first burn" << endl;
            double buffer = 0;
           
            //if (networkB.burnedSites < networkB.clusterSites){
                //continue;
                
            //}
            //if (networkB.loops < maxCloops || networkB.loops < 1){
                //continue;
            //}
            if (networkB.loops >= maxCloops){
                provenP1.push_back(p1p2.P1[i]);
                provenP2.push_back(p1p2.P2[i]);
                maxCloops = networkB.loops;
                results.infinite = true;
            }
            
            /*
            else {
                 cout << "shortestX " << networkB.shortestX << " shortestY " << networkB.shortestY << " shortestZ " << networkB.shortestZ << endl;
            
                if (networkB.shortestX < 2*(lmpi.xhi-lmpi.xlo - buffer)  && networkB.shortestY < 2*(lmpi.yhi - lmpi.ylo - buffer) && networkB.shortestZ < 2*(lmpi.zhi - lmpi.zlo - buffer)){
                    cout << "shortestX " << networkB.shortestX << " shortestY " << networkB.shortestY << " shortestZ " << networkB.shortestZ << endl;
            
                    cout << "skipping this one, not infinite distance " << endl;
                    continue;
                }
                else{
                    cout << "infinite distance" << endl;
                    results.infinite = true;
                } 
            }
            */
            //cout << "burningTime " << networkB.burningTime << " shortestPath " << networkB.shortestPath << endl;
            //cout << "burnedSites " << networkB.burnedSites << " totalSites " << networkB.clusterSites << endl;
            //cout << "loops " << networkB.loops << endl;
            //cout << "second burn" << endl;
            networkB.burn2Val = burnVal;
            t = networkB.secondBurn(p1p2.P1[i], p1p2.P2[i]);
            //cout << "sitesE " << networkB.sitesE << endl;
            bool newE = false;
            double eMass = 0;
            for (int i = 0; i < networkB.elasticBackbone.size(); i++){
                if (find(eaBeads.begin(), eaBeads.end(), networkB.elasticBackbone[i]) == eaBeads.end() && networkB.elasticBackbone.size() > 1){
                    newE = true;
                    eaBeads.push_back(networkB.elasticBackbone[i]);
                }
                eMass += network.masses[networkB.elasticBackbone[i]];
            }
            if (networkB.sitesE > 1){
                activeEchains++;
                eaStrands.push_back(networkB.sitesE);
                eaMass.push_back(eMass);
            }
        
            
            //cout << "sitesE " << networkB.sitesE << " loopsE " << networkB.loopsE << " maxE " << networkB.maxE << " cuttingsE " << networkB.cuttingsE << endl;
            /*
            if (networkB.sitesE >= maxSitesE){
                maxSitesE = networkB.sitesE;
            }
            else {
                continue;
            }
            */
            //cout << "third burn" << endl;
            networkB.burn3Val = burnVal;
            t = networkB.thirdBurn(p1p2.P1[i], p1p2.P2[i]);
           //cout << "sitesB " << networkB.sitesB << endl;
            
            GraphC backboneG(network.graph.size());
            backboneG.atominds = networkB.backbone;
            
            //cout << " backboneG.atominds.size() " << backboneG.atominds.size() << endl;
            for (int ii = 0; ii < networkB.Bedgelist1.size(); ii++){
                backboneG.addEdge(networkB.Bedgelist1[ii], networkB.Bedgelist2[ii]);
            }

        
            backboneG.xcoords = network.xcoords;
            backboneG.ycoords = network.ycoords;
            backboneG.zcoords = network.zcoords;
            backboneG.masses = network.masses;
            backboneG.xLinker = network.xLinker;

            Burning backbone;
            backbone.cluster = backboneG.graph;
            backbone.xcoords = network.xcoords;
            backbone.ycoords = network.ycoords;
            backbone.zcoords = network.zcoords;
            backbone.edges = backboneG.edges;
            backbone.clusterSites = backboneG.atominds.size();

            

            backbone.burnVal = burnVal;
             //cout << "first burn again" << endl;
            backbone.firstBurnComplete(p1p2.P1[i], p1p2.P2[i]);
           
            
            //cout << "burnedSites " << backbone.burnedSites << " totalSites " << backbone.clusterSites << endl;
            //cout << "loops " << backbone.loops << endl;
            
            indices.push_back(i);
            clusterLoops.push_back(networkB.loops);
            clusterEsites.push_back(networkB.sitesE);
            clusterShortest.push_back(networkB.shortestPath);
            clusterBsites.push_back(networkB.sitesB);
            //if (backbone.burnedSites > maxT){
            //    maxT = backbone.burnedSites;
            //    maxInd = i;
           //}
        }

       
        if (maxCloops >= pastLoops){
            results.loopsC = maxCloops;
            pastLoops = results.loopsC;
        }
        else {
            results.loopsC = pastLoops;
        }
        results.avgActiveLen = 0;
        results.avgActiveMass = 0;
        results.provenP1 = provenP1;
        results.provenP2 = provenP2;
        results.eActiveChains = activeEchains;
        results.eActiveBeads = eaBeads.size();
        results.elasticBeads = eaBeads;
        for (int i = 0; i < eaStrands.size(); i++){
            results.avgActiveLen += eaStrands[i];
            results.avgActiveMass += eaMass[i];
        }
        results.avgActiveLen = results.avgActiveLen/results.eActiveChains;
        results.avgActiveMass = results.avgActiveMass/results.eActiveChains;

        }

        

   
    return results;
    
} 

int main(){
	string path = "Example";
        string postfix = ".data";
        vector<string> datafiles;
        for (int i = 1; i <= 1; i++){
        string name = path + postfix;
        datafiles.push_back(name);
        }
	string outfile = "burning_info_Example.txt";
	string atomfile = "atom_info_Example.txt";
    string Lclusterfile = "largest_cluster_Example.txt";
	ofstream mf;
	ofstream af;
    ofstream sf;
    vector <int> provenP1;
    vector <int> provenP2;
    int pastLoops = 0;
    for (int i = 0; i < datafiles.size(); i++){
        Results results = analyzeOne(datafiles[i],provenP1, provenP2, pastLoops);
        pastLoops = results.loopsC; 
	cout << datafiles[i] << endl;
	mf.open(outfile, ofstream::app);
        if (!results.bigEnough || !results.infinite){
            
            mf << "Filename nClusters gAtoms massCluster1 nCluster1 massCluster2 nCluster2 bigEnough Cedges TotalMass d1Ends d2Ends dxLinkEnds danglingEnds" << endl;
            mf << datafiles[i] << " " << results.nClusters << " " << results.massCluster1 << " " << results.nCluster1 << " " << results.massCluster2 << " " << results.nCluster2 << " " << results.bigEnough << " " << results.Cedges << " " << fixed << setprecision(2) << results.totalMass << " "<<  results.dPrimaryEnds << " " << results.dSecondaryEnds << " " << results.dXlinkEnds << " "<<  results.danglingEnds << endl;
        }
        else if (results.bigEnough && results.infinite){
           

            mf << "Filename nClusters gAtoms massCluster1 nCluster1 massCluster2 nCluster2 bigEnough Cedges TotalMass d1Ends d2Ends dxLinkEnds danglingEnds Cloops eActiveJuncs eActiveChains eActiveBeads avgActiveLen avgActiveMass" << endl;
            mf << datafiles[i] << " " << results.nClusters << " " << results.gAtoms << " " << results.massCluster1 << " " << results.nCluster1 << " " << results.massCluster2 << " " << results.nCluster2 << " " << results.bigEnough << " " << results.Cedges << " " <<  fixed << setprecision(2) << results.totalMass << " "<<  results.dPrimaryEnds << " " << results.dSecondaryEnds << " " << results.dXlinkEnds << " " << results.danglingEnds << " ";
           
            mf << results.loopsC << " " << results.aXlinks << " " << results.eActiveChains << " " << results.eActiveBeads << " " << results.avgActiveLen << " " << results.avgActiveMass << endl;
        }    
        
    	mf.flush();
        mf.close();

  
    
    if (results.bigEnough && results.infinite){
        af.open(atomfile, ofstream::app);
        af << datafiles[i] << " Elastic_Backbone: ";
        //cout << datafiles[i] << " Elastic_Backbone: ";
        for (int i = 0; i < results.elasticBeads.size(); i++){
            af << results.elasticBeads[i] << " ";
            //cout << results.eBackAtomInds[i] << " ";
        }
        af << "" << endl;
        //cout << " " << endl;
        //if (results.bigEnough){
        //af << datafiles[i] << " P1 " << results.provenP1[provenP1.size()-1] << " P2 "  << results.provenP2[provenP2.size()-1] << endl;
        //}

        af.flush();
        af.close();

    }
        sf.open(Lclusterfile, ofstream::app);
        sf << datafiles[i] << " Largest Cluster: ";
        for (int i = 0; i <results.LclusterAInds.size(); i++){
            sf << results.LclusterAInds[i] << " ";
        }
        sf << "" << endl;
        sf << datafiles[i] << " Cluster Masses: ";
        for (int i = 0; i < results.AclusterMasses.size(); i++){
            sf << results.AclusterMasses[i] << " ";
        }
        sf << "" << endl;
        
        sf.flush();
        sf.close();

    }    
    
    return 0;
}

