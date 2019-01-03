#ifndef _TREENODE_H_
#define _TREENODE_H_

#include "./Utilities.h"


#define G 6.673e-11
#define TOL 1e-12

class TreeNode
{
public:
	TreeNode(double x1, double x2, double y1, double y2, double z1, double z2);
	~TreeNode();
	void addBody(const body_t& body);
	bool containsBody(const body_t& body) const;
	void diagnostics(int level, int& maxlevel, int& nbodies, int& nnodes) const;
    void computeCoM();
    void getCoM(body_t& com) const;
    void computeForceOnBody(const body_t& body, double theta, double3_t& F) const;
    int getBodyCount() const {return number_of_bodies_;}
    void prune();
    void LETBodies(const domain_t& domain, double theta, std::vector<body_t>& bodies);
    
private:
	
	void spawnChildren();
	
	int number_of_bodies_;
	body_t body_, com_;
	double xmin_, xmax_, ymin_, ymax_, zmin_, zmax_;	//bounding box
	
	TreeNode* children_[8];
};

#endif