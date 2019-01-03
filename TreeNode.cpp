#include <iostream>
#include "TreeNode.h"
#include <math.h>
#include <cassert>

using namespace std;

TreeNode::TreeNode(double x1, double x2, double y1, double y2, double z1, double z2)
{
    xmin_ = x1;
    ymin_ = y1;
    zmin_ = z1;
    xmax_ = x2;
    ymax_ = y2;
    zmax_ = z2;
    
    number_of_bodies_ = 0;
    for (unsigned int ichild = 0; ichild < 8; ichild++) children_[ichild] = NULL;
}

TreeNode::~TreeNode()
{
    for (unsigned int ichild = 0; ichild < 8; ichild++)
    {
        if (children_[ichild] != NULL) delete children_[ichild];
        children_[ichild] = NULL;
    }
}

void TreeNode::spawnChildren()
{
    double dx = 0.5*(xmax_ - xmin_);
    double dy = 0.5*(ymax_ - ymin_);
    double dz = 0.5*(zmax_ - zmin_);
    
    for (unsigned int ichild = 0; ichild < 8; ichild++)
    {
        assert(children_[ichild] == NULL);
    }
    
    children_[0] = new TreeNode(xmin_, xmin_ + dx, ymin_, ymin_ + dy, zmin_, zmin_ + dz);
    children_[1] = new TreeNode(xmin_ + dx, xmax_, ymin_, ymin_ + dy, zmin_, zmin_ + dz);
    children_[2] = new TreeNode(xmin_, xmin_ + dx, ymin_ + dy, ymax_, zmin_, zmin_ + dz);
    children_[3] = new TreeNode(xmin_ + dx, xmax_, ymin_ + dy, ymax_, zmin_, zmin_ + dz);
    children_[4] = new TreeNode(xmin_, xmin_ + dx, ymin_, ymin_ + dy, zmin_ + dz, zmax_);
    children_[5] = new TreeNode(xmin_ + dx, xmax_, ymin_, ymin_ + dy, zmin_ + dz, zmax_);
    children_[6] = new TreeNode(xmin_, xmin_ + dx, ymin_ + dy, ymax_, zmin_ + dz, zmax_);
    children_[7] = new TreeNode(xmin_ + dx, xmax_, ymin_ + dy, ymax_, zmin_ + dz, zmax_);
}

bool TreeNode::containsBody(const body_t& body) const
{
    if (body.r[0] >= xmin_ && body.r[0] <= xmax_ && body.r[1] >= ymin_ && body.r[1] <= ymax_ && body.r[2] >= zmin_ && body.r[2] <= zmax_) return true;
    return false;
}


void TreeNode::addBody(const body_t& body)
{
    if (number_of_bodies_ == 0)
    {
        body_ = body;
        number_of_bodies_ = 1;
    }
    else if (number_of_bodies_ == 1)
    {
        spawnChildren();

        bool found_child = false;
        for (unsigned int ichild = 0; ichild < 8; ichild++)
        {
            if (children_[ichild]->containsBody(body_))
            {
                children_[ichild]->addBody(body_);
                ichild = 8; //break loop
                found_child = true;
            }
        }
        if (!found_child)
        {
            std::cerr << "Could not find a child for body ( " << body_.r[0] << ", " << body_.r[1] << ", " << body_.r[2] << " ) " << std::endl;
            std::cerr << "In box ( " << xmin_ << ", " << ymin_ << ", " << zmin_ << " ) to ( " << xmax_ << ", " << ymax_ << ", " << zmax_ << " ) " << std::endl;
        }

        found_child = false;
        for (unsigned int ichild = 0; ichild < 8; ichild++)
        {
            if (children_[ichild]->containsBody(body))
            {
                children_[ichild]->addBody(body);
                ichild = 8; //break loop
                found_child = true;
            }
        }
        if (!found_child)
        {
            std::cerr << "Could not find a child for body ( " << body_.r[0] << ", " << body_.r[1] << ", " << body_.r[2] << " ) " << std::endl;
            std::cerr << "In box ( " << xmin_ << ", " << ymin_ << ", " << zmin_ << " ) to ( " << xmax_ << ", " << ymax_ << ", " << zmax_ << " ) " << std::endl;
        }
        number_of_bodies_++;
    }
    else
    {
        bool found_child = false;
        for (unsigned int ichild = 0; ichild < 8; ichild++)
        {
            if (children_[ichild]->containsBody(body))
            {
                children_[ichild]->addBody(body);
                ichild = 8;
                found_child = true;
            }
        }
        number_of_bodies_++;
        
        if (!found_child)
        {
            std::cerr << "Could not find a child for body ( " << body_.r[0] << ", " << body_.r[1] << ", " << body_.r[2] << " ) " << std::endl;
            std::cerr << "In box ( " << xmin_ << ", " << ymin_ << ", " << zmin_ << " ) to ( " << xmax_ << ", " << ymax_ << ", " << zmax_ << " ) " << std::endl;
        }

    }
}

void TreeNode::prune()
{
    for (unsigned int ichild = 0; ichild < 8; ichild++)
    {
        if (children_[ichild] != NULL)
        {
            if (children_[ichild] -> getBodyCount() < 1)
            {
                delete children_[ichild];
                children_[ichild] = NULL;
            }
            else
            {
                children_[ichild]->prune();
            }
        }
    }
}


void TreeNode::diagnostics(int level, int& maxlevel, int& nbodies, int& nnodes) const
{
    if (level > maxlevel) maxlevel=level;
    nnodes++;

    bool isleaf = true;
    for (unsigned int ichild = 0; ichild < 8; ichild++)
    {
        if (children_[ichild] != NULL)
        {
            children_[ichild]->diagnostics(level+1, maxlevel, nbodies, nnodes);
            isleaf = false;
        }
    }
    if (isleaf)
    {
        nbodies+= number_of_bodies_;
    }
}

void TreeNode::getCoM(body_t& com) const
{
    com = com_;
}




void TreeNode::computeCoM()
{
    if (number_of_bodies_ == 0)
    {
        com_.r[0] = 0.5*(xmax_ + xmin_);
        com_.r[1] = 0.5*(ymax_ + ymin_);
        com_.r[2] = 0.5*(zmax_ + zmin_);
        com_.m = 0.0;
    }
    else if (number_of_bodies_ == 1)
    {
        com_ = body_;
    }
    else
    {
		com_.r[0] = 0.0;
		com_.r[1] = 0.0;
		com_.r[2] = 0.0;
		com_.m = 0.0;
		body_t tmp;
		bool all_null = true;
		for (int ichild = 0; ichild < 8; ichild++) {
			if (children_[ichild] != NULL)
			{
				all_null = false;
				children_[ichild]->computeCoM();
				children_[ichild]->getCoM(tmp);
				com_.m += tmp.m;
				com_.r[0] += tmp.r[0] * tmp.m;
				com_.r[1] += tmp.r[1] * tmp.m;
				com_.r[2] += tmp.r[2] * tmp.m;
			}
		}
		assert(all_null == false);
		com_.r[0] /= com_.m;
		com_.r[1] /= com_.m;
		com_.r[2] /= com_.m;
    }
}





void TreeNode::computeForceOnBody(const body_t& body, double theta, double3_t& F) const
{
	
    if (number_of_bodies_ > 1)
    {
		double D = xmax_ - xmin_;
		
		double R = sqrt((body.r[0] - com_.r[0])*(body.r[0] - com_.r[0]) + (body.r[1] - com_.r[1])*(body.r[1] - com_.r[1]) + (body.r[2] - com_.r[2])*(body.r[2] - com_.r[2]));
		if (D/R< theta)
		{
			double tmp = G * body.m*com_.m / (R*R*R);
			F.r[0] -= tmp * (body.r[0] - com_.r[0]);
			F.r[1] -= tmp * (body.r[1] - com_.r[1]);
			F.r[2] -= tmp * (body.r[2] - com_.r[2]);
			return;

		}
		else {
			for (int ichild = 0; ichild < 8; ichild++) {
				if (children_[ichild] != NULL)
				{
					children_[ichild]->computeForceOnBody(body,theta,F);
					
				}
			}


		}
 
    }
    else if (number_of_bodies_ == 1)
    {
        //note com is equal to body in this case by our convention.
        double R = sqrt((body.r[0] - com_.r[0])*(body.r[0] - com_.r[0]) + (body.r[1] - com_.r[1])*(body.r[1] - com_.r[1]) + (body.r[2] - com_.r[2])*(body.r[2] - com_.r[2]));
        
        if (R > 1e-14)
        {
            double tmp = G*body.m*com_.m/(R*R*R);
            F.r[0] -= tmp*(body.r[0] - com_.r[0]);
            F.r[1] -= tmp*(body.r[1] - com_.r[1]);
            F.r[2] -= tmp*(body.r[2] - com_.r[2]);
        }
    }
    else return;
}


 int a = 0;
 int b = 0;
void TreeNode::LETBodies(const domain_t& domain, double theta, std::vector<body_t>& bodies)
{
	
    if (number_of_bodies_ == 0) return;
    else if (number_of_bodies_ == 1)
    {
        bodies.push_back(body_);
		b++;
        return;
    }
    
    //no need for an else here as the above conditions return.
    
    //------------------------------------------------------------------------
    //The following code will compute the minimum distance between the current
    //node and the given domain.
    //------------------------------------------------------------------------
    
    double R2 = 0;
        
    if (domain.max[0] < xmin_) R2 += (domain.max[0] - xmin_)*(domain.max[0] - xmin_);
    else if (domain.min[0] > xmax_) R2 += (domain.min[0] - xmax_)*(domain.min[0] - xmax_);
        
    if (domain.max[1] < ymin_) R2 += (domain.max[1] - ymin_)*(domain.max[1] - ymin_);
    else if (domain.min[1] > ymax_) R2 += (domain.min[1] - ymax_)*(domain.min[1] - ymax_);
        
    if (domain.max[2] < zmin_) R2 += (domain.max[2] - zmin_)*(domain.max[2] - zmin_);
    else if (domain.min[2] > zmax_) R2 += (domain.min[2] - zmax_)*(domain.min[2] - zmax_);
        
    double R = sqrt(R2);
	double D = xmax_-xmin_;
	//cout << "R= " << R << " D=" << D <<"  "<< endl;
	
	if (D / R <= theta)
	{
		//cout << "R= "<<R<<" D="<<D << endl;
		//cout << "xmin_: " << xmin_ << endl;
		//cout <<"R: "<< D/R << endl;
		//cout << "R= " << R << " D=" << D << "  D/R=" << D / R << endl;
		computeCoM();
		bodies.push_back(com_);
		a++;
	}
	else {
	
		//bool all_null = true;
		for (int ichild = 0; ichild < 8; ichild++) {
			if (children_[ichild] != NULL)
			{
				//all_null = false;
				children_[ichild]->LETBodies(domain, theta, bodies);

			}
		}
		
	
	}
	//cout << "sendCOM leaf-node size = " << b << endl;
	//cout <<"sendCOM size = "<< a << endl;
        
}

