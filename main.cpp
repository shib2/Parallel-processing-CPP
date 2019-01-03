#include <iostream>
#include <vector>
#include <mpi.h>
#include <string>
#include <ctime>
#include <cstdlib>
#include "TreeNode.h"
#include "Utilities.h"

using namespace std;
using std::string;


int main(int argc, char** argv)
{
	//----------------------------
   // MPI Setup
   //----------------------------

	int rank, nproc;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//----------------------------
	// Parse Command Line
	//----------------------------

	if (argc < 3)
	{
		if (rank == 0) {
			std::cerr << "Input arguments are not complete. Required arguments are " << std::endl;
			std::cerr << "points per processor" << std::endl;
			std::cerr << "theta" << std::endl;
		}
		MPI_Finalize();
		exit(1);
	}
	int points = atoi(argv[1]);
//	int points = 50000;
	double theta = strtod(argv[2], NULL);
//	double theta = 0.5;
	if (rank == 0) {
		cout << "There are " << points << " points on each processor. " << endl;
		cout << "Theta is " << theta << endl;
	}
	double t_start= (double)clock() / (double)CLOCKS_PER_SEC;
	srand(clock()+rank);
	
	std::vector<body_t> point;
	std::vector<double3_t> recv;
	for (int i = 0; i < points; i++) {
		body_t tmp;
		double3_t tmp1;
		
		tmp.r[0] = (double)rand() / (double)RAND_MAX;
		tmp1.r[0] = tmp.r[0];
		tmp.r[1] = (double)rand() / (double)RAND_MAX;
		tmp1.r[1] = tmp.r[1];
		tmp.r[2] = (double)rand() / (double)RAND_MAX;
		tmp1.r[2] = tmp.r[2];
		double min = 1;
		double max = 1000;

		tmp.m = (min + 1) + (((double)rand()) / (double)RAND_MAX) * (max - (min + 1));
		point.push_back(tmp);
		recv.push_back(tmp1);
	
	}

	if (nproc > 1) {
		std::vector<std::vector<int> > points_in_orb_domains(nproc);


		

		ORB(nproc, recv, points_in_orb_domains);
		std::vector<std::vector<body_t> > sendbuff(nproc);
		std::vector<std::vector<body_t> > recvbuff(nproc);
		MPI_Barrier(MPI_COMM_WORLD);

		for (int iproc = 0; iproc < nproc; iproc++)
		{
			sendbuff[iproc].resize(0);

			for (int j = 0; j < points_in_orb_domains[iproc].size(); j++)
			{

				sendbuff[iproc].push_back(point[points_in_orb_domains[iproc][j]]);
			}
		}

		MPI_Alltoall_vecvecT(sendbuff, recvbuff);
		MPI_Barrier(MPI_COMM_WORLD);
		point.resize(0);

		for (int iproc = 0; iproc < nproc; iproc++)
		{
			for (int i = 0; i < recvbuff[iproc].size(); i++)
			{
				point.push_back(recvbuff[iproc][i]);
			}
		}
	}
	domain_t domain;
	double dimmin;
	double dimmax;
	getPointExtent(point, domain, dimmin, dimmax, true);
	TreeNode * root = new TreeNode(dimmin - TOL, dimmax + TOL, dimmin - TOL, dimmax + TOL, dimmin - TOL, dimmax
		+ TOL);
	for (int i = 0; i < point.size(); i++) {
		root->addBody(point[i]);
	}
	int level = 0;
	int maxlevel = 0;
	int nnodes = 0;
	int numberofbodies = 0;
	root->diagnostics(level, maxlevel, numberofbodies, nnodes);
	std::cout << "I am rank" << rank << ". My ntMax Level = " << maxlevel << std::endl;
	std::cout << "I am rank" << rank << ". My ntN = " << numberofbodies << std::endl;
	std::cout << "I am rank" << rank << ". My ntNode Count = " << nnodes << std::endl;
	
	root->computeCoM();

	if (nproc > 1) {

		std::vector<std::vector<domain_t> > sendbuff1(nproc);
		std::vector<std::vector<domain_t> > recv_domian(nproc);
		for (int i = 0; i < nproc; i++)
		{
			sendbuff1[i].push_back(domain);
		}

		MPI_Alltoall_vecvecT(sendbuff1, recv_domian);
		std::vector<std::vector<body_t> > send_COM(nproc);
		std::vector<std::vector<body_t> > recv_COM(nproc);
		MPI_Barrier(MPI_COMM_WORLD);

		for (int i = 0; i < nproc; i++) {
			if (i != rank)
				root->LETBodies(recv_domian[i][0], theta, send_COM[i]);
		}

		MPI_Alltoall_vecvecT(send_COM, recv_COM);
		MPI_Barrier(MPI_COMM_WORLD);
		
		for (int i = 0; i < nproc; i++) {
			int count = 0;
			if (i != rank) {
				for (int j = 0; j < recv_COM[i].size(); j++) {

					root->addBody(recv_COM[i][count]);
					count++;

				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		root->computeCoM();

	}
	body_t com;
	root->getCoM(com);
	std::cout << "I am rank " << rank << ". Tree Center of Mass = ( " << com.r[0] << ", " << com.r[1] << ", " <<
		com.r[2] << " )" << " with mass " << com.m << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);

	std::vector<double3_t> F(point.size());
	for (int i = 0; i < point.size(); i++) {
		root->computeForceOnBody(point[i], theta, F[i]);
	}
	double Fxmin = fabs(F[0].r[0]);
	double Fxmax = fabs(F[0].r[0]);
	double Fymin = fabs(F[1].r[0]);
	double Fymax = fabs(F[1].r[0]);
	double Fzmin = fabs(F[2].r[0]);
	double Fzmax = fabs(F[2].r[0]);

	for (int i = 1; i < F.size(); i++) {
		if (fabs(F[i].r[0]) < Fxmin)
			Fxmin = fabs(F[i].r[0]);
		if (fabs(F[i].r[0]) > Fxmax)
			Fxmax = fabs(F[i].r[0]);
		if (fabs(F[i].r[1]) < Fymin)
			Fymin = fabs(F[i].r[1]);
		if (fabs(F[i].r[1]) > Fymax)
			Fymax = fabs(F[i].r[1]);
		if (fabs(F[i].r[2]) < Fzmin)
			Fzmin = fabs(F[i].r[2]);
		if (fabs(F[i].r[2]) > Fzmax)
			Fzmax = fabs(F[i].r[2]);

	}
	double Fxmin_ = Fxmin;
	double Fxmax_ = Fxmax;
	double Fymin_ = Fymin;
	double Fymax_ = Fymax;
	double Fzmin_ = Fzmin;
	double Fzmax_ = Fzmax;
	if (nproc > 1) {
		MPI_Allreduce(&Fxmin, &Fxmin_, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&Fymin, &Fymin_, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&Fzmin, &Fzmin_, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&Fxmax, &Fxmax_, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		MPI_Allreduce(&Fymax, &Fymax_, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		MPI_Allreduce(&Fzmax, &Fzmax_, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		
	}

	MPI_Barrier(MPI_COMM_WORLD);
	cout << "I am rank " << rank << ". The Fxmin is " << Fxmin_<<". The Fxmax is " << Fxmax_ << ". The Fymin is " << Fymin_ << ". The Fymax is " << Fymax_ << ". The Fzmin is " << Fzmin_ << ". The Fzmax is " << Fzmax_ << endl;

	
	
	double t_stop = (double)clock() / (double)CLOCKS_PER_SEC;
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0)
	std::cout << "Overall Time : " << (double)(t_stop-t_start) << std::endl;
	MPI_Finalize();
}