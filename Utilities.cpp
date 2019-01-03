#include <iostream>
#include <vector>
#include <math.h>
#include <cassert>
#include <mpi.h>
#include <queue>
#include "./Utilities.h"
#include <sstream>
#include <fstream>


using namespace std;

//-----------------------------------------------------
// Usual Parallel Partitioning Code
//-----------------------------------------------------
void parallelRange(int globalstart, int globalstop, int irank, int nproc, int& localstart, int& localstop, int& localcount)
{
	int nvals = globalstop - globalstart + 1;
	int divisor = nvals / nproc;
	int remainder = nvals % nproc;
	int offset;
	if (irank < remainder) offset = irank;
	else offset = remainder;

	localstart = irank * divisor + globalstart + offset;
	localstop = localstart + divisor - 1;
	if (remainder > irank) localstop += 1;
	localcount = localstop - localstart + 1;
}

//-----------------------------------------------------
// Extract a single coordinate list from a double3_t list
//-----------------------------------------------------
std::vector<double> getSingleCoordinateListFromPoints(const std::vector<double3_t>& points, int dim)
{
	std::vector<double> coordinate_list(points.size());
	if (dim > 3 || dim < 0)
	{
		std::cerr << "Requested dimension " << dim << " is out of bounds at line " << __LINE__ << " of file " << __FILE__ << " for function " << __FUNCTION__ << std::endl;
		return coordinate_list;
	}

	for (unsigned int ipoint = 0; ipoint < points.size(); ipoint++)
	{
		coordinate_list[ipoint] = points[ipoint].r[dim];
	}

	//Want to see the list?
	/*
	for (unsigned int ipoint = 0; ipoint < coordinate_list.size(); ipoint++)
	{
		std::cout << coordinate_list[ipoint] << std::endl;
	}
	*/
	return coordinate_list;
}


//-----------------------------------------------------
// Implementation of ORB.
//
// Inputs:
//
// P the number of domains to produce.
//
// points a list of type double3_t to partition.
//
// Output:
// points_in_orb_domains is a vector of vectors where
// points_in_orb_domains[iproc] stores the local element
// indexes of all processors in subdomain iproc (i.e.,
// destined for iproc).
//-----------------------------------------------------
void ORB(int P, const std::vector<double3_t>& points, std::vector<std::vector<int> >& points_in_orb_domains)
{
	int rank, nproc;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	points_in_orb_domains.resize(0);

	std::vector<int> point_indexes_in_domain(points.size());
	for (unsigned int ipoint = 0; ipoint < points.size(); ipoint++)
	{
		point_indexes_in_domain[ipoint] = ipoint;
	}



	int dim = 0;

	std::queue<double> weight_queue;
	std::queue<std::vector<int> > points_in_domain_queue;
	std::queue<int> dim_queue;

	weight_queue.push(P);
	dim_queue.push(dim);
	points_in_domain_queue.push(point_indexes_in_domain);

	std::cout << "Rank " << rank << " initial point size " << std::endl;



	while (!weight_queue.empty())
	{
		int weight = weight_queue.front();
		weight_queue.pop();

		dim = dim_queue.front();
		dim_queue.pop();

		std::vector<int> domain_point_indexes = points_in_domain_queue.front();
		points_in_domain_queue.pop();

		int weight00 = (int)floor(weight / 2.0);
		int weight01 = (int)ceil(weight / 2.0);

		std::vector<double3_t> domain_points(domain_point_indexes.size());
		for (unsigned int ipoint = 0; ipoint < domain_point_indexes.size(); ipoint++)
		{
			domain_points[ipoint] = points[domain_point_indexes[ipoint]];
		}

		std::vector<double> coords = getSingleCoordinateListFromPoints(domain_points, dim);
		std::vector<double> sorted_coords;

		parallelBucketSort(coords, sorted_coords);

		int n_local_points_in_domain = domain_point_indexes.size();

		int M0;

		MPI_Allreduce(&n_local_points_in_domain, &M0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		if (rank == 0) std::cout << "Total points across all procs at this ORB level = " << M0 << std::endl;

		int M00 = (int)((double)M0*(double)weight00 / (double)weight);

		//get the M00th entry of the list of sorted coordinates that is distributed across all processors

//        usleep(30000*rank);
//        std::cerr << "Rank " << rank << " has M0 = " << M0 << " and M00 = " << M00 << std::endl;

		std::vector<int> n_local_sorted_coords_vec(nproc, sorted_coords.size());
		std::vector<int> processor_sorted_coords_counts(nproc, 1);

		MPI_Alltoall(&n_local_sorted_coords_vec[0], 1, MPI_INT, &processor_sorted_coords_counts[0], 1, MPI_INT, MPI_COMM_WORLD);


		//        usleep(30000*rank);


		//        std::cerr << "Rank " << rank << " has the following sorted point distribution:  " << std::endl;
		//        for (int iproc = 0; iproc < nproc; iproc++)
		//        {
		//            std::cerr << processor_sorted_coords_counts[iproc] << std::endl;
		//        }

		int proc_with_weighted_median = -1;
		int index_of_weighted_median_on_proc = -1;
		int running_sum = 0;
		for (int iproc = 0; iproc < nproc; iproc++)
		{
			if (running_sum + processor_sorted_coords_counts[iproc] >= M00)
			{
				proc_with_weighted_median = iproc;
				index_of_weighted_median_on_proc = M00 - running_sum;
				break;
			}
			else
			{
				running_sum += processor_sorted_coords_counts[iproc];
			}
		}

		//        usleep(30000*rank);
		//        std::cerr << "Rank " << rank << " has proc with median = " << proc_with_weighted_median << std::endl;
		//        std::cerr << "Rank " << rank << " has local index of median = " << index_of_weighted_median_on_proc << std::endl;
		//        std::cerr << "Rank " << rank << " has running sum plus median index = " << running_sum +index_of_weighted_median_on_proc << std::endl;



		double local_weighted_median = 0;

		if (rank == proc_with_weighted_median)
		{
			local_weighted_median = sorted_coords[index_of_weighted_median_on_proc];
		}

		double weighted_median;
		MPI_Allreduce(&local_weighted_median, &weighted_median, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		std::vector<int> domain_00_point_indexes(0);
		std::vector<int> domain_01_point_indexes(0);

		for (unsigned int ipoint = 0; ipoint < domain_point_indexes.size(); ipoint++)
		{
			if (points[domain_point_indexes[ipoint]].r[dim] <= weighted_median)
			{
				domain_00_point_indexes.push_back(domain_point_indexes[ipoint]);
			}
			else
			{
				domain_01_point_indexes.push_back(domain_point_indexes[ipoint]);
			}
		}

		assert(domain_00_point_indexes.size() + domain_01_point_indexes.size() == domain_point_indexes.size());


		int n00 = domain_00_point_indexes.size();
		int n01 = domain_01_point_indexes.size();
		int n0 = domain_point_indexes.size();

		int n00_global, n01_global, n0_global;

		MPI_Allreduce(&n0, &n0_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&n00, &n00_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&n01, &n01_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		//        if (rank == 0)
		//        {
		//            std::cout << "Global n0 = " << n0_global << std::endl;
		//            std::cout << "Global n00 = " << n00_global << std::endl;
		//            std::cout << "Global n01 = " << n01_global << std::endl;
		//            std::cout << "Sum = " << n00_global + n01_global << std::endl;
		//        }


		if (weight00 == 1)
		{
			points_in_orb_domains.push_back(domain_00_point_indexes);
		}
		else
		{
			weight_queue.push(weight00);
			dim_queue.push((dim + 1) % 3);
			points_in_domain_queue.push(domain_00_point_indexes);
		}

		if (weight01 == 1)
		{
			points_in_orb_domains.push_back(domain_01_point_indexes);
		}
		else
		{
			weight_queue.push(weight01);
			dim_queue.push((dim + 1) % 3);
			points_in_domain_queue.push(domain_01_point_indexes);
		}


		//       usleep(1000000*rank);
		        
		  //      std::cerr << "On rank " << rank << std::endl;
		    //    std::cerr << "After ORB split there are " << points_in_orb_domains.size() << " completed subdomains " << std::endl;
		      //  std::cerr << "They each have size: " << std::endl;
		   //    for (int idom = 0; idom < points_in_orb_domains.size(); idom++)
		     //  {
		       //     std::cerr << "\t" << points_in_orb_domains[idom].size() << std::endl;
		       //}

	}

}


//-----------------------------------------------------
// Implementation of Parallel Bucket Sort.
//
// Inputs:
//
// The values to sort.
//
// Outputs:
//
// A subset of the sorted values. Note that the entries
// in sorted values are not the same as values to sort
// on any given processor. They are a subset of the
// total set of sorted values where rank 0 will contain
// the lowest sorted numbers.
//-----------------------------------------------------
void parallelBucketSort(const std::vector<double>& values_to_sort, std::vector<double>& sorted_values)
{

	int rank, nproc;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int nbuckets = rank;
	double local_max_val = std::numeric_limits<double>::min();
	double local_min_val = std::numeric_limits<double>::max();

	for (unsigned int ival = 0; ival < values_to_sort.size(); ival++)
	{
		if (values_to_sort[ival] > local_max_val) local_max_val = values_to_sort[ival];
		if (values_to_sort[ival] < local_min_val) local_min_val = values_to_sort[ival];
	}

	double max_val, min_val;

	MPI_Allreduce(&local_max_val, &max_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&local_min_val, &min_val, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	std::vector<double> bucket_starts(nproc);
	std::vector<double> bucket_stops(nproc);

	std::vector<int> value_buckets(values_to_sort.size());

	for (unsigned int ival = 0; ival < value_buckets.size(); ival++)
	{
		value_buckets[ival] = (int)floor(nproc*(values_to_sort[ival] - min_val) / (max_val - min_val));

		if (value_buckets[ival] < 0) value_buckets[ival] = 0;
		else if (value_buckets[ival] >= nproc) value_buckets[ival] = nproc - 1;
		else if (value_buckets[ival] >= 0 && value_buckets[ival] < nproc)
		{
			//everything is ok
		}
		else
		{
			std::cerr << "Something has gone wrong in parallelBucketSort" << std::endl;
			exit(1);
		}
	}

	std::vector<std::vector<double> > values_to_send(nproc, std::vector<double>(0));

	for (unsigned int ival = 0; ival < value_buckets.size(); ival++)
	{
		values_to_send[value_buckets[ival]].push_back(values_to_sort[ival]);
	}


	std::vector<std::vector<double> > received_values;

	MPI_Alltoall_vecvecT(values_to_send, received_values);

	sorted_values.resize(0);
	for (unsigned int iproc = 0; iproc < nproc; iproc++)
	{
		for (unsigned int ival = 0; ival < received_values[iproc].size(); ival++)
		{
			sorted_values.push_back(received_values[iproc][ival]);
		}
	}

	sort(sorted_values.begin(), sorted_values.end());

	//    std::ostringstream converter;
	//    converter << "Sorted_values_rank_" << rank << "_procs_" << nproc << ".txt";
	//    std::ofstream out_to_file(converter.str().c_str());
	//    
	//    for (unsigned int ival = 0; ival < sorted_values.size(); ival++)
	//    {
	//        out_to_file << sorted_values[ival] << std::endl;
	//    }
	//    out_to_file.close();

}




//-----------------------------------------------------
// For a list of bodies, compute the domain that
// contains them (tight bound) in domain.
// Also return the minimum and maximum location over all
// three dimensions in dimmin and dimmax. The global flag
// can be used to make the domains local (false) or
// global (true).
//-----------------------------------------------------
void getPointExtent(const std::vector<body_t>& bodies, domain_t& domain, double& dimmin, double& dimmax, bool global)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    double local_xmin = HUGE_VAL;
    double local_xmax = -HUGE_VAL;
    double local_ymin = HUGE_VAL;
    double local_ymax = -HUGE_VAL;
    double local_zmin = HUGE_VAL;
    double local_zmax = -HUGE_VAL;
    
    for (unsigned int ibody = 0; ibody < bodies.size(); ibody++)
    {
        if (bodies[ibody].r[0] < local_xmin) local_xmin = bodies[ibody].r[0];
        if (bodies[ibody].r[0] > local_xmax) local_xmax = bodies[ibody].r[0];
        if (bodies[ibody].r[1] < local_ymin) local_ymin = bodies[ibody].r[1];
        if (bodies[ibody].r[1] > local_ymax) local_ymax = bodies[ibody].r[1];
        if (bodies[ibody].r[2] < local_zmin) local_zmin = bodies[ibody].r[2];
        if (bodies[ibody].r[2] > local_zmax) local_zmax = bodies[ibody].r[2];
    }
    
    double local_dimmin = min(min(local_xmin,local_ymin),local_zmin);
    double local_dimmax = max(max(local_xmax,local_ymax),local_zmax);
    
    double xmin, xmax, ymin, ymax, zmin, zmax;
    
    if (global == true)
    {
        MPI_Allreduce(&local_xmin, &xmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local_ymin, &ymin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local_zmin, &zmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        
        MPI_Allreduce(&local_xmax, &xmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&local_ymax, &ymax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&local_zmax, &zmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
        MPI_Allreduce(&local_dimmin, &dimmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local_dimmax, &dimmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    else
    {
        xmin = local_xmin;
        xmax = local_xmax;
        ymin = local_ymin;
        ymax = local_ymax;
        zmin = local_zmin;
        zmax = local_zmax;
        dimmin = local_dimmin;
        dimmax = local_dimmax;
    }
    
    domain.min[0] = xmin;
    domain.min[1] = ymin;
    domain.min[2] = zmin;
    
    domain.max[0] = xmax;
    domain.max[1] = ymax;
    domain.max[2] = zmax;
    
}



