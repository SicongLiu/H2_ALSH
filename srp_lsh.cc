#include "headers.h"

// -----------------------------------------------------------------------------
SRP_LSH::SRP_LSH(					// constructor
	int n,								// cardinality of dataset
	int d,								// dimensionality of dataset
	int K,								// number of hash tables
	int L,								// number of hash layers
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_ = n;
	dim_   = d;
	K_     = K;
	L_     = L;
	data_  = data;

	// -------------------------------------------------------------------------
	//  build hash tables (bulkloading)
	// -------------------------------------------------------------------------
	gen_random_vectors();
	bulkload();
}

// -----------------------------------------------------------------------------
SRP_LSH::~SRP_LSH()					// destructor
{
	for(int l=0; l< L_; l++)
	{
		if (proj_[l] != NULL) {
			for (int i = 0; i < K_; ++i)
			{
				delete[] proj_[l][i];
				proj_[l][i] = NULL;
			}
			delete[] proj_[l];
			proj_[l] = NULL;
		}

		if (hash_code_[l] != NULL)
		{
			for (int i = 0; i < n_pts_; ++i)
			{
				delete[] hash_code_[l][i];
				hash_code_[l][i] = NULL;
			}
			delete[] hash_code_[l];
			hash_code_[l] = NULL;
		}
	}

}

// -----------------------------------------------------------------------------
void SRP_LSH::gen_random_vectors()	// generate random projection vectors
{
	proj_ = new float*[L_];
	for(int l=0; l < L_; l++)
	{
		proj_[l] = new float*[K_];
		for (int i = 0; i < K_; ++i)
		{
			proj_[l][i] = new float[dim_];
			for (int j = 0; j < dim_; ++j)
			{
				proj_[l][i][j] = gaussian(0.0f, 1.0f);
			}
		}
	}

}

// -----------------------------------------------------------------------------
void SRP_LSH::bulkload()			// bulkloading
{
	hash_code_ = new bool*[L_];
	for(int l=0; l< L_; l++)
	{
		hash_code_[l] = new bool*[n_pts_];
		for (int i = 0; i < n_pts_; ++i)
		{
			hash_code_[l][i] = new bool[K_];
			get_proj_vector(data_[i], hash_code_[l][i]);
		}
	}

}

// -----------------------------------------------------------------------------
void SRP_LSH::get_proj_vector(		// get vector after random projection
	const float *data,					// input data 
	bool *hash_code)					// hash code of input data (return)
{
	for (int i = 0; i < K_; ++i)
	{
		float sum = calc_inner_product(dim_, proj_[i], data);

		if (sum >= 0) hash_code[i] = true;
		else hash_code[i] = false;
	}
}

// -----------------------------------------------------------------------------
int SRP_LSH::kmc(					// c-k-AMC search
	int   top_k,						// top-k value
	const float *query,					// input query
	MaxK_List *list)					// top-k MC results (return)
{
	bool *mc_query = new bool[K_];
	get_proj_vector(query, mc_query);

	bool *isMarked = new bool[n_pts_];
	int match = 0;
	for (int i = 0; i < n_pts_; ++i)
	{
		int matched_dim = 0;
		for(int l = 0; l < L_; l++)
		{
			for (int j = 0; j < K_; ++j)
			{
				if (hash_code_[i][j] == mc_query[j])
				{
					++matched_dim;
				}
			}

			/**
			 * Check point by Sicong
			 * Stop condition needs to double check
			 * */
			if(matched_dim == K_ && !isMarked[i])
			{
				list->insert(match, i);
				isMarked[i] = true;
				match++;
			}
		}
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] mc_query;
	mc_query = NULL;

	delete[] isMarked;
	isMarked = NULL;
	return 0;
}
