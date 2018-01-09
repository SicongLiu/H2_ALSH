#include "headers.h"

// -----------------------------------------------------------------------------
L2_ALSH2::L2_ALSH2()				// constructor
{
	n_pts_ = -1;
	dim_ = -1;
	m_ = -1;
	U_ = -1.0f;
	appr_ratio_ = -1.0f;
	data_ = NULL;

	M_ = -1.0f;
	l2_alsh2_dim_ = -1;
	l2_alsh2_data_ = NULL;
	
	lsh_ = NULL;
}

// -----------------------------------------------------------------------------
L2_ALSH2::~L2_ALSH2()				// destructor
{
	if (l2_alsh2_data_ != NULL) {
		for (int i = 0; i < n_pts_; i++) {
			delete[] l2_alsh2_data_[i]; l2_alsh2_data_[i] = NULL;
		}
		delete[] l2_alsh2_data_; l2_alsh2_data_ = NULL;
	}

	if (lsh_ != NULL) {
		delete lsh_; lsh_ = NULL;
	}
}

// -----------------------------------------------------------------------------
void L2_ALSH2::init(				// init the parameters
	int n,								// number of data
	int qn,								// number of queries
	int d,								// dimension of data
	int m,								// additional dimension of data
	float U,							// scale factor for data
	float ratio,						// approximation ratio
	float **data,						// input data
	float **query)						// input query
{
	n_pts_ = n;
	dim_ = d;
	m_ = m;
	U_ = U;
	appr_ratio_ = ratio;
	data_ = data;

	M_ = -1.0f;
	l2_alsh2_dim_ = d + 2 * m;	
	l2_alsh2_data_ = NULL;
	lsh_ = NULL;

	pre_processing(qn, query);
}

// -----------------------------------------------------------------------------
int L2_ALSH2::pre_processing(		// pre-processing of data
	int qn,								// number of queries
	float** query)						// input query
{
	// -------------------------------------------------------------------------
	//  calculate the norm of data and find the maximum norm of data and query
	// -------------------------------------------------------------------------
	float* norm = new float[n_pts_];
	for (int i = 0; i < n_pts_; i++) {
		norm[i] = 0.0f;
	}

	M_ = MINREAL;
	for (int i = 0; i < n_pts_; i++) {
		norm[i] = 0.0f;
		for (int j = 0; j < dim_; j++) {
			norm[i] += data_[i][j] * data_[i][j];
		}
		norm[i] = sqrt(norm[i]);

		if (norm[i] > M_) M_ = norm[i];
	}

	float tmp_norm = -1.0f;
	for (int i = 0; i < qn; i++) {
		tmp_norm = 0.0f;
		for (int j = 0; j < dim_; j++) {
			tmp_norm += query[i][j] * query[i][j];
		}
		tmp_norm = sqrt(tmp_norm);

		if (tmp_norm > M_) M_ = tmp_norm;
	}

	// -------------------------------------------------------------------------
	//  construct new data and indexing
	// -------------------------------------------------------------------------
	float scale = U_ / M_;
	int exponent = -1;

	printf("Construct L2_ALSH2 Data: ");
	l2_alsh2_data_ = new float*[n_pts_];
	for (int i = 0; i < n_pts_; i++) {
		l2_alsh2_data_[i] = new float[l2_alsh2_dim_];

		norm[i] = norm[i] * scale;
		for (int j = 0; j < l2_alsh2_dim_; j++) {
			if (j < dim_) {
				l2_alsh2_data_[i][j] = data_[i][j] * scale;
			}
			else if (j < dim_ + m_) {
				exponent = (int)pow(2.0f, j - dim_ + 1);
				l2_alsh2_data_[i][j] = pow(norm[i], exponent);
			}
			else {
				l2_alsh2_data_[i][j] = 0.5f;
			}
		}
	}
	printf("finish!\n\n");

	// -------------------------------------------------------------------------
	//  indexing the new data using qalsh
	// -------------------------------------------------------------------------
	if (indexing()) return 1;

	display_params();				// display parameters
	delete[] norm; norm = NULL;		// Release space

	return 0;
}

// -----------------------------------------------------------------------------
void L2_ALSH2::display_params()		// display parameters
{
	printf("Parameters of L2_ALSH2:\n");
	printf("    n = %d\n", n_pts_);
	printf("    d = %d\n", dim_);
	printf("    m = %d\n", m_);
	printf("    U = %.2f\n", U_);
	printf("    c = %.2f\n", appr_ratio_);
	printf("    M = %.2f\n\n", M_);
}

// -----------------------------------------------------------------------------
int L2_ALSH2::indexing()			// indexing the new data
{
	lsh_ = new QALSH_Col(l2_alsh2_data_, n_pts_, l2_alsh2_dim_, appr_ratio_);
	
	if (lsh_ != NULL) return 0;
	else return 1;
}

// -----------------------------------------------------------------------------
int L2_ALSH2::kmip(					// top-k approximate mip search
	float* query,						// input query
	int top_k,							// top-k value
	MaxK_List* list)					// top-k mip results
{
	int num_of_verf = 0;			// num of verification (NN and MIP calc)

	// -------------------------------------------------------------------------
	//  Construct L2_ALSH2 query
	// -------------------------------------------------------------------------
	float norm_q = 0.0f;			// calc norm of query
	for (int i = 0; i < dim_; i++) {
		norm_q += query[i] * query[i];
	}
	norm_q = sqrt(norm_q);

	int exponent = -1;
	float scale = U_ / M_;
	float* l2_alsh2_query = new float[l2_alsh2_dim_];

	norm_q = norm_q * scale;
	for (int i = 0; i < l2_alsh2_dim_; i++) {
		if (i < dim_) {
			l2_alsh2_query[i] = query[i] * scale;
		}
		else if (i < dim_ + m_) {
			l2_alsh2_query[i] = 0.5f;
		}
		else {
			exponent = (int)pow(2.0f, i - dim_ - m_ + 1);
			l2_alsh2_query[i] = pow(norm_q, exponent);
		}
	}

	// -------------------------------------------------------------------------
	//  Perform knn search via qalsh
	// -------------------------------------------------------------------------
	int l2_alsh2_top_k = top_k;
	MinK_List *nn_list = new MinK_List(l2_alsh2_top_k);

	num_of_verf += lsh_->knn(MAXREAL, l2_alsh2_query, l2_alsh2_top_k, nn_list);

	// -------------------------------------------------------------------------
	//  Compute inner product for candidates returned by qalsh
	// -------------------------------------------------------------------------
	for (int i = 0; i < l2_alsh2_top_k; i++) {
		int id = nn_list->ith_smallest_id(i);
		float ip = calc_inner_product(data_[id], query, dim_);

		list->insert(ip, id + 1);
	}
	num_of_verf += l2_alsh2_top_k;

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] l2_alsh2_query; l2_alsh2_query = NULL;
	delete nn_list; nn_list = NULL;

	return num_of_verf;
}
