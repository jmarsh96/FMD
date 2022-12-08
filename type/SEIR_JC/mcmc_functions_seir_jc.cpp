#include <iostream>
#include <fstream>
#include <Rcpp.h>

// These are for sleep using
// std::this_thread::sleep_for(std::chrono::milliseconds(x));
#include <chrono>
#include <thread>




using namespace Rcpp;


// Helper functions
template<typename T>
void WriteVector(std::vector<T> vec, std::ofstream &myfile)
{
  for(size_t i = 0; i < vec.size(); i++)
  {
    myfile << vec[i] << " ";
  }
  myfile << std::endl;
}

template<typename T>
void PrintVector(std::vector<T> vec)
{
  for(size_t i = 0; i < vec.size(); i++)
  {
    std::cout << vec[i] << " ";
  }
  std::cout << std::endl;
}

// Determine if two vectors are equal in the sense they contain the same elements
template<typename T>
bool CompareVec(std::vector<T>& v1, std::vector<T>& v2)
{
  std::sort(v1.begin(), v1.end());
  std::sort(v2.begin(), v2.end());
  return v1 == v2;
}

// Replicates the R sugar - the following are equivalent
// WhichVec(x,vec) == which(x==vec)
template<typename T>
std::vector<int> WhichVec(T x, const std::vector<T> vec)
{
  std::vector<int> out;
  for(size_t i = 0; i < vec.size(); i++)
  {
    if(vec[i] == x) out.push_back(i);
  }
  return out;
}

// Replicates the R sugar - the following are equivalent
// WhichNotVec(x,vec) == which(x!=vec)
template<typename T>
std::vector<int> WhichNotVec(T x, const std::vector<T> vec)
{
  std::vector<int> out;
  for(size_t i = 0; i < vec.size(); i++)
  {
    if(vec[i] != x) out.push_back(i);
  }
  return out;
}

// Replicates the r sugar - the following are equivalent
// subset(index,x) == x[index]
template<typename T>
std::vector<T> subset(std::vector<int> index, std::vector<T> x)
{
  std::vector<T> out;
  for(auto idx : index)
  {
    out.push_back(x[idx]);
  }
  return out;
}

std::vector<double> subset(std::vector<int> index, std::vector<double> x)
{
  std::vector<double> out;
  for(auto idx : index)
  {
    out.push_back(x[idx]);
  }
  return out;
}

// For two vectors a and b, this is equivalent to c(a,b)
template<typename T>
std::vector<T> CombineVectors(std::vector<T> a, std::vector<T> b)
{
  std::vector<T> ab;
  ab.reserve( a.size() + b.size() ); // preallocate memory
  ab.insert( ab.end(), a.begin(), a.end() );
  ab.insert( ab.end(), b.begin(), b.end() );
  return ab;
}

std::vector<int> sort_unique(std::vector<int> vec)
{
  std::unordered_set<int> s;
  for (int i : vec)
    s.insert(i);
  vec.assign( s.begin(), s.end() );
  sort( vec.begin(), vec.end() );
  return vec;
}

double min(double x, double y)
{
  if(x < y) return x;
  return y;
}

// Sort vector x by vector y
std::vector<int> sort_two(std::vector<int> x, std::vector<double> y)
{
  std::vector<int> indices(x.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(),
            [&](int A, int B) -> bool {
              return y[A] < y[B];
            });
  return subset(indices, x);
}

// Return the index which sorts by x, then y
std::vector<int> sort_three(std::vector<int> x, std::vector<double> y, std::vector<int> z)
{
  std::vector<int> indices(x.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&](int i, int j){
    if(x[i]==x[j]) {
      if(y[i]==y[j]) {
        return z[i] < z[j];
      }
      return y[i] < y[j];
    }
    return x[i] < x[j];
  });
  return indices;
}


// random shuffler
inline int randWrapper(const int n) { return floor(unif_rand()*n); }
std::vector<int> RandomShuffle(std::vector<int> x)
{
  std::random_shuffle(x.begin(),x.end(),randWrapper);
  return x;
}




// Data class 

class Data {
public:
  // Epi data
  int N; // total number of individuals
  int n_I; // number of infected individuals
  
  // Vectors containing infection and removal times
  std::vector<double> t_e;
  std::vector<double> t_i;
  std::vector<double> t_r;
  
  // Source of infection
  std::vector<int> source;
  
  // Spatial coordinates
  std::vector<double> x;
  std::vector<double> y;
  
  // Index of infected individuals
  std::vector<int> ever_infected;
  
  // Genetic data
  int num_seqs;
  int sequence_length;
  std::vector<int> genetic_ids;
  std::vector<double> sample_times;
  std::vector<int> subtype_numbers;
  std::vector<int> gen_source;
  IntegerMatrix gen_matrix;
  
  std::vector<int> coltime_distances;
  
  
  std::vector<int> imputed_nodes;
  std::vector<int> nodes_to_impute;
  std::vector<double> times_to_impute;
  std::vector<int> subtypes_to_impute;
  std::vector<int> imputation_type;
  
  double genetic_contribution;
  
  
  // Likelihood
  double loglik;
  double removal_likelihood;
  double transmission_likelihood;
  double genetic_likelihood;
  double exposure_likelihood;
  
  std::vector<double> parameters;
  
  // Constructor
  Data(int sequence_length_,
       NumericVector t_e_,
       NumericVector t_i_,
       NumericVector t_r_,
       IntegerVector source_,
       NumericVector x_,
       NumericVector y_,
       NumericVector parameters_,
       IntegerVector genetic_ids_,
       NumericVector sample_times_,
       IntegerVector subtype_numbers_,
       IntegerMatrix gen_matrix_);
  
  // Simulation constructor
  Data(NumericVector t_e,
       IntegerVector source,
       IntegerVector genetic_ids,
       NumericVector sample_times,
       IntegerVector subtype_numbers);
  
  
  // Functions
  std::vector<int> ReturnPossibleInfectors(int target);
  
  // Likelihood functions
  double CalculateDistance(int i, int j);
  double CalculateDistanceSum(bool verbose = false);
  double CalculateDoubleSum(bool verbose = false);
  void CalculateTransmissionLikelihood(bool verbose = false);
  void CalculateRemovalLikelihood(bool verbose = false);
  void CalculateGeneticLikelihood(bool verbose = false);
  void CalculateExposureLikelihood(bool verbose = false);
  void CalculateLoglik(bool verbose = false);
  
  void UpdateParameterMetropolisHastingsUniform(int idx, 
                                                double proposal_variance,
                                                double prior_lower,
                                                double prior_upper,
                                                int &nacc, 
                                                bool verbose = false);
  
  // Genetic imputation functions
  std::vector<int> ReturnIDsBetweenNodes(int node, 
                                         std::vector<int> child_nodes);
  void ReturnNodesToImpute(bool verbose = false);
  void UpdateImputedNodes();
  std::vector<int> DetermineOrderOfImputation();
  int FindGeneticLoc(int id, double time, int subtype);
  std::vector<int> ReturnObservedChildren(int node);
  void InitialiseImputedNodes();
  void ReturnNodesToImputeInitial();
  bool DoesNodeNeedToBeImputed(int target, int subtype);
  void ReturnNodesToImputeInitial2();
  bool DoesTargetHaveSequenceGreaterThanTime(int target, double time, int subtype);
  int ReturnExteriorTarget(int node);
  void ReturnExteriorNodes();
  std::vector<int> DetermineImputationGroups(std::vector<int> order);
  void CalculateContribRevUpper(int imp_idx, Data data_can, bool verbose = false);
  void CalculateContribRevExterior(int imp_idx, Data data_can, bool verbose = false);
  void CalculateContribRevInterior(int imp_idx, Data data_can, bool verbose = false);
  void ImputeUpperNode(int imp_idx, Data data_cur,
                       IntegerMatrix &gen_matrix_can, bool verbose = false);
  void ImputeExteriorNode(int imp_idx, Data data_cur,
                          IntegerMatrix &gen_matrix_can, bool verbose = false);
  void ImputeInteriorNode(int imp_idx, Data data_cur, 
                          IntegerMatrix &gen_matrix_can, bool verbose = false);
  void ImputeColtimeDistances(Data data_cur, bool verbose = false);
  std::vector<int> ReturnImpClusterOriginator();
  int ReturnMinDistanceTriangle(int imp_idx, const IntegerMatrix gen_matrix_can);
  
  // Gen source calculations
  void CalculateGenSourceVector();
  int ReturnSourceSequence(const int sequence_loc, bool verbose = false);
  std::vector<int> ReturnTargetSwabsAtTime(const int target, const double time, 
                                           const int subtype);
  std::vector<int> ReturnSequencesAtTime(int target, double time, 
                                         int subtype);
  std::vector<int> ReturnTargetSequences(const int target);
  
  
  // Augmented data functions
  void UpdateInfectionTimeSource(int &nacc, bool verbose = false);
  void UpdateInfectionTimeNoSource(int &nacc, bool verbose = false);
  void UpdateDistanceSingle(int node, int &nacc);
  void UpdateExposureTime(int &nacc, bool verbose = false);
  void UpdateInfectionTime(int &nacc, bool verbose = false);
  
  void WriteOutputToFile(std::ofstream &myfile);
  void WriteConfigurationToFile(std::ofstream &myfile);
  
  // Tree search functions
  int CalculateDistanceBetweenNodes(int node_i, int node_j);
  int ReturnLCA(int node_i, int node_j);
  std::vector<int> ReturnPathToRoot(int node);
  std::vector<int> ReturnSurroundingNodes(int node);
  
};

// Constructor for simulation
Data::Data(NumericVector t_e_,
           IntegerVector source_,
           IntegerVector genetic_ids_,
           NumericVector sample_times_,
           IntegerVector subtype_numbers_)
{
  for(int i = 0; i < t_e_.length(); i++)
  {
    t_e.push_back(t_e_[i]);
    int cur_source = source_[i];
    if(cur_source != -1) cur_source--;
    source.push_back(cur_source);
  }
  
  for(int i = 0; i < genetic_ids_.length(); i++)
  {
    int cur_gen_id = genetic_ids_[i] - 1;
    genetic_ids.push_back(cur_gen_id);
    sample_times.push_back(sample_times_[i]);
    subtype_numbers.push_back(subtype_numbers_[i]);
    imputed_nodes.push_back(0);
  }
}

// Constructor for the data class
Data::Data(int sequence_length_,
           NumericVector t_e_,
           NumericVector t_i_,
           NumericVector t_r_,
           IntegerVector source_,
           NumericVector x_,
           NumericVector y_,
           NumericVector parameters_,
           IntegerVector genetic_ids_,
           NumericVector sample_times_,
           IntegerVector subtype_numbers_,
           IntegerMatrix gen_matrix_)
{
  // Initialise epi data
  N = t_e_.length();
  int total_infected = 0;
  
  for(int i = 0; i < N; i++)
  {
    t_e.push_back(t_e_[i]);
    t_i.push_back(t_i_[i]);
    t_r.push_back(t_r_[i]);
    x.push_back(x_[i]);
    y.push_back(y_[i]);
    
    int cur_source = source_[i];
    if(cur_source == -1)
    {
      source.push_back(-1);
    }
    else
    {
      source.push_back(cur_source - 1);
    }
    
    if(t_e_[i] != -1)
    {
      total_infected++;
      ever_infected.push_back(i);
    }
  }
  n_I = total_infected;
  
  // Genetic data
  num_seqs = genetic_ids_.length();
  sequence_length = sequence_length_;
  for(int i = 0; i < num_seqs; i++)
  {
    genetic_ids.push_back(genetic_ids_[i] - 1);
    sample_times.push_back(sample_times_[i]);
    subtype_numbers.push_back(subtype_numbers_[i]);
    imputed_nodes.push_back(0);
  }
  
  // Fill in matrices
  gen_matrix = gen_matrix_;
  
  // parameters
  for(auto i : parameters_)
  {
    parameters.push_back(i);
  }
  
  
  CalculateGenSourceVector();
  std::vector<int> imported_idx = WhichVec(-1, gen_source);
  std::vector<int> coltime_distances_(imported_idx.size());
  
  bool impute_genetic_distances = true;
  if(impute_genetic_distances)
  {
    bool valid_config = false;
    int counter = 0;
    while(!valid_config)
    {
      Data data_can((*this));
      data_can.coltime_distances = coltime_distances_;
      data_can.UpdateImputedNodes();
      data_can.CalculateLoglik();
      if(data_can.loglik != R_NegInf)
      {
        (*this) = data_can;
        valid_config = true;
        break;
      }
      if(counter > 5000) stop("cannot find valid configuration");
      counter++;
    }
  }
  else
  {
    CalculateLoglik();
  }
  
  
  
  
  std::cout << "Data initialised - loglik = " << loglik << std::endl;
}

double Data::CalculateDistance(int i, int j)
{
  double x_diff = x[i] - x[j];
  double y_diff = y[i] - y[j];
  double dist = sqrt(pow(x_diff,2)+pow(y_diff,2));
  return dist;
}

double Data::CalculateDistanceSum(bool verbose)
{
  double distance_sum = 0.0;
  for(int i = 0; i < N; i++)
  {
    double infection_time = t_e[i];
    if(infection_time != -1)
    {
      int cur_source = source[i];
      if(cur_source != -1)
      {
        double cur_dist = CalculateDistance(cur_source, i);
        distance_sum += cur_dist;
      }
    }
  }
  return distance_sum;
}

double Data::CalculateDoubleSum(bool verbose)
{
  verbose = false;
  double double_sum = 0.0;
  double kappa = parameters[3];
  if(verbose)
  {
    std::cout << "Calculating double sum" << std::endl;
  }
  for(auto i : ever_infected)
  {
    for(int j = 0; j < N; j++)
    {
      double cur_dist = CalculateDistance(i,j);
      double h_ij = exp(-1*kappa*cur_dist);
      double contribution = 0.0;
      if(t_e[j] == -1)
      {
        // j is never infected
        contribution = h_ij*(t_r[i]-t_i[i]);
      }
      else
      {
        contribution = h_ij*(std::min(t_r[i],t_e[j]) - std::min(t_i[i],t_e[j]));
      }
      double_sum += contribution;
      if(verbose)
      {
        std::cout << "i = " << i
                  << ", j = " << j
                  << ", contribution = " << contribution
                  << std::endl;
      }
    }
  }
  if(verbose)
  {
    std::cout << "double_sum = " << double_sum << std::endl;
  }
  return double_sum;
}

double CalculateMutationProbability(double rate, double t)
{
  double prob = 0.75*(1-exp(-4*rate*t));
  return prob;
}


void Data::CalculateTransmissionLikelihood(bool verbose)
{
  verbose = false;
  double transmission_likelihood_ = 0.0;
  double beta = parameters[1];
  double kappa = parameters[3];
  
  double double_sum = CalculateDoubleSum();
  double distance_sum = CalculateDistanceSum();
  
  transmission_likelihood_ = (n_I-1)*log(beta) - kappa*distance_sum -
    beta*double_sum;
  if(verbose)
  {
    std::cout << "Calculating transmission likelihood"
              << ", n_I = " << n_I
              << ", beta = " << beta
              << ", kappa = " << kappa
              << ", distance_sum = " << distance_sum
              << ", double_sum = " << double_sum
              << ", likelihood = " << transmission_likelihood_
              << std::endl;
  }
  transmission_likelihood = transmission_likelihood_;
}

bool DetermineIfImputedConfigSame(int imp_idx_cur, int imp_idx_can,
                                  Data data_cur, Data data_can)
{
  bool verbose = true;
  
  if(imp_idx_cur == -1 || imp_idx_can == -1) return false;
  int genetic_id_cur = data_cur.genetic_ids[imp_idx_cur];
  double sample_time_cur = data_cur.sample_times[imp_idx_cur];
  int subtype_cur = data_cur.subtype_numbers[imp_idx_cur];
  int type_cur = data_cur.imputed_nodes[imp_idx_cur];
  
  int genetic_id_can = data_can.genetic_ids[imp_idx_can];
  double sample_time_can = data_can.sample_times[imp_idx_can];
  int subtype_can = data_can.subtype_numbers[imp_idx_can];
  int type_can = data_can.imputed_nodes[imp_idx_can];
  
  if(verbose)
  {
    std::cout << " Determine if nodes " << imp_idx_cur << " and "
              << imp_idx_can << " are equal with IDs "
              << genetic_id_cur << " and " << genetic_id_can
              << " and times " << sample_time_cur
              << " and " << sample_time_can
              << std::endl;
    
  }
  
  if(genetic_id_cur != genetic_id_can) return false;
  if(sample_time_cur != sample_time_can) return false;
  if(subtype_cur != subtype_can) return false;
  if(type_cur != type_can) return false;
  
  std::vector<int> child_nodes_cur = WhichVec(imp_idx_cur, data_cur.gen_source);
  std::vector<int> child_nodes_can = WhichVec(imp_idx_can, data_can.gen_source);
  std::vector<int> obs_child_nodes_cur = 
    data_cur.ReturnObservedChildren(imp_idx_cur);
  std::vector<int> obs_child_nodes_can = 
    data_cur.ReturnObservedChildren(imp_idx_can);
  std::vector<int> direct_obs_child_nodes_cur;
  std::vector<int> direct_obs_child_nodes_can;
  int parent_cur = data_cur.gen_source[imp_idx_cur];
  int parent_can = data_can.gen_source[imp_idx_can];
  
  // Now check configuration specifics
  switch(type_cur)
  {
  case 1:
  {
    // Interior node, parent and obs children must be equal
    int parent_ed_cur = data_cur.genetic_ids[parent_cur];
    double parent_time_cur = data_cur.sample_times[parent_cur];
    int parent_subtype_cur = data_cur.subtype_numbers[parent_cur];
    
    int parent_ed_can = data_can.genetic_ids[parent_can];
    double parent_time_can = data_can.sample_times[parent_can];
    int parent_subtype_can = data_can.subtype_numbers[parent_can];
    
    if((parent_ed_cur != parent_ed_can) ||
       (parent_time_cur != parent_time_can) ||
       (parent_subtype_cur != parent_subtype_can))
    {
      return false;
    }
    
    return CompareVec(obs_child_nodes_cur,obs_child_nodes_can);
  }
  case 2:
    // Exterior node, observed and direct observed child nodes must be equal
    if(CompareVec(obs_child_nodes_cur,obs_child_nodes_can))
    {
      // Now look at direct observed children
      for(auto ch_node : obs_child_nodes_cur)
      {
        if(data_cur.gen_source[ch_node] == imp_idx_cur)
        {
          direct_obs_child_nodes_cur.push_back(ch_node);
        }
      }
      
      for(auto ch_node : obs_child_nodes_can)
      {
        if(data_can.gen_source[ch_node] == imp_idx_can)
        {
          direct_obs_child_nodes_can.push_back(ch_node);
        }
      }
      return CompareVec(direct_obs_child_nodes_cur,direct_obs_child_nodes_can);
    }
  case 3:
    // Upper node, only children must be equal, but we must also check
    // the ID and time since the children will be imputed
    for(auto ch_node_cur : child_nodes_cur)
    {
      int gen_id_cur = data_cur.genetic_ids[ch_node_cur];
      double samp_time_cur = data_cur.sample_times[ch_node_cur];
      int subt_cur = data_cur.subtype_numbers[ch_node_cur];
      
      bool node_found = false;
      for(auto ch_node_can : child_nodes_can)
      {
        int gen_id_can = data_can.genetic_ids[ch_node_can];
        double samp_time_can = data_can.sample_times[ch_node_can];
        int subt_can = data_can.subtype_numbers[ch_node_can];
        if((gen_id_cur == gen_id_can) &&
           (samp_time_cur == samp_time_can) &&
           (subt_cur == subt_can))
        {
          node_found = true;
          break;
        }
      }
      
      if(!node_found) return false;
    }
    return true;
  }
  stop("Code should not get this far in DetermineIfImputedConfigSame");
  return true;
}

void TriangleInequalityCheck(Data data)
{
  for(size_t j = 1; j < data.genetic_ids.size(); j++)
  {
    for(size_t i = 0; i < j; i++)
    {
      int node_i = data.genetic_ids[i];
      int node_j = data.genetic_ids[j];
      //int dist_sum = data.CalculateDistanceBetweenNodes(node_i, node_j);
      int dist_sum = 1;
      int dist = data.gen_matrix(node_i,node_j);
      if(dist_sum > dist)
      {
        std::cout << "Triangle inequality violated" << std::endl;
      }
    }
  }
}

void Data::CalculateRemovalLikelihood(bool verbose)
{
  double delta = parameters[4];
  double gamma = parameters[2];
  double removal_likelihood_ = 0.0;
  for(auto i : ever_infected)
  {
    removal_likelihood_ += R::dgamma(t_r[i]-t_i[i],delta,1/gamma,1);
  }
  removal_likelihood = removal_likelihood_;
}

void Data::CalculateExposureLikelihood(bool verbose)
{
  double zeta = parameters[5];
  double alpha = parameters[0];
  double exposure_likelihood_ = 0.0;
  if(verbose)
  {
    std::cout << "Calculaing the likelihood for the latent period with zeta = "
              << zeta << " and alpha = " << alpha << std::endl;
  }
  for(auto i : ever_infected)
  {
    double contrib = R::dgamma(t_i[i]-t_e[i],zeta,1/alpha,1);
    exposure_likelihood_ += contrib;
    if(verbose)
    {
      std::cout << " i = " << i << ", t_i[i]-t_e[i] = "
                << t_i[i]-t_e[i] << ", contribution = "
                << contrib << std::endl;
    }
  }
  exposure_likelihood = exposure_likelihood_;
}

void Data::CalculateGeneticLikelihood(bool verbose)
{
  //verbose = false;
  double genetic_likelihood_ = 0.0;
  double lambda = parameters[6]; // mutation rate
  
  double discrep_contrib = 0.0;
  double direct_contrib = 0.0;
  double coltime_contrib = 0.0;
  
  if(verbose)
  {
    std::cout << "Calculate contribution for genetic log likelihood with genetic source vector = ";
    PrintVector(gen_source);
  }
  for(size_t i = 0; i < gen_source.size(); i++)
  {
    int cur_source = gen_source[i];
    if(cur_source != -1)
    {
      
      int num_D = gen_matrix(i, cur_source);
      
      double time_diff = sample_times[i] - sample_times[cur_source];
      double mutation_prob = CalculateMutationProbability(lambda, time_diff);
      double contribution = (sequence_length-num_D)*log(1-mutation_prob) + 
        num_D*log(mutation_prob/3);
      genetic_likelihood_ += contribution;
      direct_contrib += contribution;
      if(verbose)
      {
        std::cout << "Calculate likelihood for sequence " << i
                  << " with source sequence " << cur_source << " with dist = "
                  << num_D << ", time_diff = " << time_diff
                  << ", mutation_prob = " << mutation_prob
                  << " and contribution = " << contribution << std::endl;
      }
      if(num_D < 0) stop("Negative distances in the matrix");
      
    }
    
    bool calculate_observation = true;
    if(calculate_observation && imputed_nodes[i] == 0)
    {
      int person = genetic_ids[i];
      genetic_likelihood_ -= log(t_r[person]-t_i[person]);
    }
  }
  
  // Calculate the contribution of the discrepancy from the imputed nodes
  if(verbose) std::cout << "Calculate discrepancy" << std::endl;
  bool calculate_discrepancy = true;
  if(calculate_discrepancy && false)
  {
    std::vector<int> imputed_idx = WhichNotVec(0, imputed_nodes);
    
    for(auto node : imputed_idx)
    {
      std::vector<int> surrounding_nodes = ReturnSurroundingNodes(node);
      if(verbose)
      {
        std::cout << " current idx = " << node << ", surrounding nodes = ";
        PrintVector(surrounding_nodes);
      }
      for(size_t j = 1; j < surrounding_nodes.size(); j++)
      {
        for(size_t i = 0; i < j; i++)
        {
          int node_i = surrounding_nodes[i];
          int node_j = surrounding_nodes[j];
          if(imputed_nodes[node_i] == 0 && imputed_nodes[node_j] == 0)
          {
            int dist = gen_matrix(node_i, node_j);
            int dist_sum = CalculateDistanceBetweenNodes(node_i, node_j);
            int discrepancy = dist_sum - dist;
            if(discrepancy < 0)
            {
              //genetic_likelihood_ -= 1e16;
              discrepancy = discrepancy*-1;
            }
            if(verbose)
            {
              std::cout << "For nodes " << node_i << " and " << node_j
                        << ", dist = " << dist << ", dist_sum = "
                        << dist_sum << ", discrepancy = "
                        << discrepancy;
            }
            
            double contribution = R::dgeom(discrepancy, 0.1, 1);
            genetic_likelihood_ += contribution;
            if(verbose)
            {
              std::cout << ", contribution = " << contribution << std::endl;
            }
          }
        }
      }
    }
  }
  
  if(calculate_discrepancy)
  {
    std::vector<int> imp_clust_originators = ReturnImpClusterOriginator();
    
    if(verbose)
    {
      std::cout << "Calculate discrepancy with the following cluster originators = ";
      PrintVector(imp_clust_originators);
    }
    
    for(auto imp_node : imp_clust_originators)
    {
      std::vector<int> obs_nodes = ReturnObservedChildren(imp_node);
      int parent = gen_source[imp_node];
      if(parent != -1)
      {
        if(imputed_nodes[parent] == 0) obs_nodes.push_back(parent);
      }
      
      if(verbose)
      {
        std::cout << " Current imp node = " << imp_node
                  << ", obs_nodes = ";
        PrintVector(obs_nodes);
      }
      
      for(size_t j = 1; j < obs_nodes.size(); j++)
      {
        for(size_t i = 0; i < j; i++)
        {
          int node_i = obs_nodes[i];
          int node_j = obs_nodes[j];
          
          
          int dist = gen_matrix(node_i, node_j);
          int dist_sum = CalculateDistanceBetweenNodes(node_i, node_j);
          int discrepancy = dist_sum - dist;
          if(discrepancy < 0)
          {
            genetic_likelihood_ -= 1e16;
            //discrepancy = discrepancy*-1;
          }
          if(verbose)
          {
            std::cout << "For nodes " << node_i << " and " << node_j
                      << ", dist = " << dist << ", dist_sum = "
                      << dist_sum << ", discrepancy = "
                      << discrepancy;
          }
          
          double contribution = R::dgeom(discrepancy, 0.9, 1);
          genetic_likelihood_ += contribution;
          discrep_contrib += contribution;
          if(verbose)
          {
            std::cout << ", contribution = " << contribution << std::endl;
          }
        }
      }
    }
  }
  
  
  
  
  // Calculate contribution of the colonisation time distances
  if(verbose)
  {
    std::cout << "Calculate contribution of coltime distances" << std::endl;
  }
  std::vector<int> import_idx = WhichVec(-1, gen_source);
  for(size_t i = 0; i < import_idx.size(); i++)
  {
    int current_idx = import_idx[i];
    int current_distance = coltime_distances[i];
    double current_time = sample_times[current_idx];
    int current_person = genetic_ids[current_idx];
    double col_time = t_e[current_person];
    double time_diff = current_time - col_time;
    double mutation_prob = CalculateMutationProbability(lambda, time_diff);
    double contribution = R::dbinom(current_distance, sequence_length, 
                                    mutation_prob, 1);
    genetic_likelihood_ += contribution;
    coltime_contrib += contribution;
    if(verbose)
    {
      std::cout << "Current coltime distance = " << current_distance
                << ", time diff = " << time_diff
                << ", binomial contributon = " << R::dbinom(current_distance, sequence_length, mutation_prob, 1) 
                << std::endl;
    }
  }
  

  if(verbose)
  {
    std::cout << "Direct contribution = " << direct_contrib
              << ", discrepancy contribution = " << discrep_contrib
              << ", coltime contribution = " << coltime_contrib
              << ", genetic ll = " << genetic_likelihood_
              << std::endl;
  }
  
  genetic_likelihood = genetic_likelihood_;
}

void Data::CalculateLoglik(bool verbose)
{
  CalculateTransmissionLikelihood(verbose);
  CalculateRemovalLikelihood(verbose);
  CalculateExposureLikelihood(verbose);
  //verbose = true;
  CalculateGeneticLikelihood(verbose);
  loglik = transmission_likelihood + removal_likelihood + genetic_likelihood +
    exposure_likelihood;
  if(verbose)
  {
    std::cout << "Calculated log likelihood"
              << ", transmission loglik = " << transmission_likelihood
              << ", removal loglik = " << removal_likelihood
              << ", exposure loglik = " << exposure_likelihood
              << ", genetic loglik = " << genetic_likelihood
              << std::endl;
  }
}

void Data::WriteConfigurationToFile(std::ofstream &myfile)
{
  myfile << loglik << std::endl;
  
  for(size_t i = 0; i < t_e.size(); i++)
  {
    myfile << t_e[i] << " ";
  }
  myfile << std::endl;
  
  for(size_t i = 0; i < t_e.size(); i++)
  {
    myfile << t_i[i] << " ";
  }
  myfile << std::endl;
  
  for(size_t i = 0; i < source.size(); i++)
  {
    myfile << source[i] << " ";
  }
  myfile << std::endl;
  
  for(size_t i = 0; i < genetic_ids.size(); i++)
  {
    myfile << genetic_ids[i] << " ";
  }
  myfile << std::endl;
  
  for(size_t i = 0; i < sample_times.size(); i++)
  {
    myfile << sample_times[i] << " ";
  }
  myfile << std::endl;
  
  for(size_t i = 0; i < imputed_nodes.size(); i++)
  {
    myfile << imputed_nodes[i] << " ";
  }
  myfile << std::endl;
  
  for(size_t i = 0; i < gen_source.size(); i++)
  {
    myfile << gen_source[i] << " ";
  }
  myfile << std::endl;
  
  for(int i = 0; i < gen_matrix.nrow(); i++)
  {
    for(int j = 0; j < gen_matrix.ncol(); j++)
    {
      myfile << gen_matrix(i,j) << " ";
    }
    myfile << std::endl;
  }
  myfile << std::endl;
}

void Data::WriteOutputToFile(std::ofstream &myfile)
{
  for(size_t i = 0; i < parameters.size(); i++)
  {
    myfile << parameters[i] << " ";
  }
  
  myfile << loglik << " ";
  
  for(size_t i = 0; i < t_e.size(); i++)
  {
    myfile << t_e[i] << " ";
  }
  
  for(size_t i = 0; i < t_e.size(); i++)
  {
    myfile << t_i[i] << " ";
  }
  
  for(size_t i = 0; i < source.size(); i++)
  {
    myfile << source[i] << " ";
  }
  
  myfile << std::endl;
}

int SampleVector(std::vector<int> x)
{
  double U = R::runif(0.0,1.0);
  int length = x.size();
  return x[floor(U*length)];
}


std::vector<int> Data::ReturnTargetSwabsAtTime(const int target, 
                                               const double time, 
                                               const int subtype)
{
  std::vector<int> out;
  for(size_t i = 0; i < genetic_ids.size(); i++)
  {
    if(genetic_ids[i] == target && sample_times[i] == time && 
       subtype_numbers[i] == subtype)
    {
      out.push_back(i);
    }
  }
  return out;
}

std::vector<int> Data::ReturnSequencesAtTime(int target, double time, 
                                             int subtype)
{
  
  // Return all swabs from the target at that time
  std::vector<int> sequences_at_time = ReturnTargetSwabsAtTime(target, 
                                                               time, 
                                                               subtype);
  
  // Now check the swabs of either infections or the infector (or both if an importation)
  if(time == t_e[target])
  {
    // The time is at the time of colonisation for the target, therefore check in the source of colonisation
    int target_source = source[target];
    if(target_source != -1)
    {
      // Only check if the source if the current target is NOT an importation
      // First check source swabs
      //IntegerVector source_swabs = ReturnTargetSwabsAtTime(data, target_source, subtype, time);
      std::vector<int> source_swabs = ReturnTargetSwabsAtTime(target_source, 
                                                              time, 
                                                              subtype);
      if(source_swabs.size() > 0)
      {
        sequences_at_time = CombineVectors(sequences_at_time, source_swabs);
      }
      
      // Now check source infections, only those who are colonised at that time
      std::vector<int> source_colonisations = WhichVec(target_source, source);
      std::vector<double> source_colonisations_t_e = 
        subset(source_colonisations, t_e);
      //IntegerVector source_colonisations_t_c = as<IntegerVector>(t_c[source_colonisations]);
      for(size_t i = 0; i < source_colonisations.size(); i++)
      {
        int current_source_colonisation = source_colonisations[i];
        double current_source_colonisation_time = source_colonisations_t_e[i];
        if(current_source_colonisation_time == time)
        {
          std::vector<int> current_source_colonisation_swabs = 
            ReturnTargetSwabsAtTime(current_source_colonisation, time, subtype);
          
          if(current_source_colonisation_swabs.size() > 0)
          {
            sequences_at_time = CombineVectors(sequences_at_time, 
                                               current_source_colonisation_swabs);
          }
        }
      }
    }
  }
  
  // Now check colonisations from the target at the time
  std::vector<int> target_colonisations = WhichVec(target, source);
  std::vector<double> target_colonisations_t_e = 
    subset(target_colonisations, t_e);
  for(size_t i = 0; i < target_colonisations.size(); i++)
  {
    int current_target_colonisation = target_colonisations[i];
    double current_target_colonisation_time = target_colonisations_t_e[i];
    if(current_target_colonisation_time == time)
    {
      std::vector<int> current_target_colonisation_swabs = 
        ReturnTargetSwabsAtTime(current_target_colonisation, time, subtype);
      if(current_target_colonisation_swabs.size() > 0)
      {
        sequences_at_time = CombineVectors(sequences_at_time, 
                                           current_target_colonisation_swabs);
      }
    }
  }
  
  // Now filter only unique entries and order them
  return sort_unique(sequences_at_time);
}

// Return the sequences for a specific individual, and also return the sequences
// that are from other infections at the time of exposure
std::vector<int> Data::ReturnTargetSequences(const int target)
{
  bool verbose = false;
  std::vector<int> target_sequences = WhichVec(target, genetic_ids);
  std::vector<int> target_infections = WhichVec(target, source);
  if(verbose)
  {
    std::cout << "Return target sequences for ID = " << target << std::endl;
    std::cout << " Target sequences = ";
    PrintVector(target_sequences);
    std::cout << " Target infections = ";
    PrintVector(target_infections);
  }
  for(auto cur_infection : target_infections)
  {
    std::vector<int> cur_inf_gen_idx = WhichVec(cur_infection, genetic_ids);
    if(verbose)
    {
      std::cout << "  Current infection = " << cur_infection
                << " at infection time " << t_e[cur_infection] << std::endl;
    }
    for(auto idx : cur_inf_gen_idx)
    {
      double time = sample_times[idx];
      if(verbose)
      {
        std::cout << "   Current index = " << idx << " at time = " << time;
      }
      if(time == t_e[cur_infection])
      {
        target_sequences.push_back(idx);
        if(verbose)
        {
          std::cout << ", add to target sequences" << std::endl;
        }
      }
      else
      {
        if(verbose)
        {
          std::cout << ", not at time of exposure" << std::endl;
        }
      }
    }
  }
  if(verbose)
  {
    std::cout << " Target sequences = ";
    PrintVector(target_sequences);
  }
  return target_sequences;
}

int Data::ReturnSourceSequence(const int sequence_loc, bool verbose)
{
  verbose = false;
  int current_target = genetic_ids[sequence_loc];
  double current_time = sample_times[sequence_loc];
  int current_subtype = subtype_numbers[sequence_loc];
  if(verbose)
  {
    std::cout << "Return the source for sequence " << sequence_loc
              << " from individual " << current_target
              << " at time " << current_time
              << " with subtype " << current_subtype << std::endl;
  }
  
  // Check if there are multiple sequences sampled at the same time
  std::vector<int> sequences_at_time = ReturnSequencesAtTime(current_target, 
                                                             current_time, 
                                                             current_subtype);
  
  if(sequences_at_time.size() > 1)
  {
    
    std::vector<int> sequence_locs = WhichVec(sequence_loc, sequences_at_time);
    int sequence_loc_idx = sequence_locs[0];
    if(verbose)
    {
      std::cout << "There are multiple sequences at the same time = ";
      PrintVector(sequences_at_time);
      std::cout << "Sequence loc = " << sequence_loc_idx
                << ", return sequence " 
                << sequences_at_time[sequence_loc_idx - 1] << std::endl;
    }
    
    if(sequence_loc_idx > 0)
    {
      return sequences_at_time[sequence_loc_idx - 1];
    }
  }
  
  // Look backwards to find the source sequence
  while(current_target != -1)
  {
    // Return all sequences from the target
    std::vector<int> target_sequences = ReturnTargetSequences(current_target);
    
    if(verbose)
    {
      std::cout << " Current target = " << current_target
                << ", target sequences = ";
      PrintVector(target_sequences);
    }
    
    // Return all sequences that are at the time of exposure of infections
    //std::vector<int> target_enf_seqs;
    
    
    // Find the biggest sequence such that the sampling time is less than the
    // current time
    int source_sequence = -1;
    double source_sequence_time = -1;
    
    for(auto cur_seq : target_sequences)
    {
      //std::cout << " cur_seq = " << cur_seq << std::endl;
      if(cur_seq < 0) stop("Invalid sequence index");
      double cur_time = sample_times[cur_seq];
      int cur_subtype = subtype_numbers[cur_seq];
      if((cur_time <= current_time) &&
         (cur_time > source_sequence_time) &&
         (cur_subtype == current_subtype) &&
         (cur_seq != sequence_loc))
      {
        if(verbose)
        {
          std::cout << " Found a previous source sequence at loc "
                    << cur_seq << " at time " << cur_time
                    << std::endl;
        }
        source_sequence = cur_seq;
        source_sequence_time = cur_time;
      }
    }
    //std::cout << " source sequence = " << source_sequence << std::endl;
    
    if(source_sequence != -1)
    {
      // A source sequence has been found, check for duplicates and then return
      
      if(verbose)
      {
        std::cout << " Source sequence found - ID = " << source_sequence
                  << " at time " << source_sequence_time
                  << " with subtype " << current_subtype << std::endl;
      }
      
      sequences_at_time = ReturnSequencesAtTime(genetic_ids[source_sequence], 
                                                source_sequence_time, 
                                                current_subtype);
      if(sequences_at_time.size() > 1)
      {
        std::vector<int> sequence_locs = WhichVec(source_sequence, 
                                                  sequences_at_time);
        int sequence_loc_idx = sequence_locs[0];
        if(sequence_loc_idx > 0)
        {
          return sequences_at_time[sequence_loc_idx - 1];
        }
      }
      else
      {
        return source_sequence;
      }
    }
    else
    {
      // No source sequence has been found, change source and try again
      //current_time = t_e[current_target] + 0.000001;
      current_time = t_e[current_target];
      current_target = source[current_target];
      
      //std::cout << " Look at the source, current target = " << current_target 
      //          << std::endl;
      //std::cout << " current time = " << current_time << std::endl;
      
    }
  }
  // If the function gets this far, there is no source sequence
  return -1;
}

// [[Rcpp::export]]
IntegerVector ReturnGenSourceVector_R(List data)
{
  //Rcout << "asd" << std::endl << std::flush;
  IntegerVector genetic_ids = data["genetic_ids"];
  //return genetic_ids;
  NumericVector sample_times = data["sample_times"];
  IntegerVector subtype_numbers = data["subtype_numbers"];
  IntegerVector source = data["source"];
  NumericVector branching_times = data["branching_times"];
  
  Data temp_data(branching_times, source, genetic_ids, sample_times, subtype_numbers);
  
  //std::cout << "Data initialised" << std::endl;
  temp_data.CalculateGenSourceVector();
  IntegerVector out(temp_data.gen_source.size());
  for(size_t i = 0; i < temp_data.gen_source.size(); i++)
  {
    out[i] = temp_data.gen_source[i];
  }
  return out;
}

void Data::CalculateGenSourceVector()
{
  bool verbose = false;
  //std::cout << "In function call Data::CalculateGenSourceVector()" << std::endl;
  num_seqs = genetic_ids.size();
  std::vector<int> gen_source_out(num_seqs);
  
  if(verbose)
  {
    std::cout << "num_seqs = " << num_seqs << std::endl;
    std::cout << "genetic ids = ";
    PrintVector(genetic_ids);
    std::cout << "sample times = ";
    PrintVector(sample_times);
  }
  //std::this_thread::sleep_for(std::chrono::milliseconds(5000));
  for(int i = 0; i < num_seqs; i++)
  {
    //std::cout << "Returning source for sequence i = " << i << std::endl;
    gen_source_out[i] = ReturnSourceSequence(i);
  }
  gen_source = gen_source_out;
}


// Determines the type of imputation
// 1 - Interior node
// 2 - Exterior node
// 3 - Upper
// Main idea: first do a breadth first search and identify all nodes
// that consist in an "upper"group, i.e. not enough observed information
std::vector<int> Data::DetermineImputationGroups(std::vector<int> order)
{
  std::vector<int> idx_to_update;
  // Loop through and check if the nodes are exterior
  
  for(auto node : order)
  {
    int parent = gen_source[node];
    if(parent == -1)
    {
      // check if there are any observed children
      std::vector<int> child_nodes = WhichVec(node, gen_source);
      bool observed_child_found = false;
      for(auto ch_node : child_nodes)
      {
        if(imputed_nodes[ch_node] == 0)
        {
          observed_child_found = true;
        }
      }
      if(!observed_child_found)
      {
        // No observed child has been found, we need to set all nodes
        // to be of imputation group 1
        
        // start at the parent node and breadth first search until all have
        // been found
        std::vector<int> Q = {node}; // should probably use std queue
        idx_to_update.push_back(node);
        while(Q.size() != 0)
        {
          std::vector<int> Q_new;
          for(auto current_node : Q)
          {
            std::vector<int> child_nodes = WhichVec(current_node, gen_source);
            for(auto current_child : child_nodes)
            {
              if(imputed_nodes[current_child] == 1)
              {
                idx_to_update.push_back(current_child);
                Q_new.push_back(current_child);
                //current_nodes_to_impute.push_back(current_child);
              }
            }
          }
          Q = Q_new;
        }
      }
    }
  }
  // create a vector of the groups
  std::vector<int> imputation_groups;
  for(auto imp_node : order)
  {
    std::vector<int> loc = WhichVec(imp_node, idx_to_update);
    int cur_grp = (loc.size() == 0) ? 1 : 3;
    if(cur_grp == 1)
    {
      int parent = gen_source[imp_node];
      cur_grp = (parent == -1) ? 2 : 1;
    }
    imputation_groups.push_back(cur_grp);
  }
  return imputation_groups;
}

// returns a vector of idx of nodes that start imputation cluster
std::vector<int> Data::ReturnImpClusterOriginator()
{
  std::vector<int> imputed_originators;
  
  std::vector<int> imputed_idx = WhichNotVec(0, imputed_nodes);
  for(auto idx : imputed_idx)
  {
    int parent = gen_source[idx];
    if(parent == -1)
    {
      imputed_originators.push_back(idx);
    }
    else
    {
      if(imputed_nodes[parent] == 0)
      {
        imputed_originators.push_back(idx);
      }
    }
  }
  return imputed_originators;
}


std::vector<int> Data::DetermineOrderOfImputation()
{
  // Main idea: start at each imputed node (gs=-1) and breadth first search 
  // until all imputed nodes have been found
  std::vector<int> nodes_to_impute_;
  std::vector<int> imported_idx = WhichVec(-1, gen_source);
  for(size_t i = 0; i < imported_idx.size(); i++)
  {
    //std::vector<int> current_nodes_to_impute;
    // Q is the queue (see wiki for BFS)
    int root = imported_idx[i];
    if(imputed_nodes[root] == 1)
    {
      nodes_to_impute_.push_back(root);
      //current_nodes_to_impute.push_back(root);
    }
    std::vector<int> Q = {root};
    while(Q.size() != 0)
    {
      std::vector<int> Q_new;
      for(size_t j = 0; j < Q.size(); j++)
      {
        int current_node = Q[j];
        std::vector<int> child_nodes = WhichVec(current_node, gen_source);
        for(size_t k = 0; k < child_nodes.size(); k++)
        {
          int current_child = child_nodes[k];
          if(imputed_nodes[current_child] == 1)
          {
            nodes_to_impute_.push_back(current_child);
            //current_nodes_to_impute.push_back(current_child);
          }
          Q_new.push_back(current_child);
        }
      }
      Q = Q_new;
    }
  }
  return nodes_to_impute_;
}

int Data::FindGeneticLoc(int id, double time, int subtype)
{
  for(size_t i = 0; i < genetic_ids.size(); i++)
  {
    if(id == genetic_ids[i] &&
       time == sample_times[i] &&
       subtype == subtype_numbers[i])
    {
      return i;
    }
  }
  return -1;
}

bool DetermineIfNodesAreEqual(std::vector<int> nodes_cur, std::vector<int> nodes_can,
                              Data data_cur, Data data_can)
{
  // First check if lengths are the same
  if(nodes_cur.size() == nodes_can.size())
  {
    // Look through all the current nodes and check that the same sequence appears
    // in the candidate nodes
    size_t counter = 0;
    for(size_t i = 0; i < nodes_cur.size(); i++)
    {
      int cur_node = nodes_cur[i];
      int cur_id = data_cur.genetic_ids[cur_node];
      double cur_time = data_cur.sample_times[cur_node];
      int cur_subtype = data_cur.subtype_numbers[cur_node];
      
      for(size_t j = 0; j < nodes_can.size(); j++)
      {
        int can_node = nodes_can[j];
        int can_id = data_can.genetic_ids[can_node];
        double can_time = data_can.sample_times[can_node];
        int can_subtype = data_can.subtype_numbers[can_node];
        
        if(cur_id == can_id && cur_time == can_time && cur_subtype == can_subtype)
        {
          // match has been found
          counter++;
          break;
        }
      }
    }
    if(counter == nodes_cur.size()) return true;
  }
  return false;
}

bool DetermineIfNodeAndObsChildNodesAreSame(int node_cur, int node_can, Data data_cur, Data data_can)
{
  bool verbose = false;
  if(verbose)
  {
    std::cout << "--------------------" << std::endl;
    std::cout << "Determine if node and observed children are the same" << std::endl;
  }
  if(data_cur.genetic_ids[node_cur] != data_can.genetic_ids[node_can] ||
     data_cur.sample_times[node_cur] != data_can.sample_times[node_can] ||
     data_cur.subtype_numbers[node_cur] != data_can.subtype_numbers[node_can])
  {
    return false;
  }
  
  std::vector<int> observed_child_nodes_cur = data_cur.ReturnObservedChildren(node_cur);
  std::vector<int> observed_child_nodes_can = data_can.ReturnObservedChildren(node_can);
  
  
  
  bool obs_child_nodes_equal = DetermineIfNodesAreEqual(observed_child_nodes_cur, observed_child_nodes_can, 
                                                        data_cur, data_can);
  if(verbose)
  {
    if(obs_child_nodes_equal)
    {
      std::cout << "Nodes are equal - observed children cur = ";
      PrintVector(observed_child_nodes_cur);
      std::cout << "observed children can = ";
      PrintVector(observed_child_nodes_can);
    }
    else
    {
      std::cout << "Nodes are not equal - observed children cur = ";
      PrintVector(observed_child_nodes_cur);
      std::cout << "observed children can = ";
      PrintVector(observed_child_nodes_can);
    }
    
  }
  return obs_child_nodes_equal;
}

bool DetermineIfConfigurationSame(int node_cur, int node_can, Data data_cur, Data data_can)
{
  bool verbose = false;
  if(verbose)
  {
    std::cout << "Determine if cur node " << node_cur << " and can node " << node_can 
              << " have the same configuration" << std::endl;
    std::cout << "data_cur.gen_source[node_cur] = " << data_cur.gen_source[node_cur]
              << ", data_can.gen_source[node_can] = " << data_can.gen_source[node_can] << std::endl;
  }
  // check if they are both imported
  if(data_cur.gen_source[node_cur] == -1 && data_can.gen_source[node_can] == -1)
  {
    // both nodes are imported, do a quick check
    bool current_nodes_equal = DetermineIfNodeAndObsChildNodesAreSame(node_cur, node_can,
                                                                      data_cur, data_can);
    return current_nodes_equal;
  }
  
  // otherwise find same root
  bool root_found = false;
  while(!root_found)
  {
    bool current_nodes_equal = DetermineIfNodeAndObsChildNodesAreSame(node_cur, node_can, data_cur, data_can);
    if(verbose)
    {
      std::cout << "node cur = " << node_cur << ", node can = " << node_can << ", nodes equal = "
                << current_nodes_equal << std::endl;
    }
    if(!current_nodes_equal) return false;
    
    // Nodes are equal, update nodes and try again
    node_cur = data_cur.gen_source[node_cur];
    node_can = data_can.gen_source[node_can];
    if(verbose)
    {
      std::cout << "  updated node cur = " << node_cur << ", updated node can = " << node_can << std::endl;
    }
    
    if(node_cur == -1 && node_can == -1)
    {
      // previous node are root, we can stop
      root_found = true;
    }
    else if(node_cur != -1 && node_can != -1)
    {
      // both nodes are valid nodes
      // check if the nodes are observed
      if(node_cur == node_can)
      {
        // nodes are equal
        if(data_cur.imputed_nodes[node_cur] == 0) return true;
      }
    }
    else
    {
      // one of the nodes is a root and the other isn't, 
      return false;
    }
  }
  return true;
}

bool NodesEqual(int node_i, int node_j,
                Data data_i, Data data_j)
{
  int gen_id_i = data_i.genetic_ids[node_i];
  int gen_id_j = data_j.genetic_ids[node_j];
  if(gen_id_i != gen_id_j) return false;
  
  double samp_time_i = data_i.sample_times[node_i];
  double samp_time_j = data_j.sample_times[node_j];
  if(samp_time_i != samp_time_j) return false;
  
  int subtype_i = data_i.subtype_numbers[node_i];
  int subtype_j = data_j.subtype_numbers[node_j];
  if(subtype_i != subtype_j) return false;
  
  int imp_type_i = data_i.imputed_nodes[node_i];
  int imp_type_j = data_j.imputed_nodes[node_j];
  if(imp_type_i != imp_type_j) return false;
  
  return true;
}

// Will check if two vectors of idx contain the same nodes
bool MultipleNodesEqual(std::vector<int> nodes_i, std::vector<int> nodes_j,
                        Data data_i, Data data_j)
{
  bool verbose = false;
  
  if(verbose)
  {
    std::cout << "Check if nodes_i = ";
    for(auto i : nodes_i)
    {
      std::cout << i << " ";
    }
    std::cout << " and nodes_j = ";
    for(auto j : nodes_j)
    {
      std::cout << j << " ";
    }
    std::cout << std::endl;
  }
  
  // Check vectors are same length
  if(nodes_i.size() != nodes_j.size()) return false;
  
  // Look at each node in set i, and try and find the corresponding node
  // in set j
  size_t num_nodes = nodes_i.size();
  
  // set to 1 if already found, otherwise zero
  std::vector<int> already_searched_j(num_nodes);
  for(auto node_i : nodes_i)
  {
    size_t counter = 0;
    for(size_t k = 0; k < num_nodes; k++)
    {
      if(already_searched_j[k] == 0)
      {
        int node_j = nodes_j[k];
        bool nodes_equal = NodesEqual(node_i, node_j, data_i, data_j);
        if(verbose)
        {
          std::cout << "Comparing node_i = " << node_i
                    << " and node_j = " << node_j
                    << ", are equal = " << nodes_equal
                    << std::endl;
        }
        if(nodes_equal)
        {
          already_searched_j[k] = 1;
          break;
        }
      }
      counter++;
    }
    if(counter == num_nodes) return false;
  }
  return true;
}

// Given two imputed nodes, are the configurations the same
bool SameConfiguration(int imp_idx_cur, int imp_idx_can,
                       Data data_cur, Data data_can)
{
  bool verbose = false;
  // either of the imp idxes may be -1, check first
  if(imp_idx_cur == -1) return false;
  if(imp_idx_can == -1) return false;
  
  if(verbose)
  {
    std::cout << "Determine if nodes " << imp_idx_cur << " and " << imp_idx_can
              << " have the same configuration" << std::endl;
  }
  
  // Check that both nodes are imputed
  if(data_cur.imputed_nodes[imp_idx_cur] == 0 ||
     data_can.imputed_nodes[imp_idx_can] == 0)
  {
    std::cout << std::endl << std::endl << "*** Error ***" << std::endl;
    std::cout << "imp_idx_cur = " << imp_idx_cur
              << ", imp_idx_can = " << imp_idx_can
              << ", data_cur.imputed_nodes[imp_idx_cur] = "
              << data_cur.imputed_nodes[imp_idx_cur] 
              << ", data_can.imputed_nodes[imp_idx_can] = "
              << data_can.imputed_nodes[imp_idx_can]
              << std::endl;
    stop("One of the nodes are observed");
  }
  
  
  
  // Check if the nodes themselves are the same
  bool imp_nodes_equal = NodesEqual(imp_idx_cur, imp_idx_can,
                                    data_cur, data_can);
  if(!imp_nodes_equal) return false;
  if(verbose)
  {
    std::cout << " Nodes themselves are equal, now determine other parts" << std::endl;
  }
  
  // By this point the nodes themselves are the same, now we should
  // determine the configuration depending on the type of node
  
  int imputation_type = data_cur.imputed_nodes[imp_idx_cur];
  int imputation_type_can = data_can.imputed_nodes[imp_idx_can];
  if(imputation_type != imputation_type_can) return false;
  std::vector<int> child_nodes_cur = WhichVec(imp_idx_cur, data_cur.gen_source);
  std::vector<int> child_nodes_can = WhichVec(imp_idx_can, data_can.gen_source);
  
  switch(imputation_type)
  {
  case 1:
  {
    // Interior node
    
    // Determine the path of imputed nodes
    int cur_node = imp_idx_cur;
    int cur_parent = data_cur.gen_source[cur_node];
    std::vector<int> path_cur = {cur_node};
    while(data_cur.imputed_nodes[cur_parent] > 0)
    {
      path_cur.push_back(cur_parent);
      cur_node = cur_parent;
      cur_parent = data_cur.gen_source[cur_node];
    }
    
    int can_node = imp_idx_can;
    int can_parent = data_can.gen_source[can_node];
    std::vector<int> path_can = {can_node};
    while(data_can.imputed_nodes[can_parent] > 0)
    {
      path_can.push_back(can_parent);
      can_node = can_parent;
      can_parent = data_can.gen_source[can_node];
    }
    
    if(verbose)
    {
      std::cout << " Looking an an interior node with path_cur = ";
      for(auto i : path_cur)
      {
        std::cout << i << " ";
      }
      std::cout << " and path_can = ";
      for(auto i : path_can)
      {
        std::cout << i << " ";
      }
      std::cout << std::endl;
    }
    
    // Now start at the nodes in the current and candidate path, and check that both
    // the nodes are the same, the parent is the same, and the observed children
    // are the same
    
    if(path_can.size() != path_cur.size()) return false;
    
    for(size_t i = 0; i < path_can.size(); i++)
    {
      // Check nodes are equal
      cur_node = path_cur[i];
      can_node = path_can[i];
      bool imp_nodes_equal = NodesEqual(cur_node, can_node,
                                        data_cur, data_can);
      if(!imp_nodes_equal) return false;
      
      // Check parent nodes are equal
      cur_parent = data_cur.gen_source[cur_node];
      can_parent = data_can.gen_source[can_node];
      bool parent_nodes_equal = NodesEqual(cur_parent, can_parent,
                                           data_cur, data_can);
      if(!parent_nodes_equal) return false;
      
      // Check observed child nodes are equal
      std::vector<int> obs_ch_nodes_cur = data_cur.ReturnObservedChildren(cur_node);
      std::vector<int> obs_ch_nodes_can = data_can.ReturnObservedChildren(can_node);
      bool obs_child_nodes_equal = MultipleNodesEqual(obs_ch_nodes_cur, obs_ch_nodes_can,
                                                      data_cur, data_can);
      if(!obs_child_nodes_equal) return false;
    }
    return true;
  }
  case 2:
  {
    // Exterior node
    // observed and direct observed child nodes must be the same
    //return false;
    std::vector<int> obs_ch_nodes_cur = data_cur.ReturnObservedChildren(imp_idx_cur);
    std::vector<int> obs_ch_nodes_can = data_can.ReturnObservedChildren(imp_idx_can);
    bool obs_child_nodes_equal = MultipleNodesEqual(obs_ch_nodes_cur, obs_ch_nodes_can,
                                                    data_cur, data_can);
    if(!obs_child_nodes_equal) return false;
    
    // now check the direct nodes
    std::vector<int> dir_obs_ch_nodes_cur;
    std::vector<int> dir_obs_ch_nodes_can;
    
    for(size_t i = 0; i < obs_ch_nodes_cur.size(); i++)
    {
      int ch_node_cur = obs_ch_nodes_cur[i];
      int ch_node_can = obs_ch_nodes_can[i];
      if(data_cur.gen_source[ch_node_cur] == imp_idx_cur)
      {
        dir_obs_ch_nodes_cur.push_back(ch_node_cur);
      }
      
      if(data_can.gen_source[ch_node_can] == imp_idx_can)
      {
        dir_obs_ch_nodes_can.push_back(ch_node_can);
      }
    }
    
    bool dir_obs_child_nodes_equal = MultipleNodesEqual(dir_obs_ch_nodes_cur,
                                                        dir_obs_ch_nodes_can,
                                                        data_cur, data_can);
    return dir_obs_child_nodes_equal;
  }
  case 3:
    // Upper node
    // Child nodes should also be equal
    bool child_nodes_equal = MultipleNodesEqual(child_nodes_cur, child_nodes_can,
                                                data_cur, data_can);
    return child_nodes_equal;
  }
  return false;
}

// Given a node, return all the children that are not imputed
std::vector<int> Data::ReturnObservedChildren(int node)
{
  bool child_nodes_found = false;
  std::vector<int> observed_child_nodes;
  std::vector<int> nodes_to_search = {node};
  
  while(!child_nodes_found)
  {
    std::vector<int> new_nodes_to_search;
    for(size_t i = 0; i < nodes_to_search.size(); i++)
    {
      int current_node = nodes_to_search[i];
      std::vector<int> child_nodes = WhichVec(current_node, gen_source);
      for(size_t j = 0; j < child_nodes.size(); j++)
      {
        int current_child = child_nodes[j];
        if(imputed_nodes[current_child] > 0)
        {
          // The current child is imputed, therefore add to the list of nodes to search
          new_nodes_to_search.push_back(current_child);
        }
        else
        {
          observed_child_nodes.push_back(current_child);
        }
      }
    }
    nodes_to_search = new_nodes_to_search;
    if(nodes_to_search.size() == 0)
    {
      child_nodes_found = true;
    }
  }
  
  std::sort(observed_child_nodes.begin(), observed_child_nodes.end());
  
  return observed_child_nodes;
}

bool Data::DoesTargetHaveSequenceGreaterThanTime(int target, double time, int subtype)
{
  for(size_t i = 0; i < genetic_ids.size(); i++)
  {
    if(target == genetic_ids[i] && subtype == subtype_numbers[i] && sample_times[i] > time)
    {
      return true;
    }
  }
  return false;
}

// Given a target and a subtype, determine if they need a node imputed
// at the time of colonisation
bool Data::DoesNodeNeedToBeImputed(int target, int subtype)
{
  double target_col_time = t_e[target];
  
  // First check at the time of colonisation
  for(size_t i = 0; i < genetic_ids.size(); i++)
  {
    if(target == genetic_ids[i] && target_col_time == sample_times[i] && subtype == subtype_numbers[i])
    {
      return false;
    }
  }
  
  bool forward_sequence_found_target = (*this).DoesTargetHaveSequenceGreaterThanTime(target, target_col_time, subtype);
  
  if(!forward_sequence_found_target)
  {
    std::vector<int> primary_targets = WhichVec(target, source);
    bool still_searching_secondary = true;
    while(still_searching_secondary)
    {
      if(primary_targets.size() == 0)
      {
        // no one found in this branch, node does not need to be imputed
        return false;
      }
      std::vector<int> secondary_targets;
      for(size_t i = 0; i<primary_targets.size(); i++)
      {
        int current_target = primary_targets[i];
        forward_sequence_found_target = (*this).DoesTargetHaveSequenceGreaterThanTime(current_target, target_col_time, 
                                         subtype);
        if(!forward_sequence_found_target)
        {
          // nothing found for this target, store his children to check later
          std::vector<int> current_target_children = WhichVec(current_target, source);
          for(size_t j = 0; j < current_target_children.size(); j++)
          {
            secondary_targets.push_back(current_target_children[j]);
          }
        }
        else
        {
          // found a sequence, we can exit early
          still_searching_secondary = false;
          break;
        }
      }
      primary_targets = secondary_targets;
    }
  }
  
  if(!forward_sequence_found_target) return false;
  
  // now do the same but with the infector
  bool forward_sequence_found_infector = false;
  
  // check if the infector has a swab at the time of infection, OR another infection with a swab at that time
  int infector = source[target];
  for(size_t i = 0; i < genetic_ids.size(); i++)
  {
    if(infector == genetic_ids[i] && target_col_time == sample_times[i] && subtype == subtype_numbers[i])
    {
      return false;
    }
  }
  
  std::vector<int> infector_sources = WhichVec(infector, source);
  std::vector<int> primary_targets;
  for(size_t i = 0; i < infector_sources.size(); i++)
  {
    int current_source = infector_sources[i];
    if(current_source != target)
    {
      double current_source_col_time = t_e[current_source];
      // only look at individuals who are colonised at the same time or after the target
      if(current_source_col_time >= target_col_time)
      {
        if(current_source_col_time == target_col_time)
        {
          for(size_t j = 0; j < genetic_ids.size(); j++)
          {
            if(current_source == genetic_ids[j] && target_col_time == sample_times[j] && subtype == subtype_numbers[j])
            {
              return false;
            }
          }
        }
        primary_targets.push_back(current_source);
      }
    }
  }
  
  
  
  
  forward_sequence_found_infector = (*this).DoesTargetHaveSequenceGreaterThanTime(infector, target_col_time, 
                                     subtype);
  
  // look through children
  bool still_searching_secondary = true;
  if(forward_sequence_found_infector)
  {
    still_searching_secondary = false;
  }
  while(still_searching_secondary)
  {
    if(primary_targets.size()==0)
    {
      // no one found in this branch, node does not need to be imputed
      return false;
    }
    std::vector<int> secondary_targets;
    for(size_t i = 0; i<primary_targets.size(); i++)
    {
      int current_target = primary_targets[i];
      forward_sequence_found_infector = (*this).DoesTargetHaveSequenceGreaterThanTime(current_target, target_col_time, 
                                         subtype);
      if(!forward_sequence_found_infector)
      {
        // nothing found for this target, store his children to check later
        std::vector<int> current_target_children = WhichVec(current_target, source);
        for(size_t j = 0; j<current_target_children.size(); j++)
        {
          secondary_targets.push_back(current_target_children[j]);
        }
      }
      else
      {
        // found a sequence, we can exit early
        still_searching_secondary = false;
        break;
      }
    }
    primary_targets = secondary_targets;
  }
  
  if(forward_sequence_found_infector && forward_sequence_found_target)
  {
    return true;
  }
  return false;
}




int Data::ReturnExteriorTarget(int node)
{
  int current_id = genetic_ids[node];
  int current_subtype = subtype_numbers[node];
  int current_source = source[current_id];
  
  while(current_source != -1)
  {
    bool impute_node = DoesNodeNeedToBeImputed(current_id, current_subtype);
    if(impute_node)
    {
      return current_id;
      // This should put the imputed node in the source rather than the target
      // it should make no difference but this is to ensure consistency with
      // the way nodes are imputed in the first place
      //eturn current_source;
    }
    else
    {
      current_id = current_source;
      current_source = source[current_id];
    }
  }
  return -1;
}

void Data::ReturnExteriorNodes()
{
  std::vector<int> no_parent_idx = WhichVec(-1, gen_source);
  bool verbose = false;
  if(verbose)
  {
    std::cout << " Attempt to return exterior nodes" << std::endl;
    std::cout << " gen source = ";
    PrintVector(gen_source);
    std::cout << " no parent idx = ";
    PrintVector(no_parent_idx);
  }
  
  
  std::vector<int> subtype_numbers_cur(subtype_numbers);
  // only include genetic data with no parent sequence
  genetic_ids = subset(no_parent_idx, genetic_ids);
  sample_times = subset(no_parent_idx, sample_times);
  subtype_numbers = subset(no_parent_idx, subtype_numbers);
  
  size_t num_seqs = genetic_ids.size();
  //std::vector<int> nodes_to_impute_(num_seqs);
  //std::vector<int> times_to_impute_(num_seqs);
  //std::vector<int> subtypes_to_impute_(num_seqs);
  std::vector<int> nodes_to_impute_;
  std::vector<double> times_to_impute_;
  std::vector<int> subtypes_to_impute_;
  
  std::vector<int> exterior_targets(num_seqs);
  
  //#pragma omp parallel for schedule(dynamic) if(num_seqs > 200)
  for(size_t i = 0; i < num_seqs; i++)
  {
    int exterior_node = ReturnExteriorTarget(i);
    exterior_targets[i] = exterior_node;
    //int current_subtype = subtype_numbers[i];
    //std::cout << "node " << i << " has exterior node " << exterior_node  << ", from ID = "
    //          << source[exterior_node] << "at time t = " << t_c[exterior_node] 
    //          << ", with subtype = " << current_subtype << std::endl;
  }
  
  // Now check if there are duplicate entries
  for(size_t i = 0; i < num_seqs; i++)
  {
    int cur_exterior_target = exterior_targets[i];
    if(cur_exterior_target != -1)
    {
      int cur_id = source[cur_exterior_target];
      double cur_time = t_e[cur_exterior_target];
      int cur_subtype = subtype_numbers[i];
      
      bool duplicate_found = false;
      for(size_t j = 0; j < nodes_to_impute_.size(); j++)
      {
        if(cur_id == nodes_to_impute_[j] &&
           cur_time == times_to_impute_[j] &&
           cur_subtype == subtypes_to_impute_[j])
        {
          duplicate_found = true;
          break;
        }
      }
      
      if(!duplicate_found)
      {
        nodes_to_impute_.push_back(cur_id);
        times_to_impute_.push_back(cur_time);
        subtypes_to_impute_.push_back(cur_subtype);
      }
    }
  }
  
  if(verbose)
  {
    std::cout << " Exterior nodes to impute = ";
    PrintVector(nodes_to_impute_);
    std::cout << " Exterior times to impute = ";
    PrintVector(times_to_impute_);
    std::cout << " Exterior subtypes to impute = ";
    PrintVector(subtypes_to_impute_);
  }
  
  // ensure there are no duplicate entries
  
  nodes_to_impute = nodes_to_impute_;
  times_to_impute = times_to_impute_;
  subtypes_to_impute = subtypes_to_impute_;
}



// returns surrounding nodes for an imputed node
std::vector<int> Data::ReturnSurroundingNodes(int node)
{
  if(imputed_nodes[node] == 0)
  {
    stop("trying to return surrounding nodes for an observed node");
  }
  
  // first add observed children
  std::vector<int> surrounding_nodes = ReturnObservedChildren(node);
  
  // then check if parent needs to be added
  int parent = gen_source[node];
  if(parent != -1)
  {
    surrounding_nodes.push_back(parent);
  }
  return surrounding_nodes;
}

// with root at position zero
std::vector<int> Data::ReturnPathToRoot(int node)
{
  bool verbose = false;
  if(verbose)
  {
    std::cout << "Attempt to return path to root for node " << node
              << std::endl;
  }
  std::vector<int> path = {node};
  int current_node = gen_source[node];
  while(current_node != -1)
  {
    path.push_back(current_node);
    current_node = gen_source[current_node];
  }
  std::reverse(path.begin(),path.end());
  if(verbose)
  {
    std::cout << "Path = ";
    PrintVector(path);
  }
  return path;
}

// Returns the lowest common ancestor between nodes i and j
int Data::ReturnLCA(int node_i, int node_j)
{
  bool verbose = false;
  if(verbose)
  {
    std::cout << "---Attempt to return LCA for nodes " << node_i
              << " and " << node_j << std::endl;
  }
  std::vector<int> path_i = ReturnPathToRoot(node_i);
  std::vector<int> path_j = ReturnPathToRoot(node_j);
  
  //std::reverse(path_i.begin(),path_i.end());
  //std::reverse(path_j.begin(),path_j.end());
  
  if(verbose)
  {
    std::cout << "Path i = ";
    PrintVector(path_i);
    std::cout << "Path j = ";
    PrintVector(path_j);
  }
  
  
  int min_len = std::min(path_i.size(),path_j.size());
  for(int k = 0; k < min_len; k++)
  {
    if(verbose)
    {
      std::cout << "k = " << k << ", path_i[k] = " << path_i[k] << ", path_j[k] = "
                << path_j[k] << std::endl;
    }
    if(path_i[k] != path_j[k])
    {
      return path_i[k-1];
    }
  }
  return -1;
}

int Data::CalculateDistanceBetweenNodes(int node_i, int node_j)
{
  bool verbose = false;
  //if(node_i == 3 && node_j == 8) verbose = true;
  if(verbose)
  {
    std::cout << "Attempt to calculate the distance between nodes "
              << node_i << " and " << node_j << " with gen source vector ";
    PrintVector(gen_source);
  }
  std::vector<int> path_i = ReturnPathToRoot(node_i);
  std::vector<int> path_j = ReturnPathToRoot(node_j);
  int dist_sum = 0;
  bool same_pathway = false;
  if(verbose)
  {
    std::cout << "path_i = ";
    PrintVector(path_i);
    std::cout << "path_j = ";
    PrintVector(path_j);
  }
  if(path_i.size() > path_j.size())
  {
    // look for j in path_i
    for(auto i : path_i)
    {
      if(node_j == i)
      {
        if(verbose)
        {
          std::cout << "Node j in contained in path_i" << std::endl;
        }
        same_pathway = true;
        break;
      }
    }
  }
  else
  {
    // look for i in path_j
    for(auto j : path_j)
    {
      if(node_i == j)
      {
        if(verbose)
        {
          std::cout << "Node i in contained in path_j" << std::endl;
        }
        same_pathway = true;
        break;
      }
    }
  }
  
  if(same_pathway)
  {
    // They are in the same pathway so just add nodes along the path
    if(path_i.size() > path_j.size())
    {
      // j must be contained in path_i
      int target = node_i;
      int parent = gen_source[node_i];
      dist_sum += gen_matrix(target, parent);
      while(parent != node_j)
      {
        target = parent;
        parent = gen_source[target];
        dist_sum += gen_matrix(target, parent);
      }
      return dist_sum;
    }
    else
    {
      // i must be contained in path_j
      int target = node_j;
      int parent = gen_source[node_j];
      dist_sum += gen_matrix(target, parent);
      while(parent != node_i)
      {
        target = parent;
        parent = gen_source[target];
        dist_sum += gen_matrix(target, parent);
      }
      return dist_sum;
    }
  }
  else
  {
    // Find the lowest common ancestor
    int LCA = ReturnLCA(node_i, node_j);
    if(verbose)
    {
      std::cout << "LCA = " << LCA << std::endl;
    }
    // Sum in each path to the LCA
    int target = node_i;
    int parent = gen_source[node_i];
    if(verbose)
    {
      std::cout << " target = " << target << ", parent = "
                << parent << ", dist = " << gen_matrix(target, parent)
                << std::endl;
    }
    dist_sum += gen_matrix(target, parent);
    while(parent != LCA)
    {
      target = parent;
      parent = gen_source[target];
      dist_sum += gen_matrix(target, parent);
      if(verbose)
      {
        std::cout << " target = " << target << ", parent = "
                  << parent << ", dist = " << gen_matrix(target, parent)
                  << std::endl;
      }
    }
    
    target = node_j;
    parent = gen_source[node_j];
    dist_sum += gen_matrix(target, parent);
    if(verbose)
    {
      std::cout << " target = " << target << ", parent = "
                << parent << ", dist = " << gen_matrix(target, parent)
                << std::endl;
    }
    while(parent != LCA)
    {
      target = parent;
      parent = gen_source[target];
      dist_sum += gen_matrix(target, parent);
      if(verbose)
      {
        std::cout << " target = " << target << ", parent = "
                  << parent << ", dist = " << gen_matrix(target, parent)
                  << std::endl;
      }
    }
    return dist_sum;
  }
  stop("Error calculating distance sum, should not get this far");
}


// Return the nodes to impute when first initialised
void Data::ReturnNodesToImputeInitial()
{
  bool verbose = false;
  if(verbose)
  {
    std::cout << "Attempt to return the initial nodes to impute"
              << std::endl;
  }
  std::vector<int> nodes_to_impute_;
  std::vector<double> times_to_impute_;
  std::vector<int> subtypes_to_impute_;
  
  
  // Store the IDs of all colonisation events
  std::vector<int> colonisation_events;
  std::vector<double> colonisation_times;
  for(int i = 0; i < N; i++)
  {
    if(source[i] >= 0)
    {
      colonisation_events.push_back(i);
      colonisation_times.push_back(t_e[i]);
    }
  }
  
  if(verbose)
  {
    std::cout << "colonisation_events = ";
    PrintVector(colonisation_events);
    std::cout << "colonisation_times = ";
    PrintVector(colonisation_times);
  }
  
  // Sort the colonisation events by times in increasing order
  colonisation_events = sort_two(colonisation_events, colonisation_times);
  
  // Create temp data
  Data data_temp(*this);
  
  std::vector<int> unique_subtype_numbers = sort_unique(subtype_numbers);
  
  for(size_t i = 0; i < colonisation_events.size(); i++)
  {
    for(size_t j = 0; j < unique_subtype_numbers.size(); j++)
    {
      int current_subtype = unique_subtype_numbers[j];
      int current_colonisation_event = colonisation_events[i];
      bool impute_node = data_temp.DoesNodeNeedToBeImputed(current_colonisation_event, current_subtype);
      std::cout << "Current event = " << current_colonisation_event << ", current subtype = "
                << current_subtype << ", impute node = " << impute_node << std::endl;
      if(impute_node)
      {
        // Check for a duplicate
        bool duplicate_node = false;
        for(size_t i = 0; i < nodes_to_impute.size(); i++)
        {
          if(nodes_to_impute[i] == data_temp.source[current_colonisation_event] &&
             times_to_impute[i] == data_temp.t_e[current_colonisation_event] &&
             subtypes_to_impute[i] == current_subtype)
          {
            duplicate_node = true;
            break;
          }
        }
        
        if(!duplicate_node)
        {
          // If the node is not a duplicate
          data_temp.genetic_ids.push_back(source[current_colonisation_event]);
          data_temp.sample_times.push_back(t_e[current_colonisation_event]);
          data_temp.subtype_numbers.push_back(current_subtype);
          
          nodes_to_impute_.push_back(source[current_colonisation_event]);
          times_to_impute_.push_back(t_e[current_colonisation_event]);
          subtypes_to_impute_.push_back(current_subtype);
          
          if(verbose)
          {
            std::cout << "Node to impute from genetic id = "
                      << source[current_colonisation_event] << " at time "
                      << t_e[current_colonisation_event] << " with subtype "
                      << current_subtype << std::endl;
          }
        }
      }
    }
  }
  
  
  
  // Sort the nodes to ensure a unique ordering
  std::vector<int> sorted_indices = sort_three(nodes_to_impute_, times_to_impute_, subtypes_to_impute_);
  
  if(verbose)
  {
    std::cout << "Before sort, times_to_impute_ = ";
    PrintVector(times_to_impute_);
  }
  nodes_to_impute = subset(sorted_indices, nodes_to_impute_);
  times_to_impute = subset(sorted_indices, times_to_impute_);
  subtypes_to_impute = subset(sorted_indices, subtypes_to_impute_);
  
  if(verbose)
  {
    std::cout << "nodes_to_impute = ";
    PrintVector(nodes_to_impute);
    std::cout << "times_to_impute = ";
    PrintVector(times_to_impute);
    std::cout << "subtypes_to_impute";
    PrintVector(subtypes_to_impute);
  }
}


void Data::InitialiseImputedNodes()
{
  bool verbose = false;
  if(verbose)
  {
    std::cout << "Attempt to initialise the imputed nodes"
              << std::endl;
  }
  
  ReturnNodesToImputeInitial();
  if(verbose)
  {
    std::cout << "Returned initial nodes to impute" << std::endl;
  }
  
  int observed_length = genetic_ids.size();
  // Update genetic data with imputed
  if(nodes_to_impute.size() > 0)
  {
    genetic_ids = CombineVectors(genetic_ids, nodes_to_impute);
    sample_times = CombineVectors(sample_times, times_to_impute);
    subtype_numbers = CombineVectors(subtype_numbers, subtypes_to_impute);
  }
  //CalculateGenSourceVector_parallel();
  
  if(verbose)
  {
    std::cout << "Updated genetic vectors" << std::endl;
  }
  
  // Update the imputed_nodes vector 
  int updated_length = genetic_ids.size();
  std::vector<int> imputed_nodes_(updated_length);
  
  
  for(int i = observed_length; i < updated_length; i++)
  {
    imputed_nodes_[i] = 1;
  }
  imputed_nodes = imputed_nodes_;
  if(verbose)
  {
    std::cout << "Imputed nodes length = " << imputed_nodes.size() << ", imputed nodes = ";
    PrintVector(imputed_nodes);
    std::cout << "genetic_ids = ";
    PrintVector(genetic_ids);
    std::cout << "sample_times = ";
    PrintVector(sample_times);
    std::cout << "subtype_numbers = ";
    PrintVector(subtype_numbers);
  }
  
  
  // Create a new matrix and fill in the entries
  
  IntegerMatrix gen_mat_out(updated_length, updated_length);
  for(int i = 0; i < updated_length; i++)
  {
    for(int j = 0; j < updated_length; j++)
    {
      if(i < observed_length && j < observed_length)
      {
        gen_mat_out(i,j) = gen_matrix(i,j);
      }
      else
      {
        gen_mat_out(i,j) = 0;
      }
      
      //gen_mat_out(j,i) = genetic_matrix(j,i);
    }
  }
  
  if(verbose)
  {
    std::cout << "Initialise new matrix and copy observed entries" << std::endl;
    std::cout << "gen Matrix rows = " << gen_mat_out.nrow() << ", columns = " << gen_mat_out.ncol() << std::endl;
  }
  
  CalculateGenSourceVector();
  if(verbose)
  {
    std::cout << "Genetic source vector calculated" << std::endl;
  }
  
  // Copy the parent of all imputed nodes
  std::vector<int> order_of_nodes_to_impute = DetermineOrderOfImputation();
  if(verbose)
  {
    std::cout << "Order of imputation calculated - ";
    PrintVector(order_of_nodes_to_impute);
  }
  
  for(size_t ii = 0; ii < order_of_nodes_to_impute.size(); ii++)
  {
    int node_to_impute = order_of_nodes_to_impute[ii];
    int parent_node = gen_source[node_to_impute];
    std::vector<int> observed_child_nodes = ReturnObservedChildren(node_to_impute);
    if(verbose)
    {
      std::cout << "Imputing node " << node_to_impute << " with parent " << parent_node
                << " and observed child nodes ";
      PrintVector(observed_child_nodes);
    }
    
    if(parent_node != -1)
    {
      // Interior node, copy the parent
      
      // Find the first observed parent
      int previous_parent = parent_node;
      int observed_parent = parent_node;
      bool observed_parent_found = false;
      while(!observed_parent_found)
      {
        std::cout << " observed_parent = " << observed_parent << ", imputed_nodes[observed_parent] = " 
                  << imputed_nodes[observed_parent] << std::endl;
        if(observed_parent != -1)
        {
          if(imputed_nodes[observed_parent] == 0)
          {
            observed_parent_found = true;
          }
          else
          {
            previous_parent = observed_parent;
            observed_parent = gen_source[observed_parent];
          }
        }
        else
        {
          observed_parent = previous_parent;
          observed_parent_found = true;
        }
        
      }
      if(verbose) std::cout << "Observed parent = " << observed_parent << std::endl;
      gen_mat_out(node_to_impute, parent_node) = 0;
      
      gen_mat_out(node_to_impute, parent_node) = 0;
      
      
      for(size_t i = 0; i < observed_child_nodes.size(); i++)
      {
        int current_observed_child = observed_child_nodes[i];
        //gen_mat_out(node_to_impute, current_observed_child) = genetic_matrix(observed_parent, current_observed_child);
        //gen_mat_out(current_observed_child, node_to_impute) = genetic_matrix(observed_parent, current_observed_child);
        gen_mat_out(node_to_impute, current_observed_child) = gen_mat_out(observed_parent, current_observed_child);
        gen_mat_out(current_observed_child, node_to_impute) = gen_mat_out(observed_parent, current_observed_child);
        
      }
    }
    else
    {
      // Exterior node, copy smalled observed child node
      int min_i = observed_child_nodes[0];
      int min_j = observed_child_nodes[1];
      int min_distance = gen_mat_out(min_i, min_j);
      
      
      
      // Look for smallest distance
      for(size_t j = 1; j < observed_child_nodes.size(); j++)
      {
        for(size_t i = 0; i < j; i++)
        {
          int cur_i = observed_child_nodes[i];
          int cur_j = observed_child_nodes[j];
          
          int current_distance = gen_matrix(cur_i, cur_j);
          
          if(current_distance < min_distance)
          {
            min_distance = current_distance;
            min_i = cur_i;
            min_j = cur_j;
          }
          
          
        }
      }
      
      // Copy the entries of the min child
      for(int i = 0; i < updated_length; i++)
      {
        gen_mat_out(i, node_to_impute) = gen_mat_out(i, min_i);
        gen_mat_out(node_to_impute, i) = gen_mat_out(i, min_i);
        
        
      }
    }
  }
  
  if(verbose)
  {
    std::cout << "Immediate distances copied" << std::endl;
  }
  
  std::vector<int> imported_idx = WhichVec(-1, gen_source);
  std::vector<int> coltime_distances_(imported_idx.size());
  
  coltime_distances = coltime_distances_;
  
  
  
  gen_matrix = gen_mat_out;
}

// if true we need to calculate the contribution
// otherwise false we don't
bool CalculateContribution(int imp_idx_cur,
                           Data data_cur, Data data_can)
{
  bool verbose = false;
  int current_id = data_cur.genetic_ids[imp_idx_cur];
  double current_time = data_cur.sample_times[imp_idx_cur];
  int current_subtype = data_cur.subtype_numbers[imp_idx_cur];
  
  int imp_idx_can = data_can.FindGeneticLoc(current_id, current_time, current_subtype);
  
  if(verbose)
  {
    std::cout << "--Determine if the sequence at idx " << imp_idx_cur 
              << " from ID = " << current_id
              << " at time t = " << current_time
              << " with imp_idx_can = " << imp_idx_can
              << " exists in the other data set"
              << std::endl;
  }
  
  
  // Loc is not found in the candidate data set, need to calculate contribution
  if(imp_idx_can == -1) return true;
  
  bool same_config = DetermineIfConfigurationSame(imp_idx_cur, imp_idx_can, 
                                                  data_cur, data_can);
  return !same_config;
}


// Calculates the contribution of the imputed idx
// for the reverse process where idx is an upper node
void Data::CalculateContribRevUpper(int imp_idx, Data data_can, bool verbose)
{
  // distances were simulated from the mutation model
  std::vector<int> child_nodes = WhichVec(imp_idx, gen_source);
  if(verbose)
  {
    std::cout << " Removing upper node " << imp_idx << " from ID = "
              << genetic_ids[imp_idx] << " at t = "
              << sample_times[imp_idx] << " with child nodes = ";
    PrintVector(child_nodes);
  }
  
  int gen_id_cur = genetic_ids[imp_idx];
  double samp_time_cur = sample_times[imp_idx];
  int subtype_cur = subtype_numbers[imp_idx];
  
  int imp_idx_can = data_can.FindGeneticLoc(gen_id_cur, samp_time_cur, 
                                            subtype_cur);
  bool same_config = SameConfiguration(imp_idx, imp_idx_can,
                                       (*this), data_can);
  
  if(same_config)
  {
    if(verbose)
    {
      std::cout << " Same configuration as the candidate data set, do not calculate contribution"
                << std::endl << std::endl;
    }
    return;
  }
  
  
  for(auto ch_node : child_nodes)
  {
    double time_diff = sample_times[ch_node] - sample_times[imp_idx];
    double mutation_rate = parameters[6];
    double mut_prob = CalculateMutationProbability(mutation_rate, time_diff);
    int draw = gen_matrix(imp_idx, ch_node);
    double contribution = R::dbinom(draw, sequence_length, mut_prob, 1);
    genetic_contribution += contribution;
    if(verbose)
    {
      std::cout << "  Current child node = " << ch_node << " from ID = "
                << genetic_ids[ch_node] << " at t = " 
                << sample_times[ch_node] << ", draw = " << draw
                << ", contribution = " << contribution
                << std::endl;
    }
  }
  if(verbose) std::cout << std::endl;
}

// Calculates the contribution of the imputed idx
// for the reverse process, where the idx is an exterior node
void Data::CalculateContribRevExterior(int imp_idx, Data data_can, bool verbose)
{
  std::vector<int> observed_child_nodes = ReturnObservedChildren(imp_idx);
  int min_distance = sequence_length;
  
  std::vector<int> direct_observed_child_nodes;
  for(auto node : observed_child_nodes)
  {
    if(gen_source[node] == imp_idx)
    {
      direct_observed_child_nodes.push_back(node);
    }
  }
  
  // node i should be directly connected, j will be an observed child
  int min_i = -1;
  int min_j = -1;
  
  int gen_id_cur = genetic_ids[imp_idx];
  double samp_time_cur = sample_times[imp_idx];
  int subtype_cur = subtype_numbers[imp_idx];
  int imp_idx_can = data_can.FindGeneticLoc(gen_id_cur, samp_time_cur, 
                                            subtype_cur);
  
  
  if(verbose)
  {
    std::cout << " Removing exterior node " << imp_idx 
              << " from ID = " << genetic_ids[imp_idx]
              << " at t = " << sample_times[imp_idx] 
              << ", observed child nodes = ";
    PrintVector(observed_child_nodes);
    std::cout << " direct observed child nodes = ";
    for(auto node : direct_observed_child_nodes)
    {
      std::cout << node << " ";
    }
    std::cout << std::endl;
  }
  
  bool same_config = SameConfiguration(imp_idx, imp_idx_can,
                                       (*this), data_can);
  if(same_config)
  {
    if(verbose)
    {
      std::cout << " Same configuration as the candidate data set, do not calculate contribution"
                << std::endl << std::endl;
    }
    return;
  }
  
  for(size_t i = 0; i < direct_observed_child_nodes.size(); i++)
  {
    for(size_t j = 0; j < observed_child_nodes.size(); j++)
    {
      int node_i = direct_observed_child_nodes[i];
      int node_j = observed_child_nodes[j];
      if(node_i != node_j)
      {
        int current_distance = gen_matrix(node_i, node_j);
        if(verbose && false)
        {
          std::cout << "node_i = " << node_i << ", node_j = "
                    << node_j << ", current_distance = "
                    << current_distance << std::endl;
        }
        // Check the T distances
        if(current_distance < min_distance)
        {
          min_distance = current_distance;
          min_i = node_i;
          min_j = node_j;
        }
      }
    }
  }
  
  // Calculate time ratios 
  double upper = sample_times[min_i] - sample_times[imp_idx];
  double lower = sample_times[min_i] + sample_times[min_j] - 
    2*sample_times[imp_idx];
  double time_ratio = upper/lower;
  if(time_ratio < 0 || time_ratio > 1)
  {
    std::cout << "Time ratio = " << time_ratio << std::endl;
    stop("invalid time ratio");
  }
  
  
  
  int observed_distance = gen_matrix(imp_idx, min_i);
  double contribution = R::dbinom(observed_distance, min_distance, time_ratio, 1);
  genetic_contribution += contribution;
  
  if(verbose)
  {
    std::cout << "  Different config, calculate contribution - "
              << " min_i = " << min_i << ", parent_distance = "
              << observed_distance << ", time_ratio = "
              << time_ratio << ", contribution = "
              << contribution << ", observed child nodes = ";
    PrintVector(observed_child_nodes);
    std::cout << std::endl;
  }
}

// Calculates the contribution of the imputed idx
// for the reverse process, where the idx is an interior node
void Data::CalculateContribRevInterior(int imp_idx, Data data_can, bool verbose)
{
  // CHeck if the node exists in the can data set
  int gen_id_cur = genetic_ids[imp_idx];
  double samp_time_cur = sample_times[imp_idx];
  int subtype_cur = subtype_numbers[imp_idx];
  int imp_idx_can = data_can.FindGeneticLoc(gen_id_cur, samp_time_cur, subtype_cur);
  
  bool same_config = SameConfiguration(imp_idx, imp_idx_can, (*this), data_can);
  int parent_node = gen_source[imp_idx];
  std::vector<int> observed_child_nodes = ReturnObservedChildren(imp_idx); 
  
  if(verbose)
  {
    std::cout << " Removing interior node " << imp_idx << " from ID = "
              << genetic_ids[imp_idx] << ", at t = "
              << sample_times[imp_idx] << ", with parent node = "
              << parent_node << " and observed child nodes = ";
    PrintVector(observed_child_nodes);
  }
  
  
  if(same_config)
  {
    if(verbose)
    {
      std::cout << "  Same configuration as the candidate data set, do not calculate contribution"
                << std::endl << std::endl;
    }
    return;
  }
  
  
  int min_child = observed_child_nodes[0];
  double min_time = sample_times[min_child];
  int min_distance = gen_matrix(parent_node, min_child);
  

  
  for(size_t i = 1; i < observed_child_nodes.size(); i++)
  {
    int current_child = observed_child_nodes[i];
    int current_distance = gen_matrix(parent_node, current_child);
    if(current_distance <= min_distance)
    {
      min_distance = current_distance;
      double current_time = sample_times[current_child];
      if(current_time < min_time)
      {
        min_time = current_time;
        min_child = current_child;
      }
    }
  }
  
  // Now perform the triangle inequality checks
  int triangle_min = ReturnMinDistanceTriangle(imp_idx, gen_matrix);
  min_distance = std::min(min_distance, triangle_min);
  
  double upper = sample_times[imp_idx] - sample_times[parent_node];
  double lower = sample_times[min_child] - sample_times[parent_node];
  double time_ratio = upper/lower;
  if(time_ratio < 0 || time_ratio > 1)
  {
    std::cout << "Upper = " << upper << ", lower = " << lower
              << ", ratio = " << upper/lower << std::endl;
    stop("Invalid time ratio");
  }
  
  int observed_distance = gen_matrix(imp_idx, parent_node);
  double contribution = R::dbinom(observed_distance, min_distance, time_ratio, 1);
  genetic_contribution += contribution;
  
  if(verbose)
  {
    std::cout << "  Different config, min_distance = "
              << min_distance << ", parent_distance = "
              << observed_distance << ", time_ratio = "
              << time_ratio << ", contribution = "
              << contribution
              << std::endl << std::endl;
  }
}

IntegerMatrix InitialiseCandidateMatrix(Data data_can)
{
  size_t length_can = data_can.genetic_ids.size();
  IntegerMatrix gen_matrix_can(length_can, length_can);
  // Fill in entries with -1
  for(int i = 0; i < gen_matrix_can.nrow(); i++)
  {
    for(int j = 0; j < gen_matrix_can.ncol(); j++)
    {
      if(data_can.imputed_nodes[i] == 0 && data_can.imputed_nodes[j] == 0)
      {
        gen_matrix_can(i,j) = data_can.gen_matrix(i,j);
      }
      else
      {
        gen_matrix_can(i,j) = -1;
      }
    }
  }
  return gen_matrix_can;
}

void Data::UpdateImputedNodes()
{
  bool verbose = false;
  if(verbose)
  {
    std::cout << "Attempt to update the imputed nodes" << std::endl;
  }
  genetic_contribution = 0.0;
  
  // Create a copy of the data which contains the current values
  Data data_cur((*this));
  
  ReturnNodesToImpute();
  
  
  
  if(verbose)
  {
    std::cout << "Nodes to impute cur = ";
    PrintVector(data_cur.nodes_to_impute);
    std::cout << "Nodes to impute can = ";
    PrintVector(nodes_to_impute);
    std::cout << "Times to impute cur = ";
    PrintVector(data_cur.times_to_impute);
    std::cout << "Times to impute can = ";
    PrintVector(times_to_impute);
    std::cout << "subtypes to impute cur = ";
    PrintVector(data_cur.subtypes_to_impute);
    std::cout << "subtypes to impute can = ";
    PrintVector(subtypes_to_impute);
    std::cout << "gen source cur = ";
    PrintVector(data_cur.gen_source);
    std::cout << "gen source can = ";
    PrintVector(gen_source);
    std::cout << "imputed nodes cur = ";
    PrintVector(data_cur.imputed_nodes);
    std::cout << "imputed nodes can = ";
    PrintVector(imputed_nodes);
    
    std::cout << "*** Look at the reverse process ***" << std::endl;
    
  }
  
  
  // Look at the reverse process
  
  std::vector<int> idx_to_impute_cur = WhichNotVec(0, data_cur.imputed_nodes);
  for(auto imp_idx_cur : idx_to_impute_cur)
  {
    //int current_id = data_cur.genetic_ids[imp_idx_cur];
    //double current_time = data_cur.sample_times[imp_idx_cur];
    //int current_subtype = data_cur.subtype_numbers[imp_idx_cur];
    int imputation_type = data_cur.imputed_nodes[imp_idx_cur];
    
    switch(imputation_type) 
    {
    case 1:
      // Interior node
      data_cur.CalculateContribRevInterior(imp_idx_cur, (*this), verbose);
      break;
    case 2:
      // Exterior node
      data_cur.CalculateContribRevExterior(imp_idx_cur, (*this), verbose);
      break;
    case 3:
      // Upper node
      data_cur.CalculateContribRevUpper(imp_idx_cur, (*this), verbose);
      break;
    default:
      stop("Imputation type undefined");
    }
  }
  genetic_contribution += data_cur.genetic_contribution;
  // Look at the forward process
  if(verbose)
  {
    std::cout << "*** Look at the forward process ***" << std::endl;
  }
  
  IntegerMatrix gen_mat_can = InitialiseCandidateMatrix((*this));
  
  std::vector<int> idx_to_impute_can = WhichNotVec(0, imputed_nodes);
  for(auto imp_idx_can : idx_to_impute_can)
  {
    //int current_id = genetic_ids[imp_idx_can];
    //double current_time = sample_times[imp_idx_can];
    //int current_subtype = subtype_numbers[imp_idx_can];
    int imputation_type = imputed_nodes[imp_idx_can];
    
    
    switch(imputation_type)
    {
    case 1:
      // Impute interior node
      ImputeInteriorNode(imp_idx_can, data_cur, gen_mat_can, verbose);
      break;
    case 2:
      // Impute exterior node
      ImputeExteriorNode(imp_idx_can, data_cur, gen_mat_can, verbose);
      break;
    case 3:
      // Impute upper node
      ImputeUpperNode(imp_idx_can, data_cur, gen_mat_can, verbose);
      break;
    default:
      stop("Incorrect imputation type for the forward process");
    }
  }
  
  gen_matrix = gen_mat_can;
  
  // Update coltime distances
  ImputeColtimeDistances(data_cur, verbose);
  
  return;
}

void Data::ImputeColtimeDistances(Data data_cur, bool verbose)
{
  // Look at the reverse process
  std::vector<int> import_idx_cur = WhichVec(-1, data_cur.gen_source);
  
  if(import_idx_cur.size() != data_cur.coltime_distances.size())
  {
    std::cout << "Coltime distances are not initialised correctly"
              << std::endl;
  }
  
  if(verbose)
  {
    std::cout << "*** Update the colonisation time distances ***" << std::endl;
    std::cout << "Length cur = " << import_idx_cur.size()
              << std::endl;
    std::cout << "coltime distance length = " << coltime_distances.size() << ", T_distances = ";
    PrintVector(coltime_distances);
  }
  
  for(size_t i = 0; i < import_idx_cur.size(); i++)
  {
    int idx = import_idx_cur[i];
    
    int genetic_id = data_cur.genetic_ids[idx];
    double sample_time = data_cur.sample_times[idx];
    int subtype_number = data_cur.subtype_numbers[idx];
    
    int loc_in_can = FindGeneticLoc(genetic_id, sample_time, subtype_number);
    if(verbose)
    {
      std::cout << "idx = " << idx << ", loc_in_can = "
                << loc_in_can;
    }
    
    bool calculate_contribution = true;
    if(loc_in_can != -1)
    {
      if(gen_source[loc_in_can] == -1)
      {
        calculate_contribution = false;
      }
    }
    
    if(calculate_contribution)
    {
      double infection_time = data_cur.t_e[genetic_id];
      double time_diff = sample_time - infection_time;
      double mutation_rate = parameters[6];
      double mut_prob = CalculateMutationProbability(mutation_rate, time_diff);
      
      int total_mutations = coltime_distances[i];
      
      
      
      genetic_contribution += R::dbinom(total_mutations, sequence_length, 
                                        mut_prob, 1);
      
      
      if(verbose)
      {
        std::cout << ", observed distance = " << total_mutations
                  << ", contribution = " 
                  << R::dbinom(total_mutations, 
        sequence_length, mut_prob, 1)  << std::endl;
      }
      
    }
  }
  
  
  // Forward contribution
  std::vector<int> import_idx_can = WhichVec(-1, gen_source);
  std::vector<int> coltime_distances_can(import_idx_can.size());
  
  for(size_t i = 0; i < import_idx_can.size(); i++)
  {
    int current_idx = import_idx_can[i];
    int genetic_id = genetic_ids[current_idx];
    double sample_time = sample_times[current_idx];
    int subtype_number = subtype_numbers[current_idx];
    
    
    
    int loc_in_cur = data_cur.FindGeneticLoc(genetic_id, sample_time, 
                                             subtype_number);
    
    bool simulate_distance = true;
    
    if(loc_in_cur != -1)
    {
      // The node is found in the current data set, check that it is in fact an imported node
      int current_idx_gen_source = data_cur.gen_source[loc_in_cur];
      if(current_idx_gen_source == -1)
      {
        // The node in cur data set is also imported, copy
        int cur_import_idx = WhichVec(loc_in_cur, import_idx_cur)[0];
        //if(verbose) std::cout << "current_idx_gen_source = " << current_idx_gen_source 
        //                      << ", cur import idx = " << cur_import_idx << std::endl;
        //std::cout << "coltime_distances[cur_import_idx] = " << coltime_distances[cur_import_idx] << std::endl;
        coltime_distances_can[i] = data_cur.coltime_distances[cur_import_idx];
        //std::cout << "coltime_distances_can[i] = " << coltime_distances_can[i] << std::endl;
        simulate_distance = false;
      }
    }
    
    if(simulate_distance)
    {
      //std::cout << "simulate coltime distance for i = " << i << std::endl;
      double exposure_time = t_e[genetic_id];
      double time_diff = sample_time - exposure_time;
      double mutation_rate = parameters[6];
      double mut_prob = CalculateMutationProbability(mutation_rate, time_diff);
      
      int total_mutations = R::rbinom(sequence_length, mut_prob);
      
      
      
      coltime_distances_can[i] = total_mutations;
      
      //std::cout << "forward Draw = " << draw << ", prob = " << mutation_probability << ", contribution = "
      //      << R::dbinom(draw, sequence_length, mutation_probability, 1) << std::endl;
      
      genetic_contribution -= R::dbinom(total_mutations, sequence_length, 
                                        mut_prob, 1);
      
      
      if(verbose)
      {
        //std::cout << "Simulated distance = " << draw << ", prob = " << mutation_probability << ", contribution = "
        //          << R::dbinom(draw, sequence_length, mutation_probability, 1) << std::endl;
      }
    }
  }
  coltime_distances = coltime_distances_can;
}

void Data::ImputeUpperNode(int imp_idx, Data data_cur,
                           IntegerMatrix &gen_matrix_can, bool verbose)
{
  // For now always update these
  std::vector<int> child_nodes = WhichVec(imp_idx, gen_source);
  
  if(verbose)
  {
    std::cout << " Adding upper node " << imp_idx << " from ID = "
              << genetic_ids[imp_idx] << " at t = " 
              << sample_times[imp_idx] << " with child nodes = ";
    PrintVector(child_nodes);
  }
  
  
  int gen_id_can = genetic_ids[imp_idx];
  double samp_time_can = sample_times[imp_idx];
  int subtype_can = subtype_numbers[imp_idx];
  
  int imp_idx_cur = data_cur.FindGeneticLoc(gen_id_can, samp_time_can,
                                            subtype_can);
  
  bool same_config = SameConfiguration(imp_idx_cur, imp_idx, data_cur, 
                                       (*this));
  
  if(same_config)
  {
    // Do not calculate contribution, copy from current matrix
    for(auto ch_node : child_nodes)
    {
      int child_id_can = genetic_ids[ch_node];
      double samp_time_can = sample_times[ch_node];
      int subtype_can = subtype_numbers[ch_node];
      int child_idx_cur = data_cur.FindGeneticLoc(child_id_can, samp_time_can,
                                                  subtype_can);
      int cur_dist = data_cur.gen_matrix(imp_idx_cur, child_idx_cur);
      gen_matrix_can(imp_idx, ch_node) = cur_dist;
      gen_matrix_can(ch_node, imp_idx) = cur_dist;

      if(verbose)
      {
        std::cout << "  Current child node = " << ch_node << " from ID = "
                  << genetic_ids[ch_node] << " at t = "
                  << sample_times[ch_node] 
                  << ", found in current data set - copy distance = "
                  << cur_dist
                  << std::endl;
      }
    }
  }
  else
  {
    // Node is not found in the current data set, simulate new distances
    for(auto ch_node : child_nodes)
    {
      double time_diff = sample_times[ch_node] - sample_times[imp_idx];
      double mutation_rate = parameters[6];
      double mut_prob = CalculateMutationProbability(mutation_rate, time_diff);
      int draw = R::rbinom(sequence_length, mut_prob);
      //draw = 0;
      double contribution = R::dbinom(draw, sequence_length, mut_prob, 1);
      genetic_contribution -= contribution;
      if(verbose)
      {
        std::cout << "  Current child node = " << ch_node << " from ID = "
                  << genetic_ids[ch_node] << " at t = "
                  << sample_times[ch_node] 
                  << ", mut_prob = " << mut_prob
                  << ", draw = " << draw
                  << ", contribution = " << contribution
                  << std::endl;
      }
      
      gen_matrix_can(imp_idx, ch_node) = draw;
      gen_matrix_can(ch_node, imp_idx) = draw;
    }
  }
  
  
  
  if(verbose) std::cout << std::endl;
}



void Data::ImputeExteriorNode(int imp_idx, Data data_cur,
                              IntegerMatrix &gen_matrix_can, bool verbose)
{
  int current_id = genetic_ids[imp_idx];
  double current_time = sample_times[imp_idx];
  int current_subtype = subtype_numbers[imp_idx];
  
  
  
  std::vector<int> observed_child_nodes = ReturnObservedChildren(imp_idx);
  std::vector<int> direct_obs_child_nodes;
  for(auto ch_node : observed_child_nodes)
  {
    if(gen_source[ch_node] == imp_idx)
    {
      direct_obs_child_nodes.push_back(ch_node);
    }
  }
  
  if(verbose)
  {
    std::cout << " Adding exterior node "<< imp_idx << " from ID = " << current_id
              << " at t = " << current_time << " observed child nodes = ";
    PrintVector(observed_child_nodes);
    std::cout << "  direct observed child nodes = ";
    PrintVector(direct_obs_child_nodes);
  }
  
  // Error check
  if(direct_obs_child_nodes.size() == 0)
  {
    stop("No directly observed child nodes when imputing an exterior distance");
  }
  
  
  // here parent_distance = min_i distance
  // min_i is the `parent' node
  int parent_distance = -1;
  int min_i = -1; 
  int min_distance = sequence_length;
  
  int imp_idx_cur = data_cur.FindGeneticLoc(current_id,current_time,
                                            current_subtype);
  
  bool same_config = SameConfiguration(imp_idx_cur, imp_idx,
                                       data_cur, (*this));
  

  

  if(same_config)
  {
    // Node is found in the current data set, find min_i
    for(size_t i = 0; i < direct_obs_child_nodes.size(); i++)
    {
      for(size_t j = 0; j < observed_child_nodes.size(); j++)
      {
        int node_i = direct_obs_child_nodes[i];
        int node_j = observed_child_nodes[j];
        
        if(node_i != node_j)
        {
          int current_distance = gen_matrix(node_i,node_j);
          if(current_distance < min_distance)
          {
            min_distance = current_distance;
            min_i = node_i;
          }
        }
      }
    }
    
    
    
    
    parent_distance = data_cur.gen_matrix(imp_idx_cur, min_i);
    
    if(verbose)
    {
      std::cout << "  imp_idx_can = " << imp_idx 
                << ", imp_idx_cur = " << imp_idx_cur
                << ", min_i = "
                << min_i << ", parent_distance = " << parent_distance
                << std::endl;
      std::cout << "  Same configuration - do not calculate contribution"
                << std::endl << std::endl;
    }
  }
  else
  {
    // Simulate a new distance
    // min_i will be directly connected, min_j may not be
    int min_j = -1;
    
    for(size_t i = 0; i < direct_obs_child_nodes.size(); i++)
    {
      for(size_t j = 0; j < observed_child_nodes.size(); j++)
      {
        int node_i = direct_obs_child_nodes[i];
        int node_j = observed_child_nodes[j];
        
        if(node_i != node_j)
        {
          int current_distance = gen_matrix(node_i,node_j);
          if(current_distance < min_distance)
          {
            min_distance = current_distance;
            min_i = node_i;
            min_j = node_j;
          }
        }
      }
    }
    
    // Look at the time interval 
    double upper = sample_times[min_i] - sample_times[imp_idx];
    double lower = sample_times[min_i] + sample_times[min_j] - 
      2*sample_times[imp_idx];
    
    double time_ratio = upper/lower;
    if(time_ratio < 0 || time_ratio > 1)
    {
      std::cout << "T node = " << imp_idx 
                << ", observed children = ";
      PrintVector(observed_child_nodes);
      std::cout << "min_i = " << min_i << ", min_j = " << min_j << std::endl;
      std::cout << "sample_times[min_i] = " << sample_times[min_i] << std::endl;
      std::cout << "sample_times[min_j] = " << sample_times[min_j] << std::endl;
      std::cout << "sample_times[imputed_idx_can] = " << sample_times[imp_idx] << std::endl;
      std::cout << "Time ratio = " << time_ratio << std::endl;
      stop("invalid time ratio");
    }
    
    // Simulate distance from binomial
    parent_distance = R::rbinom(min_distance, time_ratio);
    
    // Calculate contribution
    double contribution = R::dbinom(parent_distance, min_distance, 
                                    time_ratio, 1);
    genetic_contribution -= contribution;
    
    
    if(verbose)
    {
      std::cout << "  Different config, calculate contribution - "
                << " min_i = " << min_i << ", parent_distance = "
                << parent_distance << ", time_ratio = "
                << time_ratio << ", contribution = "
                << contribution << ", observed child nodes = ";
      PrintVector(observed_child_nodes);
      std::cout << std::endl;
    }
  }
  
  if(min_i == -1)
  {
    stop("Could not find min_i when imputing exterior distances");
  }
  
  // Update distances in the matrix
  
  
  gen_matrix_can(imp_idx, min_i) = parent_distance;
  gen_matrix_can(min_i, imp_idx) = parent_distance;
  
  for(size_t i = 0; i < observed_child_nodes.size(); i++)
  {
    int current_child_node = observed_child_nodes[i];
    if(current_child_node != min_i)
    {
      gen_matrix_can(imp_idx, current_child_node) = gen_matrix_can(min_i, current_child_node) -
        gen_matrix_can(imp_idx, min_i);
      gen_matrix_can(current_child_node, imp_idx) = gen_matrix_can(min_i, current_child_node) -
        gen_matrix_can(imp_idx, min_i);
      
      bool print_extra = false;
      if(print_extra)
      {
        std::cout << " current child node = " << current_child_node
                  << ", gen_matrix_can(min_i, current_child_node) = "
                  << gen_matrix_can(min_i, current_child_node)
                  << ", gen_matrix_can(imp_idx, min_i) = "
                  << gen_matrix_can(imp_idx, min_i)
                  << std::endl;
      }
    }
  }
  
  std::vector<int> import_idx_can = WhichVec(-1, gen_source);
  
  // Also update the distances between this imputed node and all other imputed nodes
  for(size_t jj = 0; jj < import_idx_can.size(); jj++)
  {
    int cur_imp_node = import_idx_can[jj];
    if(cur_imp_node != imp_idx)
    {
      gen_matrix_can(imp_idx, cur_imp_node) = gen_matrix_can(min_i, cur_imp_node) -
        gen_matrix_can(imp_idx, min_i);
      gen_matrix_can(cur_imp_node, imp_idx) = gen_matrix_can(min_i, cur_imp_node) -
        gen_matrix_can(imp_idx, min_i);
    }
  }
}

// Returns the minimum possible distance as a consequence from the triangle inequality
int Data::ReturnMinDistanceTriangle(int imp_idx, const IntegerMatrix gen_matrix_can)
{
  bool verbose = false;
  //if(imp_idx == 30) verbose = true;
  int min_distance = sequence_length;
  int parent_node = gen_source[imp_idx];
  if(imputed_nodes[imp_idx]==0) stop("the node is not imputed");
  
  std::vector<int> set1;
  std::vector<int> set2;
  
  std::vector<int> direct_child_nodes = WhichVec(imp_idx, gen_source);
  if(direct_child_nodes.size() != 2) stop("incorrect number of direct child nodes");
  
  int node1 = direct_child_nodes[0];
  int node2 = direct_child_nodes[1];
  if(imputed_nodes[node1] == 0)
  {
    set1.push_back(node1);
  }
  else
  {
    // Direct child is imputed, look at obs child nodes instead
    std::vector<int> obs_ch_nodes = ReturnObservedChildren(node1);
    for(auto obs_ch_node : obs_ch_nodes)
    {
      set1.push_back(obs_ch_node);
    }
  }
  
  if(imputed_nodes[node2] == 0)
  {
    set2.push_back(node2);
  }
  else
  {
    // Direct child is imputed, look at obs child nodes instead
    std::vector<int> obs_ch_nodes = ReturnObservedChildren(node2);
    for(auto obs_ch_node : obs_ch_nodes)
    {
      set2.push_back(obs_ch_node);
    }
  }
  
  if(verbose)
  {
    std::cout << "  Perform triangle inequality checks for imputed node " << imp_idx
              << ", set1 = ";
    for(auto i : set1)
    {
      std::cout << i << " ";
    }
    std::cout << "and set2 = ";
    for(auto i : set2)
    {
      std::cout << i << " ";
    }
    std::cout << std::endl;
  }
  
  for(auto i : set1)
  {
    for(auto j : set2)
    {
      int dist = (gen_matrix_can(parent_node,i) + gen_matrix_can(parent_node,j) - gen_matrix_can(i,j))/2;
      if(verbose)
      {
        std::cout << " i = " << i << ", j = " << j
                  << ", dist = " << dist << std::endl;
        std::cout << "  gen_matrix_can(parent_node,i) = "
                  << gen_matrix_can(parent_node,i)
                  << ", gen_matrix_can(parent_node,j) = "
                  << gen_matrix_can(parent_node,j)
                  << ", gen_matrix_can(i,j) = "
                  << gen_matrix_can(i,j)
                  << std::endl;
      }
      if(dist < min_distance)
      {
        min_distance = dist;
      }
    }
  }
  
  
  return min_distance;
}



void Data::ImputeInteriorNode(int imp_idx, Data data_cur, 
                              IntegerMatrix &gen_matrix_can, bool verbose)
{
  // First check if the node appears in the current data set
  
  int current_id = genetic_ids[imp_idx];
  double current_time = sample_times[imp_idx];
  int current_subtype = subtype_numbers[imp_idx];
  int parent_node_can = gen_source[imp_idx];
  std::vector<int> observed_child_nodes = ReturnObservedChildren(imp_idx);
  
  if(verbose)
  {
    std::cout << " Adding interior node " << imp_idx << " from ID = "
              << genetic_ids[imp_idx] << " at t = " 
              << sample_times[imp_idx] << ", with parent node = "
              << parent_node_can << " and observed child nodes = ";
    PrintVector(observed_child_nodes);
  }
  
  int imp_idx_cur = data_cur.FindGeneticLoc(current_id,current_time,
                                            current_subtype);
  
  
  int parent_distance = -1;
  
  // For now do not bother partial updating
  //bool same_config = DetermineIfImputedConfigSame(imp_idx_cur, imp_idx,
  //                                                data_cur, (*this));
  bool same_config = SameConfiguration(imp_idx_cur, imp_idx,
                                       data_cur, (*this));
  if(same_config)
  {
    // Copy distances from the current matrix
    int parent_cur = data_cur.gen_source[imp_idx_cur];
    parent_distance = data_cur.gen_matrix(imp_idx_cur, parent_cur);
    if(verbose)
    {
      std::cout << "  Same config - parent_cur = " << parent_cur
                << ", parent_distance = " << parent_distance << std::endl
                << std::endl;
    }
  }
  else
  {
    // Simulate distances from the model
    int min_distance = gen_matrix_can(parent_node_can, observed_child_nodes[0]);
    double min_time = sample_times[observed_child_nodes[0]];
    for(size_t i = 1; i < observed_child_nodes.size(); i++)
    {
      int current_child_node = observed_child_nodes[i];
      
      int current_distance = gen_matrix_can(parent_node_can, current_child_node);
      if(current_distance <= min_distance)
      {
        min_distance = current_distance;
        double current_time = sample_times[current_child_node];
        min_time = std::min(min_time, current_time);
      }
    }
    
    // Now perform the triangle inequality checks
    
    int triangle_min = ReturnMinDistanceTriangle(imp_idx, gen_matrix_can);
    //std::cout << "min distance = " << min_distance << ", triangle min = "
    //          << triangle_min << std::endl;
    min_distance = std::min(min_distance, triangle_min);
    
    if(min_distance < 0)
    {
      // Violation of the triangle inequality
      min_distance = 0;
      genetic_contribution -= 9e16;
    }
    
    
    
    // Calculate time ratios
    double upper = sample_times[imp_idx] - sample_times[parent_node_can];
    double lower = min_time - sample_times[parent_node_can];
    double time_ratio = upper/lower;
    
    if(time_ratio < 0 || time_ratio > 1)
    {
      //IntegerVector gen_source_cur = ReturnGenSourceVector_Rcpp(data_cur);
      //IntegerVector gen_source_can = ReturnGenSourceVector_Rcpp((*this));
      //Rcout << gen_source_cur << std::endl;
      //Rcout << gen_source_can << std::endl;
      std::cout << std::endl << "Error" << std::endl;
      std::cout << "parent node = " << parent_node_can << ", observed children = ";
      PrintVector(observed_child_nodes);
      std::cout << "sample_times[imputed_idx_can] = " << sample_times[imp_idx] << std::endl;
      std::cout << "sample_times[parent_node] = " << sample_times[parent_node_can] << std::endl;
      std::cout << "min time = " << min_time << std::endl;
      std::cout << "Time ratio = " << time_ratio << std::endl;
      
      std::vector<int> ever_infected;
      for(int i = 0; i < N; i++)
      {
        if(t_e[i] != -1) ever_infected.push_back(i);
      }
      
      for(auto i : ever_infected) {
        std::cout << "i = " << i  << ", t_e = " 
                  << t_e[i]
                  << ", t_r = " << t_r[i] << ", source = " << source[i]
                  << std::endl;
      }
      std::cout << "genetic_ids = ";
      PrintVector(genetic_ids);
      std::cout << "sample_times = ";
      PrintVector(sample_times);
      
      stop("invalid time ratio for T");
    }
    
    // Simulate a distance between the imputed node and the parent
    parent_distance = R::rbinom(min_distance, time_ratio);
    
    // Update the genetic contribution
    double contribution = R::dbinom(parent_distance, min_distance, 
                                    time_ratio, 1);
    genetic_contribution -= contribution;
    
    if(verbose)
    {
      std::cout << "  Different config" << ", min_distance = "
                << min_distance << ", parent_distance = "
                << parent_distance << ", time_ratio = "
                << time_ratio << ", contribution = "
                << contribution
                << std::endl << std::endl;
    }
  }
  
  // Fill in matrix entries
  gen_matrix_can(parent_node_can, imp_idx) = parent_distance;
  gen_matrix_can(imp_idx, parent_node_can) = parent_distance;
  for(size_t i = 0; i < observed_child_nodes.size(); i++)
  {
    int current_child_node = observed_child_nodes[i];
    gen_matrix_can(imp_idx, current_child_node) = gen_matrix_can(parent_node_can, current_child_node) -
      gen_matrix_can(parent_node_can, imp_idx);
    gen_matrix_can(current_child_node, imp_idx) = gen_matrix_can(parent_node_can, current_child_node) -
      gen_matrix_can(parent_node_can, imp_idx);
  }
  
}


std::vector<int> Data::ReturnIDsBetweenNodes(int node, 
                                             std::vector<int> child_nodes)
{
  std::vector<int> ids_to_return;
  int node_id = genetic_ids[node];
  
  for(size_t i = 0; i < child_nodes.size(); i++)
  {
    int current_node = child_nodes[i];
    int current_node_id = genetic_ids[current_node];
    
    while(current_node_id != node_id && current_node_id != -1)
    {
      ids_to_return.push_back(current_node_id);
      current_node_id = source[current_node_id];
    }
  }
  return sort_unique(ids_to_return);
}


void Data::ReturnNodesToImpute(bool verbose)
{
  verbose = false;
  if(verbose)
  {
    std::cout << std::endl << "*** Attempt to return nodes to impute ***" 
              << std::endl;
  }
  //PrintMemberVariables();
  
  // Vectors to store the nodes to impute
  std::vector<int> genetic_ids_to_impute;
  std::vector<double> sample_times_to_impute;
  std::vector<int> subtype_numbers_to_impute;
  
  // Data temp is a temporary store which is used to calculate the genetic source
  // vector using only observed data for determining which interior nodes need
  // to be imputed
  Data data_temp((*this));
  std::vector<int> observed_index = WhichVec(0, imputed_nodes);
  data_temp.genetic_ids = subset(observed_index, data_temp.genetic_ids);
  data_temp.sample_times = subset(observed_index, data_temp.sample_times);
  data_temp.subtype_numbers = subset(observed_index, data_temp.subtype_numbers);
  
  data_temp.CalculateGenSourceVector();
  
  if(verbose)
  {
    std::cout << "genetic ids = ";
    PrintVector(genetic_ids);
    std::cout << "sample times = ";
    PrintVector(sample_times);
    std::cout << " Observed gen source = ";
    PrintVector(data_temp.gen_source);
    //std::cout << "t_c" << std::endl;
    //PrintVector(t_c);
    std::cout << "source = ";
    PrintVector(source);
  }
  
  
  // Calculate interior nodes to impute
  for(size_t i = 0; i < data_temp.gen_source.size(); i++)
  {
    std::vector<int> child_nodes = WhichVec((int)i, data_temp.gen_source);
    if(child_nodes.size() > 1)
    {
      // Current node has two or more children, therefore check if 
      // imputation is necessary
      std::vector<int> colonisation_events = 
        data_temp.ReturnIDsBetweenNodes(i, child_nodes);
      if(verbose)
      {
        std::cout << " for node " << i << ", colonisation events = ";
        PrintVector(colonisation_events);
      }
      
      // This temp data is created for calculating the gen source of the small
      // sub trees
      Data data_temp2((*this));
      int current_variant = subtype_numbers[i];
      int num_observed_nodes = 1 + child_nodes.size();
      std::vector<int> temp_genetic_ids(num_observed_nodes);
      std::vector<double> temp_sample_times(num_observed_nodes);
      std::vector<int> temp_subtype_numbers(num_observed_nodes);
      
      temp_genetic_ids[0] = genetic_ids[i];
      temp_sample_times[0] = sample_times[i];
      temp_subtype_numbers[0] = current_variant;
      
      for(int j = 1; j < num_observed_nodes; j++)
      {
        temp_genetic_ids[j] = data_temp.genetic_ids[child_nodes[j-1]];
        temp_sample_times[j] = data_temp.sample_times[child_nodes[j-1]];
        temp_subtype_numbers[j] = data_temp.subtype_numbers[child_nodes[j-1]];
      }
      
      for(size_t j = 0; j < colonisation_events.size(); j++)
      {
        int current_id = colonisation_events[j];
        int current_source = source[current_id];
        if(current_source != -1)
        {
          // Check if the node is already in the temp IDs
          int loc = -1;
          for(size_t ii = 0; ii < temp_genetic_ids.size(); ii++)
          {
            if(current_source == temp_genetic_ids[ii] &&
               t_e[current_id] == temp_sample_times[ii] &&
               current_variant == temp_subtype_numbers[ii])
            {
              loc = ii;
              break;
            }
            else if(current_id == temp_genetic_ids[ii] &&
                    t_e[current_id] == temp_sample_times[ii] &&
                    current_variant == temp_subtype_numbers[ii])
            {
              loc = ii;
              break;
            }
          }
          
          if(loc == -1)
          {
            // if loc == -1, then the node is not in the temp ids
            temp_genetic_ids.push_back(current_source);
            temp_sample_times.push_back(t_e[current_id]);
            temp_subtype_numbers.push_back(current_variant);
          }
        }
      }
      
      data_temp2.genetic_ids = temp_genetic_ids;
      data_temp2.sample_times = temp_sample_times;
      data_temp2.subtype_numbers = temp_subtype_numbers;
      
      data_temp2.CalculateGenSourceVector();
      
      if(verbose)
      {
        std::cout << std::endl << "index i = " << i << " has " 
                  << child_nodes.size() << " child nodes = ";
        PrintVector(child_nodes);
        std::cout << "IDs between nodes = ";
        PrintVector(colonisation_events);
        std::cout << "genetic ids = ";
        PrintVector(data_temp2.genetic_ids);
        std::cout << "sample times = ";
        PrintVector(data_temp2.sample_times);
        std::cout << "subtype numbers = ";
        PrintVector(data_temp2.subtype_numbers);
        std::cout << "gen source = ";
        PrintVector(data_temp2.gen_source);
      }
      
      for(size_t ii = num_observed_nodes; ii < data_temp2.gen_source.size(); ii++)
      {
        // check there is not a zero time difference between 
        // the current node and parent
        bool observed_time_diff_zero = false;
        int parent_node = data_temp2.gen_source[ii];
        double time_diff = temp_sample_times[ii] - 
          temp_sample_times[parent_node];
        if(verbose)
        {
          std::cout << " node = " << ii << ", parent = " << parent_node 
                    << ", time diff = " << time_diff << std::endl;
        }
        if(time_diff != 0)
        {
          // there is a non zero time between the node and the parent
          std::vector<int> current_emputed_child_nodes = 
            WhichVec((int)ii, data_temp2.gen_source);
          
          // check the number of child nodes is greater than one
          if(current_emputed_child_nodes.size() > 1)
          {
            // There is more than one child node, check they have positive time
            for(size_t jj = 0; jj < current_emputed_child_nodes.size(); jj++)
            {
              int current_child_node = current_emputed_child_nodes[jj];
              time_diff = temp_sample_times[current_child_node] - 
                temp_sample_times[ii];
              if(time_diff == 0)
              {
                observed_time_diff_zero = true;
                break;
              }
            }
            
            if(!observed_time_diff_zero)
            {
              if(verbose)
              {
                std::cout << " Impute ID = " << temp_genetic_ids[ii] 
                          << ", at time t = " << temp_sample_times[ii]
                          << ", with subtype = " << temp_subtype_numbers[ii] 
                          << std::endl;
                std::cout << " current_emputed_child_nodes = ";
                PrintVector(current_emputed_child_nodes);
              }
              genetic_ids_to_impute.push_back(temp_genetic_ids[ii]);
              sample_times_to_impute.push_back(temp_sample_times[ii]);
              subtype_numbers_to_impute.push_back(temp_subtype_numbers[ii]);
            }
          }
        }
      }
      
    }
  }
  
  if(verbose)
  {
    std::cout << "Interior nodes calculated" << std::endl;
    std::cout << " ids = ";
    PrintVector(genetic_ids_to_impute);
    std::cout << " times = ";
    PrintVector(sample_times_to_impute);
    std::cout << " subtypes = ";
    PrintVector(subtype_numbers_to_impute);
  }
  
  
  // Calculate exterior nodes to impute
  
  // This is a copy of the current state of the chain to calculate exterior nodes
  //Data data_temp2((*this));
  
  //std::cout << "*this gen source = ";
  //PrintVector(gen_source);
  
  
  // This will return the extiror nodes and store them in
  // nodes to impute
  // times to impute
  // subtypes to impute
  //std::cout << "data_temp gen source = ";
  //PrintVector(data_temp.gen_source);
  data_temp.ReturnExteriorNodes();
  
  if(verbose)
  {
    std::cout << " Exterior nodes calculated" << std::endl;
    std::cout << " ids = ";
    PrintVector(data_temp.nodes_to_impute);
    std::cout << " times = ";
    PrintVector(data_temp.times_to_impute);
    std::cout << " subtypes = ";
    PrintVector(data_temp.subtypes_to_impute);
  }
  
  for(size_t i = 0; i < data_temp.nodes_to_impute.size(); i++)
  {
    genetic_ids_to_impute.push_back(data_temp.nodes_to_impute[i]);
    sample_times_to_impute.push_back(data_temp.times_to_impute[i]);
    subtype_numbers_to_impute.push_back(data_temp.subtypes_to_impute[i]);
  }
  // Calculate the genetic source vector with these new nodes
  
  
  // Update genetic ids, times, subtypes and imputed nodes
  std::vector<int> genetic_ids_can;
  std::vector<double> sample_times_can;
  std::vector<int> subtype_numbers_can;
  std::vector<int> imputed_nodes_can;
  
  
  
  int observed_length = 0;
  for(size_t i = 0; i < imputed_nodes.size(); i++)
  {
    if(imputed_nodes[i] == 0)
    {
      genetic_ids_can.push_back(genetic_ids[i]);
      sample_times_can.push_back(sample_times[i]);
      subtype_numbers_can.push_back(subtype_numbers[i]);
      imputed_nodes_can.push_back(0);
      observed_length++;
    }
  }
  
  
  
  for(size_t i = 0; i < genetic_ids_to_impute.size(); i++)
  {
    genetic_ids_can.push_back(genetic_ids_to_impute[i]);
    sample_times_can.push_back(sample_times_to_impute[i]);
    subtype_numbers_can.push_back(subtype_numbers_to_impute[i]);
    imputed_nodes_can.push_back(1);
  }
  
  genetic_ids = genetic_ids_can;
  sample_times = sample_times_can;
  subtype_numbers = subtype_numbers_can;
  imputed_nodes = imputed_nodes_can;
  
  // Calculate new gen source vector
  CalculateGenSourceVector();
  
  
  
  // Determine order of imputation
  std::vector<int> order_to_impute = DetermineOrderOfImputation();
  std::vector<int> imputation_type_can = DetermineImputationGroups(order_to_impute);
  
  
  std::vector<int> nodes_to_impute_;
  std::vector<double> times_to_impute_;
  std::vector<int> subtypes_to_impute_;
  for(size_t i = 0; i < order_to_impute.size(); i++)
  {
    int idx = order_to_impute[i];
    genetic_ids[observed_length + i] = genetic_ids_can[idx];
    sample_times[observed_length + i] = sample_times_can[idx];
    subtype_numbers[observed_length + i] = subtype_numbers_can[idx];
    imputed_nodes[observed_length + i] = imputation_type_can[i];
    
    nodes_to_impute_.push_back(genetic_ids_can[idx]);
    times_to_impute_.push_back(sample_times_can[idx]);
    subtypes_to_impute_.push_back(subtype_numbers_can[idx]);
  }
  
  nodes_to_impute = nodes_to_impute_;
  times_to_impute = times_to_impute_;
  subtypes_to_impute = subtypes_to_impute_;
  
  
  CalculateGenSourceVector();
  
}

int NonUnifSampling(std::vector<int> x, std::vector<double> weights)
{
  double U = R::runif(0.0,1.0);
  double cur_val = 0.0;
  int num_samples = x.size();
  for(int i = 0; i < num_samples; i++)
  {
    cur_val += weights[i];
    if(U < cur_val)
    {
      return x[i];
    }
  }
  return -1;
}









// Return the infectors associated with the exposure time
// of the target
std::vector<int> Data::ReturnPossibleInfectors(int target)
{
  std::vector<int> possible_infectors;
  double infection_time = t_e[target];
  for(int i = 0; i < N; i++)
  {
    if(i != target)
    {
      if(t_i[i] < infection_time && t_r[i] > infection_time)
      {
        possible_infectors.push_back(i);
      }
    }
  }
  return possible_infectors;
}

void Data::UpdateParameterMetropolisHastingsUniform(int idx, 
                                                    double proposal_variance,
                                                    double prior_lower,
                                                    double prior_upper,
                                                    int &nacc, 
                                                    bool verbose)
{
  Data data_can(*this);
  data_can.parameters[idx] = R::rnorm(parameters[idx],
                                      proposal_variance);
  
  if(data_can.parameters[idx] > prior_lower &&
     data_can.parameters[idx] < prior_upper)
  {
    double logpi_cur = loglik;
    data_can.CalculateLoglik();
    double logpi_can = data_can.loglik;
    double U = R::runif(0.0,1.0);
    if(verbose)
    {
      std::cout << "Updating parameter " << idx << ", logpi can = " << logpi_can
                << ", logpi cur = " << logpi_cur << ", loglik can = " 
                << data_can.loglik << ", loglik cur = " << loglik;
    }
    if(log(U) < logpi_can - logpi_cur)
    {
      loglik = data_can.loglik;
      parameters = data_can.parameters;
      nacc++;
      if(verbose)
      {
        std::cout << " - accepted" << std::endl;
        std::cout << "  Parameter cur = " << parameters[idx] << ", parameter can = "
                  << data_can.parameters[idx] << std::endl;
      }
    }
    else
    {
      if(verbose)
      {
        std::cout << " - rejected" << std::endl;
        std::cout << "  Parameter cur = " << parameters[idx] << ", parameter can = "
                  << data_can.parameters[idx] << std::endl;
      }
    }
  }
}

void Data::UpdateInfectionTime(int &nacc, bool verbose)
{
  verbose = false;
  double log_prop_ratio = 0.0;
  Data data_can((*this));
  int target = SampleVector(ever_infected);
  int target_source = source[target];
  if(target_source != -1)
  {
    // Calculate the last possible time the target may be exposed
    
    double max_time = t_r[target];
    
    std::vector<int> target_infections = WhichVec(target, source);
    for(auto i : target_infections)
    {
      max_time = std::min(max_time, t_e[i]);
    }
    
    std::vector<int> gen_idx = WhichVec(target, genetic_ids);
    for(auto i : gen_idx)
    {
      max_time = std::min(max_time, sample_times[i]);
    }
    
    double t_i_can = R::runif(t_e[target], max_time);
    data_can.t_i[target] = t_i_can;
    
    
    // do not need to bother updating genetic distances, the genetic network is unchanged
    bool impute_genetic_distances = false;
    if(impute_genetic_distances)
    {
      data_can.UpdateImputedNodes();
      log_prop_ratio += data_can.genetic_contribution;
    }
    
    
    
    data_can.CalculateLoglik();
    if(verbose)
    {
      std::cout << "Update an infection time for i = " << target
                << " from t = " << t_i[target] << " to t = "
                << data_can.t_i[target] << " from source = "
                << source[target] << " to source = "
                << data_can.source[target] << std::endl;
    }
    
    double U = R::runif(0,1);
    if(log(U) < data_can.loglik - loglik + log_prop_ratio)
    {
      if(verbose)
      {
        std::cout << " Loglik cur = " << loglik << ", loglik can = "
                  << data_can.loglik 
                  << ", log prop ratio = " << log_prop_ratio
                  << " - accepted." << std::endl;
      }
      t_i = data_can.t_i;
      source = data_can.source;
      loglik = data_can.loglik;
      exposure_likelihood = data_can.exposure_likelihood;
      transmission_likelihood = data_can.transmission_likelihood;
      removal_likelihood = data_can.removal_likelihood;
      genetic_likelihood = data_can.genetic_likelihood;
      nacc++;
      
      genetic_ids = data_can.genetic_ids;
      sample_times = data_can.sample_times;
      subtype_numbers = data_can.subtype_numbers;
      gen_source = data_can.gen_source;
      gen_matrix = data_can.gen_matrix;
      
      
      coltime_distances = data_can.coltime_distances;
      
      imputed_nodes = data_can.imputed_nodes;
      nodes_to_impute = data_can.nodes_to_impute;
      times_to_impute = data_can.times_to_impute;
      subtypes_to_impute = data_can.subtypes_to_impute;
    }
    else
    {
      if(verbose)
      {
        std::cout << " Loglik cur = " << loglik << ", loglik can = "
                  << data_can.loglik 
                  << ", log prop ratio = " << log_prop_ratio
                  << " - rejected." << std::endl;
      }
    }
  }
}

void Data::UpdateExposureTime(int &nacc, bool verbose)
{
  verbose = false;
  double log_prop_ratio = 0.0;
  Data data_can((*this));
  int target = SampleVector(ever_infected);
  int target_source = source[target];
  if(target_source != -1)
  {
    // Calculate the last possible time the target may be exposed
    
    double alpha = parameters[0];
    double zeta = parameters[5];
    double latent_period_can = R::rgamma(zeta, 1/alpha);
    double t_e_can = t_i[target] - latent_period_can;
    
    data_can.t_e[target] = t_e_can;
    
    // Determine a new source
    std::vector<int> possible_infectors_cur = ReturnPossibleInfectors(target);
    std::vector<int> possible_infectors_can = 
      data_can.ReturnPossibleInfectors(target);
    
    // Exit early if there are no infectors at the proposed infection time
    if(possible_infectors_can.size()==0) return;
    
    // Uniformly sample from the available sources
    int source_can = SampleVector(possible_infectors_can);
    data_can.source[target] = source_can;
    double latent_period_cur = t_i[target]-t_e[target];
    double upper = R::dgamma(latent_period_cur,zeta,1/alpha,1) + 
      log(possible_infectors_can.size());
    double lower = R::dgamma(latent_period_can,zeta,1/alpha,1) + 
      log(possible_infectors_cur.size());
    log_prop_ratio = upper - lower;
    
    /*
     std::cout << "R::dgamma(latent_period_cur,zeta,1/alpha,1) = "
               << R::dgamma(latent_period_cur,zeta,1/alpha,1)
               << ", log(possible_infectors_can.size()) = "
               << log(possible_infectors_can.size())
               << ", R::dgamma(latent_period_can,zeta,1/alpha,1) = "
               << R::dgamma(latent_period_can,zeta,1/alpha,1)
               << ", log(possible_infectors_cur.size()) = "
               << log(possible_infectors_cur.size())
               << std::endl;
     */
    bool impute_genetic_distances = true;
    if(impute_genetic_distances)
    {
      data_can.UpdateImputedNodes();
      log_prop_ratio += data_can.genetic_contribution;
    }
    
    
    
    data_can.CalculateLoglik();
    if(verbose)
    {
      std::cout << "Update an exposure time for i = " << target
                << " from t = " << t_e[target] << " to t = "
                << data_can.t_e[target] << " from source = "
                << source[target] << " to source = "
                << data_can.source[target] << std::endl;
    }
    
    double U = R::runif(0,1);
    if(log(U) < data_can.loglik - loglik + log_prop_ratio)
    {
      if(verbose)
      {
        std::cout << " Loglik cur = " << loglik << ", loglik can = "
                  << data_can.loglik 
                  << ", log prop ratio = " << log_prop_ratio
                  << " - accepted." << std::endl;
      }
      t_e = data_can.t_e;
      source = data_can.source;
      loglik = data_can.loglik;
      exposure_likelihood = data_can.exposure_likelihood;
      transmission_likelihood = data_can.transmission_likelihood;
      removal_likelihood = data_can.removal_likelihood;
      genetic_likelihood = data_can.genetic_likelihood;
      nacc++;
      
      genetic_ids = data_can.genetic_ids;
      sample_times = data_can.sample_times;
      subtype_numbers = data_can.subtype_numbers;
      gen_source = data_can.gen_source;
      gen_matrix = data_can.gen_matrix;
      
      
      coltime_distances = data_can.coltime_distances;
      
      imputed_nodes = data_can.imputed_nodes;
      nodes_to_impute = data_can.nodes_to_impute;
      times_to_impute = data_can.times_to_impute;
      subtypes_to_impute = data_can.subtypes_to_impute;
    }
    else
    {
      if(verbose)
      {
        std::cout << " Loglik cur = " << loglik << ", loglik can = "
                  << data_can.loglik 
                  << ", log prop ratio = " << log_prop_ratio
                  << " - rejected." << std::endl;
      }
    }
  }
}


// [[Rcpp::export]]
void MCMC_SEIR_JC(List MCMC_options,
                  int sequence_length,
                  NumericVector t_e,
                  NumericVector t_i, 
                  NumericVector t_r,
                  IntegerVector source,
                  NumericVector x,
                  NumericVector y,
                  IntegerVector genetic_ids, 
                  NumericVector sample_times, 
                  IntegerVector subtype_numbers, 
                  IntegerMatrix gen_matrix)
{
  

  // Load data from the options
  NumericVector parameters = MCMC_options["initial_chain_state"];

  int num_augmented_updates = MCMC_options["num_aug_updates"];
  int max_iterations = MCMC_options["iterations"];
  std::string output_file = MCMC_options["output_file"];
  List debug_flags = MCMC_options["debug_flags"];
  List proposal_variance = MCMC_options["proposal_variance"];
  List prior_parameters = MCMC_options["prior_parameters"];

  // Instantiate the data
  Data data(sequence_length,
            t_e,
            t_i,
            t_r,
            source,
            x,
            y,
            parameters,
            genetic_ids,
            sample_times,
            subtype_numbers,
            gen_matrix);
  

  
  
  // Acceptance counters
  int nacc_kappa = 0;
  int nacc_kappa_prop = 0;
  int nacc_delta = 0;
  int nacc_delta_prop = 0;
  int nacc_zeta = 0;
  int nacc_zeta_prop = 0;
  int nacc_lambda = 0;
  int nacc_lambda_prop = 0;
  int nacc_exp = 0;
  int nacc_exp_prop = 0;
  int nacc_inf = 0;
  int nacc_inf_prop = 0;
  int nacc_dist = 0;
  int nacc_dist_prop = 0;
  
  // Debug file
  std::string debug_file = "debug_file.dat";
  remove(debug_file.c_str());
  std::ofstream myfile2;
  myfile2.open(debug_file.c_str());
  assert(myfile2.is_open());
  
  // Write the current state of the chain
  myfile2 << 1 << std::endl;
  data.WriteConfigurationToFile(myfile2);
  
  
  
  // Output file
  remove(output_file.c_str());
  std::ofstream myfile; // Define output stream
  myfile.open(output_file.c_str()); // Open file
  assert(myfile.is_open());
  
  // Write current state of the chain to file
  data.WriteOutputToFile(myfile);
  
  
  // Begin MCMC
  for(int iter = 1; iter < max_iterations; iter++)
  {
    // Update alpha by gibbs
    int alpha_flag = debug_flags["alpha"];
    if(alpha_flag == 0)
    {
      double prior_shape = prior_parameters["alpha_shape"];
      double prior_rate = prior_parameters["alpha_rate"];
      double zeta = data.parameters[5];
      double latent_period_sum = 0.0;
      double n_I = 0.0;
      for(auto i : data.ever_infected)
      {
        latent_period_sum += (data.t_i[i] - data.t_e[i]);
        n_I += 1;
      }
      data.parameters[0] = R::rgamma(zeta*n_I + prior_shape,
                                     1/(prior_rate + latent_period_sum));
    }
    
    // Update beta by gibbs
    int beta_flag = debug_flags["beta"];
    if(beta_flag==0)
    {
      double double_sum = data.CalculateDoubleSum();
      double prior_shape = prior_parameters["beta_shape"];
      double prior_rate = prior_parameters["beta_rate"];
      data.parameters[1] = R::rgamma(data.n_I + prior_shape,
                                     1/(prior_rate + double_sum));
    }
    
    // Update gamma by gibbs
    int gamma_flag = debug_flags["gamma"];
    if(gamma_flag==0)
    {
      double infectious_period_sum = 0.0;
      double delta = parameters[4];
      double prior_shape = prior_parameters["gamma_shape"];
      double prior_rate = prior_parameters["gamma_rate"];
      double n_R = 0.0;
      for(auto i : data.ever_infected)
      {
        infectious_period_sum += (data.t_r[i]-data.t_i[i]);
        n_R += 1;
      }
      data.parameters[2] = R::rgamma(delta*n_R + prior_shape,
                                     1/(prior_rate + infectious_period_sum));
    }
    
    // Calculate loglik of new sampled parameters
    data.CalculateLoglik();
    
    // Update kappa by MH
    int kappa_flag = debug_flags["kappa"];
    if(kappa_flag==0)
    {
      double prop_var = proposal_variance["kappa"];
      double kappa_lower = prior_parameters["kappa_lower"];
      double kappa_higher = prior_parameters["kappa_higher"];
      nacc_kappa_prop++;
      data.UpdateParameterMetropolisHastingsUniform(3, 
                                                    prop_var,
                                                    kappa_lower,
                                                    kappa_higher,
                                                    nacc_kappa);
    }
    
    // Update delta by MH
    int delta_flag = debug_flags["delta"];
    if(delta_flag==0)
    {
      double prop_var = proposal_variance["delta"];
      double delta_lower = prior_parameters["delta_lower"];
      double delta_higher = prior_parameters["delta_higher"];
      nacc_delta_prop++;
      data.UpdateParameterMetropolisHastingsUniform(4, 
                                                    prop_var,
                                                    delta_lower,
                                                    delta_higher,
                                                    nacc_delta);
    }
    
    // Update zeta by MH
    int zeta_flag = debug_flags["zeta"];
    if(zeta_flag == 0)
    {
      double prop_var = proposal_variance["zeta"];
      double zeta_lower = prior_parameters["zeta_lower"];
      double zeta_higher = prior_parameters["zeta_higher"];
      nacc_zeta_prop++;
      data.UpdateParameterMetropolisHastingsUniform(5, 
                                                    prop_var,
                                                    zeta_lower,
                                                    zeta_higher,
                                                    nacc_zeta);
    }
    
    // Update lambda by MH
    int lambda_flag = debug_flags["lambda"];
    if(lambda_flag == 0)
    {
      double prop_var = proposal_variance["lambda"];
      double lambda_lower = prior_parameters["lambda_lower"];
      double lambda_higher = prior_parameters["lambda_higher"];
      nacc_lambda_prop++;
      data.UpdateParameterMetropolisHastingsUniform(6, 
                                                    prop_var,
                                                    lambda_lower,
                                                    lambda_higher,
                                                    nacc_lambda);
    }
    
    
    // Update the augmented data
    int aug_data_flag = debug_flags["aug_data"];
    if(aug_data_flag==0)
    {
      for(int i = 0; i < num_augmented_updates; i++)
      {
        int move = floor(R::runif(0,2));
        //int move = 0;
        if(move==0)
        {
          nacc_exp_prop++;
          data.UpdateExposureTime(nacc_exp, false);
        }
        else if(move==1)
        {
          nacc_inf_prop++;
          data.UpdateInfectionTime(nacc_inf, false);
        }
      }
      
      
      
    }
    
    
    
    
    
    
    //TriangleInequalityCheck(data);
    data.WriteOutputToFile(myfile);
    myfile2 << iter + 1<< std::endl;
    data.WriteConfigurationToFile(myfile2);
    
    
    if(data.loglik < -1e10) stop("-inf loglik");
    
    if(iter%100==0)
    {
      std::cout << "Iteration " << iter << " completed." << std::endl;
      checkUserInterrupt();
    }
  }
  
  
  double kappa_prob = (double)nacc_kappa/(double)nacc_kappa_prop;
  double delta_prob = (double)nacc_delta/(double)nacc_delta_prop;
  double zeta_prob = (double)nacc_zeta/(double)nacc_zeta_prop;
  double lambda_prob = (double)nacc_lambda/(double)nacc_lambda_prop;
  double exp_prob = (double)nacc_exp/(double)nacc_exp_prop;
  double inf_prob = (double)nacc_inf/(double)nacc_inf_prop;
  double dist_prob = (double)nacc_dist/(double)nacc_dist_prop;
  
  std::cout << "Acceptance probabilities: kappa = " << kappa_prob
            << ", delta = " << delta_prob
            << ", zeta = " << zeta_prob
            << ", lambda = " << lambda_prob
            << ", exp_time = " << exp_prob
            << ", inf_time = " << inf_prob
            << ", dist_prop = " << dist_prob
            << std::endl;
  
  myfile.close();
  myfile2.close();
}



