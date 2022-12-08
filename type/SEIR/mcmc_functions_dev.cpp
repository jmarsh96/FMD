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
  
  
  
  std::vector<int> coltime_distances;
  

  

  
  
  // Likelihood
  double loglik;
  double removal_likelihood;
  double transmission_likelihood;
  double exposure_likelihood;
  
  std::vector<double> parameters;
  
  // Constructor
  Data(NumericVector t_e_,
       NumericVector t_i_,
       NumericVector t_r_,
       IntegerVector source_,
       NumericVector x_,
       NumericVector y_,
       NumericVector parameters_);
  
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
  void CalculateContribRevUpper(int imp_idx, bool verbose = false);
  void CalculateContribRevExterior(int imp_idx, bool verbose = false);
  void CalculateContribRevInterior(int imp_idx, bool verbose = false);
  void ImputeUpperNode(int imp_idx, Data data_cur,
                       IntegerMatrix &gen_matrix_can, bool verbose = false);
  void ImputeExteriorNode(int imp_idx, Data data_cur,
                          IntegerMatrix &gen_matrix_can, bool verbose = false);
  void ImputeInteriorNode(int imp_idx, Data data_cur, 
                          IntegerMatrix &gen_matrix_can, bool verbose = false);
  void ImputeColtimeDistances(Data data_cur, bool verbose = false);
  
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
  void UpdateExposureTimeUnif(int &nacc, bool verbose = false);
  
  void WriteOutputToFile(std::ofstream &myfile);
  void WriteConfigurationToFile(std::ofstream &myfile);
  
  // Tree search functions
  int CalculateDistanceBetweenNodes(int node_i, int node_j);
  int ReturnLCA(int node_i, int node_j);
  std::vector<int> ReturnPathToRoot(int node);
  std::vector<int> ReturnSurroundingNodes(int node);
  
};



// Constructor for the data class
Data::Data(NumericVector t_e_,
           NumericVector t_i_,
           NumericVector t_r_,
           IntegerVector source_,
           NumericVector x_,
           NumericVector y_,
           NumericVector parameters_)
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

  

  
  // parameters
  for(auto i : parameters_)
  {
    parameters.push_back(i);
  }
  
  
  CalculateLoglik();
  
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



void Data::CalculateLoglik(bool verbose)
{
  verbose = false;
  CalculateTransmissionLikelihood(verbose);
  CalculateRemovalLikelihood(verbose);
  CalculateExposureLikelihood(verbose);
  loglik = transmission_likelihood + removal_likelihood + exposure_likelihood;

  if(verbose)
  {
    std::cout << "Calculated log likelihood"
              << ", transmission loglik = " << transmission_likelihood
              << ", removal loglik = " << removal_likelihood
              << ", exposure loglik = " << exposure_likelihood
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
    
    // Reject early if outbreak is not possible
    if(t_e_can < t_i[target_source]) return;
    
    data_can.t_e[target] = t_e_can;
    
    // Determine a new source
    std::vector<int> possible_infectors_cur = ReturnPossibleInfectors(target);
    std::vector<double> weights_cur;
    double beta = parameters[1];
    double kappa = parameters[3];
    double weight_sum_cur = 0.0;
    for(auto i : possible_infectors_cur)
    {
      double cur_dist = CalculateDistance(i,target);
      double h_ij = exp(-1*kappa*cur_dist);
      weight_sum_cur += beta*h_ij;
    }
    
    for(auto i : possible_infectors_cur)
    {
      double cur_dist = CalculateDistance(i,target);
      double h_ij = exp(-1*kappa*cur_dist);
      weights_cur.push_back(beta*h_ij/weight_sum_cur);
    }

    
    int source_cur = source[target];

    //int cur_loc = WhichVec(source_cur, possible_infectors_cur)[0];
    std::vector<int> cur_locs = WhichVec(source_cur, possible_infectors_cur);

    int cur_loc = cur_locs[0];

    

    std::vector<int> possible_infectors_can = 
    data_can.ReturnPossibleInfectors(target);

    
    // Exit early if there are no infectors at the proposed infection time
    if(possible_infectors_can.size()==0) return;
    
    std::vector<double> weights_can;
    double weight_sum_can = 0.0;
    for(auto i : possible_infectors_can)
    {
      double cur_dist = CalculateDistance(i,target);
      double h_ij = exp(-1*kappa*cur_dist);
      weight_sum_can += beta*h_ij;
    }
    
    for(auto i : possible_infectors_can)
    {
      double cur_dist = CalculateDistance(i,target);
      double h_ij = exp(-1*kappa*cur_dist);
      weights_can.push_back(beta*h_ij/weight_sum_can);
    }

    
    
    int source_can = NonUnifSampling(possible_infectors_can, weights_can);
    int can_loc = WhichVec(source_can, possible_infectors_can)[0];

    
    log_prop_ratio += log(weights_cur[cur_loc]) - log(weights_can[can_loc]);
    for(auto i : weights_cur)
    {
      log_prop_ratio -= log(i);
    }
    for(auto i : weights_can)
    {
      log_prop_ratio += log(i);
    }
    
    
    
    
    data_can.source[target] = source_can;
    
    double latent_period_cur = t_i[target]-t_e[target];
    double upper = R::dgamma(latent_period_cur,zeta,1/alpha,1);
    double lower = R::dgamma(latent_period_can,zeta,1/alpha,1);
    log_prop_ratio += upper - lower;
    

    
    
    
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

      nacc++;
      
      
      
      

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

void Data::UpdateExposureTimeUnif(int &nacc, bool verbose)
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
    
    // Reject early if outbreak is not possible
    //if(t_e_can < t_i[target_source]) return;
    
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
      
      nacc++;
      
      
      
      
      
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
void MCMC_SEIR(List MCMC_options,
               NumericVector t_e,
               NumericVector t_i, 
               NumericVector t_r,
               IntegerVector source,
               NumericVector x,
               NumericVector y)
{
  
  std::cout << "got this far" << std::endl;
  // Load data from the options
  NumericVector parameters = MCMC_options["initial_chain_state"];
  std::cout << "got this far" << std::endl;
  int num_augmented_updates = MCMC_options["num_aug_updates"];
  int max_iterations = MCMC_options["iterations"];
  std::string output_file = MCMC_options["output_file"];
  List debug_flags = MCMC_options["debug_flags"];
  List proposal_variance = MCMC_options["proposal_variance"];
  List prior_parameters = MCMC_options["prior_parameters"];
  std::cout << "got this far" << std::endl;
  
  // Instantiate the data
  Data data(t_e,
            t_i,
            t_r,
            source,
            x,
            y,
            parameters);
  
  
  
  
  // Acceptance counters
  int nacc_kappa = 0;
  int nacc_kappa_prop = 0;
  int nacc_delta = 0;
  int nacc_delta_prop = 0;
  int nacc_zeta = 0;
  int nacc_zeta_prop = 0;

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
    

    
    
    // Update the augmented data
    int aug_data_flag = debug_flags["aug_data"];
    if(aug_data_flag==0)
    {
      for(int i = 0; i < num_augmented_updates; i++)
      {
        //int move = floor(R::runif(0,2));
        int move = 0;
        if(move==0)
        {
          nacc_exp_prop++;
          data.UpdateExposureTimeUnif(nacc_exp, false);
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
  double exp_prob = (double)nacc_exp/(double)nacc_exp_prop;
  double inf_prob = (double)nacc_inf/(double)nacc_inf_prop;
  double dist_prob = (double)nacc_dist/(double)nacc_dist_prop;
  
  std::cout << "Acceptance probabilities: kappa = " << kappa_prob
            << ", delta = " << delta_prob
            << ", zeta = " << zeta_prob
            << ", exp_time = " << exp_prob
            << ", inf_time = " << inf_prob
            << ", dist_prop = " << dist_prob
            << std::endl;
  
  myfile.close();
  myfile2.close();
}



