#include <iostream>
#include <fstream>
#include <Rcpp.h>

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

class Data {
public:
  // Epi data
  int N; // total number of individuals
  int n_I; // number of infected individuals
  
  // Vectors containing infection, removal and symptom times
  std::vector<double> t_i;
  std::vector<double> t_r;
  std::vector<double> t_s;
  
  // Source of infection
  std::vector<int> source;
  
  // Spatial coordinates
  std::vector<double> x;
  std::vector<double> y;
  
  // Index of infected individuals
  std::vector<int> ever_infected;
  
  
  // Likelihood
  double loglik;
  double removal_likelihood;
  double transmission_likelihood;
  
  std::vector<double> parameters;
  
  // Constructor
  Data(NumericVector t_i_,
       NumericVector t_r_,
       NumericVector t_s_,
       IntegerVector source_,
       NumericVector x_,
       NumericVector y_,
       NumericVector parameters_);
  
  
  // Functions
  std::vector<int> ReturnPossibleInfectors(int target);
  
  // Likelihood functions
  double CalculateDistance(int i, int j);
  double CalculateDistanceSum(bool verbose = false);
  double CalculateDoubleSum(bool verbose = false);
  void CalculateTransmissionLikelihood(bool verbose = false);
  void CalculateRemovalLikelihood(bool verbose = false);
  void CalculateLoglik(bool verbose = false);
  
  void UpdateParameterMetropolisHastingsUniform(int idx, 
                                                double proposal_variance,
                                                double prior_lower,
                                                double prior_upper,
                                                int &nacc, 
                                                bool verbose = false);
  
  // Augmented data functions
  void UpdateInfectionTime(int &nacc, bool verbose = false);
  void UpdateInfectionTimeGamma(int &nacc, bool verbose = false);
  
  void WriteOutputToFile(std::ofstream &myfile);
};

// Constructor for the data class
Data::Data(NumericVector t_i_,
           NumericVector t_r_,
           NumericVector t_s_,
           IntegerVector source_,
           NumericVector x_,
           NumericVector y_,
           NumericVector parameters_)
{
  // Initialise epi data
  N = t_i_.length();
  int total_infected = 0;
  
  for(int i = 0; i < N; i++)
  {
    t_i.push_back(t_i_[i]);
    t_r.push_back(t_r_[i]);
    t_s.push_back(t_s_[i]);
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
    
    if(t_i_[i] != -1)
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
    double infection_time = t_i[i];
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
      if(t_i[j] == -1)
      {
        // j is never infected
        contribution = h_ij*(t_r[i]-t_i[i]);
      }
      else
      {
        contribution = h_ij*(std::min(t_r[i],t_i[j]) - std::min(t_i[i],t_i[j]));
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


void Data::CalculateTransmissionLikelihood(bool verbose)
{
  verbose = false;
  double transmission_likelihood_ = 0.0;
  double beta = parameters[0];
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
  
  bool observational_contribution = true;
  if(observational_contribution)
  {
    for(auto person : ever_infected)
    {
      double inf_time = t_i[person];
      double rem_time = t_r[person];
      double symp_time = t_s[person];
      if(symp_time > inf_time && symp_time < rem_time)
      {
        transmission_likelihood_ -= log(t_r[person]-t_i[person]);
      }
      else
      {
        transmission_likelihood_ -= 1e16;
      }
      
    }
  }
  
  
  
  
  transmission_likelihood = transmission_likelihood_;
}

void Data::CalculateRemovalLikelihood(bool verbose)
{
  double delta = parameters[1];
  double gamma = parameters[2];
  double removal_likelihood_ = 0.0;
  for(auto i : ever_infected)
  {
    removal_likelihood_ += R::dgamma(t_r[i]-t_i[i],delta,1/gamma,1);
    //removal_likelihood_ -= R::pgamma(t_s[i],5,1/0.5,)
  }
  removal_likelihood = removal_likelihood_;
}

void Data::CalculateLoglik(bool verbose)
{
  CalculateTransmissionLikelihood();
  CalculateRemovalLikelihood();
  loglik = transmission_likelihood + removal_likelihood;
}

void Data::WriteOutputToFile(std::ofstream &myfile)
{
  for(size_t i = 0; i < parameters.size(); i++)
  {
    myfile << parameters[i] << " ";
  }
  
  myfile << loglik << " ";
  
  for(size_t i = 0; i < t_i.size(); i++)
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

std::vector<int> WhichVec(int x, const std::vector<int> vec)
{
  std::vector<int> out;
  for(size_t i = 0; i < vec.size(); i++)
  {
    if(vec[i] == x) out.push_back(i);
  }
  return out;
}

// Return the infectors associated with the exposure time
// of the target
std::vector<int> Data::ReturnPossibleInfectors(int target)
{
  std::vector<int> possible_infectors;
  double infection_time = t_i[target];
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

/*
void Data::UpdateInfectionTimeNoSource(int &nacc, bool verbose)
{
  
}
*/

void Data::UpdateInfectionTime(int &nacc, bool verbose)
{
  double log_prop_ratio = 0.0; // Metropolis hastings log proposal ratio
  
  // Create a candidate data set which is a copy of the data
  Data data_can((*this));
  
  // Uniformly at random choose an individual to update
  int target = SampleVector(ever_infected);
  int target_source = source[target];
  if(target_source != -1)
  {
    double last_time = t_r[target];
    
    // Check when the first offspring occurs as we cannot move the time
    // to be later than this
    std::vector<int> offspring = WhichVec(target,source);
    for(auto i : offspring)
    {
      last_time = std::min(last_time,t_i[i]);
    }
    
    // Propose an infection time
    double U = R::runif(0,1);
    data_can.t_i[target] = last_time*U;
    
    // Determine a new source
    std::vector<int> possible_infectors_cur = ReturnPossibleInfectors(target);
    std::vector<int> possible_infectors_can = 
      data_can.ReturnPossibleInfectors(target);
    
    // Exit early if there are no infectors at the proposed infection time
    if(possible_infectors_can.size()==0) return;
    
    // Uniformly sample from the available sources
    int source_can = SampleVector(possible_infectors_can);
    data_can.source[target] = source_can;
    
    log_prop_ratio = log(possible_infectors_can.size()) -
      log(possible_infectors_cur.size());
    
   
   
    data_can.CalculateLoglik();
    
    U = R::runif(0,1);
    if(log(U) < data_can.loglik - loglik + log_prop_ratio)
    {
      if(verbose)
      {
        std::cout << " Loglik cur = " << loglik << ", loglik can = "
                  << data_can.loglik << " - accepted." << std::endl;
      }
      t_i = data_can.t_i;
      source = data_can.source;
      loglik = data_can.loglik;
      nacc++;
    }
    else
    {
      if(verbose)
      {
        std::cout << " Loglik cur = " << loglik << ", loglik can = "
                  << data_can.loglik << " - rejected." << std::endl;
      }
    }
  }
}

void Data::UpdateInfectionTimeGamma(int &nacc, bool verbose)
{
  verbose = false;
  double log_prop_ratio = 0.0; // Metropolis hastings log proposal ratio
  
  // Create a candidate data set which is a copy of the data
  Data data_can((*this));
  
  // Uniformly at random choose an individual to update
  int target = SampleVector(ever_infected);
  int target_source = source[target];
  if(target_source != -1)
  {
    double last_time = t_r[target];
    
    // Check when the first offspring occurs as we cannot move the time
    // to be later than this
    std::vector<int> offspring = WhichVec(target,source);
    for(auto i : offspring)
    {
      last_time = std::min(last_time,t_i[i]);
    }
    
    last_time = std::min(last_time,t_s[target]);
    
    // Propose an infection time
    double delta = parameters[1];
    double gamma = parameters[2];
    double inf_period_can = R::rgamma(delta,1/gamma);
    
    data_can.t_i[target] = t_r[target] - inf_period_can;
    
    // Exit early if the proposed time is greater than the last time,
    // i.e. the outbreak would not be valid
    if(data_can.t_i[target] > last_time) return;
    
    
    // Determine a new source
    std::vector<int> possible_infectors_cur = ReturnPossibleInfectors(target);
    std::vector<int> possible_infectors_can = 
      data_can.ReturnPossibleInfectors(target);
    
    // Exit early if there are no infectors at the proposed infection time
    if(possible_infectors_can.size()==0) return;
    
    // Uniformly sample from the available sources
    int source_can = SampleVector(possible_infectors_can);
    data_can.source[target] = source_can;
    double inf_period_cur = t_r[target]-t_i[target];
    double upper = R::dgamma(inf_period_cur,delta,1/gamma,1) + 
      log(possible_infectors_can.size());
    double lower = R::dgamma(inf_period_can,delta,1/gamma,1) + 
      log(possible_infectors_cur.size());
    log_prop_ratio = upper - lower;
    
    if(false)
    {
      std::cout << "log prop ratio components:"
                << " R::dgamma(inf_period_cur,delta,1/gamma,1) = "
                << R::dgamma(inf_period_cur,delta,1/gamma,1)
                << ", log(possible_infectors_cur.size()) = "
                << log(possible_infectors_cur.size())
                << ", R::dgamma(inf_period_can,delta,1/gamma,1) = "
                << R::dgamma(inf_period_can,delta,1/gamma,1)
                << ", log(possible_infectors_can.size()) = "
                << log(possible_infectors_can.size())
                << std::endl;
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
void MCMC_SIR_E(List MCMC_options, 
                NumericVector t_i, 
                NumericVector t_r,
                NumericVector t_s,
                IntegerVector source,
                NumericVector x,
                NumericVector y)
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
  Data data(t_i,
            t_r,
            t_s,
            source,
            x,
            y,
            parameters);
  
  // Acceptance counters
  int nacc_beta = 0;
  int nacc_beta_prop = 0;
  int nacc_delta = 0;
  int nacc_delta_prop = 0;
  int nacc_gamma = 0;
  int nacc_gamma_prop = 0;
  int nacc_kappa = 0;
  int nacc_kappa_prop = 0;
  int nacc_inf = 0;
  int nacc_inf_prop = 0;
  
  
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
    // Update beta by gibbs
    int beta_flag = debug_flags["beta"];
    if(beta_flag==0)
    {
      
      double double_sum = data.CalculateDoubleSum();
      double prior_shape = prior_parameters["beta_shape"];
      double prior_rate = prior_parameters["beta_rate"];
      data.parameters[0] = R::rgamma(data.n_I + prior_shape,
                                     1/(prior_rate + double_sum));
      /*
      double prop_var = proposal_variance["beta"];
      nacc_beta_prop++;
      data.UpdateParameterMetropolisHastingsUniform(0, 
                                                    prop_var,
                                                    0,
                                                    0.5,
                                                    nacc_beta);
       */
    }
    
    // Update gamma by gibbs
    int gamma_flag = debug_flags["gamma"];
    if(gamma_flag==0)
    {
      
      double infectious_period_sum = 0.0;
      double delta = parameters[1];
      double prior_shape = prior_parameters["gamma_shape"];
      double prior_rate = prior_parameters["gamma_rate"];
      for(auto i : data.ever_infected)
      {
        infectious_period_sum += (data.t_r[i]-data.t_i[i]);
      }
      data.parameters[2] = R::rgamma(delta*data.n_I + prior_shape,
                                     1/(prior_rate + infectious_period_sum));
      /*
      double prop_var = proposal_variance["gamma"];
      nacc_gamma_prop++;
      data.UpdateParameterMetropolisHastingsUniform(2, 
                                                    prop_var,
                                                    0,
                                                    10,
                                                    nacc_gamma);
      */
    }
    
    // Update kappa by MH
    int kappa_flag = debug_flags["kappa"];
    if(kappa_flag==0)
    {
      double prop_var = proposal_variance["kappa"];
      nacc_kappa_prop++;
      data.UpdateParameterMetropolisHastingsUniform(3, 
                                                    prop_var,
                                                    0,
                                                    1,
                                                    nacc_kappa);
    }
    
    // Update delta by MH
    int delta_flag = debug_flags["delta"];
    if(delta_flag==0)
    {
      double prop_var = proposal_variance["delta"];
      nacc_delta_prop++;
      data.UpdateParameterMetropolisHastingsUniform(1, 
                                                    prop_var,
                                                    0,
                                                    10,
                                                    nacc_delta);
    }
    
    // Update the augmented data
    int aug_data_flag = debug_flags["aug_data"];
    if(aug_data_flag==0)
    {
      for(int i = 0; i < num_augmented_updates; i++)
      {
        int move = floor(R::runif(1,2));
        if(move==0)
        {
          nacc_inf_prop++;
          data.UpdateInfectionTime(nacc_inf);
        }
        else if(move == 1)
        {
          nacc_inf_prop++;
          data.UpdateInfectionTimeGamma(nacc_inf);
        }
      }
    }
    
    
    
    data.WriteOutputToFile(myfile);
    
    std::cout << "tranmission_ll = " << data.transmission_likelihood
              << ", removal_ll = " << data.removal_likelihood
              << std::endl;
    
    if(iter%100==0)
    {
      std::cout << "Iteration " << iter << " completed." << std::endl;
      checkUserInterrupt();
    }
  }
  
  double beta_prob = (double)nacc_beta/(double)nacc_beta_prop;
  double delta_prob = (double)nacc_delta/(double)nacc_delta_prop;
  double gamma_prob = (double)nacc_gamma/(double)nacc_gamma_prop;
  double kappa_prob = (double)nacc_kappa/(double)nacc_kappa_prop;
  double inf_prob = (double)nacc_inf/(double)nacc_inf_prop;
  
  std::cout << "Acceptance probabilities: beta = " << beta_prob
            << ", delta = " << delta_prob
            << ", gamma = " << gamma_prob
            << ", kappa = " << kappa_prob
            << ", inf_time = " << inf_prob
            << std::endl;
}



