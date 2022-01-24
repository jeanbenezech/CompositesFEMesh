"""

"""
import numpy as np
from math import exp, sqrt
from scipy import stats
from node_element_extraction import cov_mat_eval
import time


class mcmc_sim():
    
    def __init__(self, model, prior_dist_hyper_param, data):
        
        self.model = model
        self.prior_dist_hyper_param = prior_dist_hyper_param
        self.data = data
        self.dim_in_white_noise = self.model.white_noise_size
        self.dim_in_hyper_param = len(prior_dist_hyper_param)
        self.dim_out = len(model.output_node_indices)
        self.rng = np.random.default_rng()
        
    def simulate(self, draws=1000, burnin=1000, initial_white_noise=None, initial_hyper_param=None, beta_white_noise=0.3, beta_hyper_param=0.3, accept_ratio_target=0.33, tune_interval=200):
        '''
        Get the posterior samples by running a Markov chain
        '''
        num_draws = draws + burnin
        self.beta_white_noise = beta_white_noise
        self.beta_hyper_param = beta_hyper_param
        self.accept_ratio_target = accept_ratio_target
        self.tune_interval = tune_interval
        self.num_accept_white_noise = 0
        self.num_accept_white_noise_recently = 0
        self.num_accept_hyper_param = 0
        self.num_accept_hyper_param_recently = 0
        self.num_accept_white_noise_after_tune = 0
        self.num_accept_hyper_param_after_tune = 0
        self.sample_mcmc_white_noise = np.ones((num_draws, self.dim_in_white_noise))
        self.sample_mcmc_hyper_param = np.ones((num_draws, self.dim_in_hyper_param))
        self.log_likelihood_value = np.ones(num_draws)
        self.sample_mcmc_output = np.ones((num_draws, self.dim_out))
        
        ## initialization
        current_white_noise, current_hyper_param_normal, current_hyper_param_target, current_log_likelihood, current_model_output = self.initalization(initial_white_noise, initial_hyper_param)
        
        ## iterations
        start = time.time()
        for iteration in range(num_draws):
            
            # update white noise
            current_white_noise, current_log_likelihood, current_model_output, accepted_white_noise = self.step_white_noise(current_white_noise, current_hyper_param_target, current_log_likelihood, current_model_output)
            
            self.num_accept_white_noise += accepted_white_noise
            self.num_accept_white_noise_recently += accepted_white_noise
            if iteration >= burnin:
                self.num_accept_white_noise_after_tune += accepted_white_noise
            print("accept(white_noise):", self.num_accept_white_noise)
            
            # update hyper parameters
            current_hyper_param_normal, current_hyper_param_target, current_log_likelihood, current_model_output, accepted_hyper_param = self.step_hyper_param(current_hyper_param_normal, current_hyper_param_target, current_white_noise, current_log_likelihood, current_model_output)
            
            self.num_accept_hyper_param += accepted_hyper_param
            self.num_accept_hyper_param_recently += accepted_hyper_param
            if iteration >= burnin: 
                self.num_accept_hyper_param_after_tune += accepted_hyper_param
            print("accept(hyper_param):", self.num_accept_hyper_param)
            
            # update beta value based on recent acceptance ratio
            if ((iteration+1) % self.tune_interval ==0) & (iteration<burnin) :
                # recent acceptance ratio
                alpha_white_noise_recently = self.num_accept_white_noise_recently/self.tune_interval
                alpha_hyper_param_recently = self.num_accept_hyper_param_recently/self.tune_interval
                # update beta value
                self.beta_white_noise = self.tune(self.beta_white_noise, alpha_white_noise_recently)
                self.beta_hyper_param = self.tune(self.beta_hyper_param, alpha_hyper_param_recently)
                # reset recent number of accepted sample
                self.num_accept_white_noise_recently = 0
                self.num_accept_hyper_param_recently = 0
            
            # save the current sample, likelihood and posterior value
            self.sample_mcmc_white_noise[iteration,:] = current_white_noise
            self.sample_mcmc_hyper_param[iteration,:] = current_hyper_param_target
            self.log_likelihood_value[iteration] = current_log_likelihood
            self.sample_mcmc_output[iteration,:] = current_model_output
            print('Iteration: %d in %d, Elapsed time: %g' % (iteration, num_draws, time.time()-start))
        
        # return self.sample_mcmc, self.log_likelihood_value
        
    def initalization(self, initial_white_noise=None, initial_hyper_param=None):
        ''' 
        Initalization. Set the initial sample for the Markov chain
        If initial sample is None, use the random sample as initial sample
        '''
        
        # white noise
        if initial_white_noise is None:
            current_white_noise = self.rng.standard_normal(self.dim_in_white_noise) # random initial white noise
        else:
            if len(initial_white_noise) != self.dim_in_white_noise:
                raise ValueError("size of initial sample for white noise is not correct")
            current_white_noise = initial_white_noise
        
        # hyper parameters
        if initial_hyper_param is None:
            current_log_likelihood = None
            while current_log_likelihood is None:
                while True:
                    current_hyper_param_normal = self.rng.standard_normal(self.dim_in_hyper_param)  # random initial sample in standard normal space
                    
                    current_hyper_param_target = np.ones(self.dim_in_hyper_param)
                    for i in range(self.dim_in_hyper_param):
                        current_hyper_param_target[i] = self.prior_dist_hyper_param[i].ppf(stats.norm.cdf(current_hyper_param_normal[i])) # in target space
                    
                    if current_hyper_param_target[0] > current_hyper_param_target[1]: # constraints between the first two parameters
                        break
                current_log_likelihood, current_model_output = self.log_likelihood(current_white_noise, current_hyper_param_target)
        else:
            if len(initial_hyper_param) != self.dim_in_hyper_param:
                raise ValueError("size of initial sample for hyper-parameters is not correct!")
            elif initial_hyper_param[1] > initial_hyper_param[0]:
                raise ValueError("first hyper-parameter should be larger than the second one!")
            else:
                current_hyper_param_target = initial_hyper_param  # in target space
                
                current_hyper_param_normal = np.ones(self.dim_in_hyper_param) # in standard normal space
                for i in range(self.dim_in_hyper_param):
                    current_hyper_param_normal[i] = stats.norm.ppf(self.prior_dist_hyper_param[i].cdf(current_hyper_param_target[i]))
        
            # calculate the likelihood value of current sample
            current_log_likelihood, current_model_output = self.log_likelihood(current_white_noise, current_hyper_param_target)
            if current_log_likelihood is None:
                raise ValueError("initial sample of hyper-parameter cannot produce a reasonable model output!")
            
        return current_white_noise, current_hyper_param_normal, current_hyper_param_target, current_log_likelihood, current_model_output
        
    def model_run(self, white_noise, hyper_param):
        
        return self.model.model_run(Ylength=hyper_param[0], ThicknessVar=hyper_param[1], Sigma=hyper_param[2], Length=hyper_param[3], white_noise=white_noise)
        
    def log_likelihood(self, white_noise, hyper_param):
        '''
        Calculate the log likelihood value given a parameter value
        '''
        if self.model.gaussian_thickness == 1:
            # get the model output
            model_output = self.model_run(white_noise, hyper_param[:-1])  # hyper parameters except the last one
            # calculate the likelihood
            if (model_output == 1e5).all():
                return None, None
            else:
                # log_like = stats.t.logpdf(model_output - self.data, df=3, scale=hyper_param[-1]).sum()
                log_like = stats.norm.logpdf(model_output - self.data, scale=hyper_param[-1]).sum()
                return log_like, model_output
        elif self.model.gaussian_thickness == 0:
            model_output = self.model_run(white_noise, hyper_param[0:2])
            
            if (model_output == 1e5).all():
                return None, None
            else:
                cov_mat = cov_mat_eval(self.model.output_node_coords_2d, hyper_param[-2], hyper_param[-1])
                log_like = stats.multivariate_normal.logpdf(model_output-self.data, cov=cov_mat).sum()
                return log_like, model_output
        
    def get_sample_white_noise(self, current_white_noise):
        '''
        Propose a candidate sample of white noise with the pCN proposal
        '''
        return  sqrt(1-self.beta_white_noise**2)*current_white_noise + self.beta_white_noise*self.rng.standard_normal(self.dim_in_white_noise)
        
    def get_sample_hyper_param(self, current_hyper_param_normal):
        '''
        Propose a candidate sample of hyper parameters with the pCN proposal
        '''
        return sqrt(1-self.beta_hyper_param**2)*current_hyper_param_normal + self.beta_hyper_param*self.rng.standard_normal(self.dim_in_hyper_param)
        
    def step_white_noise(self, current_white_noise, current_hyper_param_target, current_log_likelihood, current_model_output):
        '''
        Get a candidate sample and determine whether to accept it or not
        '''
        
        # get a candidate sample which can produce a reasonable model output (not None)
        candidate_log_likelihood = None
        while candidate_log_likelihood is None:
            # propose a candidate sample
            candidate_white_noise = self.get_sample_white_noise(current_white_noise)
            
            # calculate the loglikelihood of candidate sample
            candidate_log_likelihood, candidate_model_output = self.log_likelihood(candidate_white_noise, current_hyper_param_target)
        
        # calculate the ratio
        if candidate_log_likelihood > current_log_likelihood: # definitely accept
            ratio = 1
        else:
            ratio = exp(candidate_log_likelihood - current_log_likelihood)
            
        print("ratio(white_noise):", ratio)
        accepted_white_noise = 0
        if ratio > self.rng.random(): # accept the candidate sample
            current_white_noise = candidate_white_noise
            current_log_likelihood = candidate_log_likelihood
            current_model_output = candidate_model_output
            accepted_white_noise = 1
            
        return current_white_noise, current_log_likelihood, current_model_output, accepted_white_noise
        
    def step_hyper_param(self, current_hyper_param_normal, current_hyper_param_target, current_white_noise, current_log_likelihood, current_model_output):
        '''
        Get a candidate sample and determine whether to accept it or not
        '''
        
        # get a candidate sample which can produce a reasonable model output (not None)
        candidate_log_likelihood = None
        while candidate_log_likelihood is None:
            # propose a candidate sample
            while True:
                candidate_hyper_param_normal = self.get_sample_hyper_param(current_hyper_param_normal)
                
                candidate_hyper_param_target = np.ones(self.dim_in_hyper_param)  # in target dpace
                for i in range(self.dim_in_hyper_param):
                    candidate_hyper_param_target[i] = self.prior_dist_hyper_param[i].ppf(stats.norm.cdf(candidate_hyper_param_normal[i]))
                    
                if candidate_hyper_param_target[0] > candidate_hyper_param_target[1]: # restriction between first two parameters
                    break
            
            # calculate the loglikelihood of candidate sample
            candidate_log_likelihood, candidate_model_output = self.log_likelihood(current_white_noise, candidate_hyper_param_target)
        
        # calculate the ratio
        if candidate_log_likelihood > current_log_likelihood:
            ratio = 1
        else:
            ratio = np.exp(candidate_log_likelihood - current_log_likelihood)
            
        print("ratio(hyper_param):", ratio)
        accepted_hyper_param = 0
        if ratio > self.rng.random(): # accept the candidate sample
            current_hyper_param_normal = candidate_hyper_param_normal
            current_hyper_param_target = candidate_hyper_param_target
            current_log_likelihood = candidate_log_likelihood
            current_model_output = candidate_model_output
            accepted_hyper_param = 1
            
        return current_hyper_param_normal, current_hyper_param_target, current_log_likelihood, current_model_output, accepted_hyper_param
    
    def tune(self, beta, alpha_recently):
        '''
        Tune the jumping factor beta in the proposal
        '''
        # parameter tuning 
        alpha_dicrepancy = alpha_recently/self.accept_ratio_target
        if alpha_dicrepancy>1.5:
            beta = beta*1.25
        elif alpha_dicrepancy>1.2:
            beta = beta*1.1
        elif alpha_dicrepancy<0.5:
            beta = beta*0.75
        elif alpha_dicrepancy<0.8:
            beta = beta*0.9
        
        if beta>1:
            beta = 0.9999
        
        return beta
        
if __name__ == '__main__':
    
    from c_spar_shell import c_spar_shell_model
    import node_element_extraction
    import os
    
    c_spar_shell_model()  # get initial '.msh' file for node extraction
    node_indices, node_coords = node_element_extraction.get_output_nodes()
    model_c_spar = c_spar_shell_model(gaussian_thickness=1, corner_thickness=1, output_node_indices=node_indices, output_node_coords_2d=node_coords)
    
    
    prior_dist_hyper_param = []
    prior_dist_hyper_param.append(stats.lognorm(s=0.101, scale=np.exp(1.844))) # Thickness
    prior_dist_hyper_param.append(stats.lognorm(s=0.236, scale=np.exp(0.549))) # ThicknessVar
    #prior_dist_hyper_param.append(stats.lognorm(s=0.644, scale=np.exp(-2.414))) # Sigma_Cov_Thickness
    prior_dist_hyper_param.append(stats.expon(scale=0.087)) # Sigma_Cov_Thickness
    prior_dist_hyper_param.append(stats.lognorm(s=0.236, scale=np.exp(2.159))) # Length_Cov_Thickness
    prior_dist_hyper_param.append(stats.lognorm(s=0.495, scale=np.exp(-4.658))) # Sigma_Error
    
    
    data_disp_select_z = np.loadtxt(os.path.join('.','data','data_disp_select_z_20_all.csv'))
    data_disp_select_z_s = data_disp_select_z[0:3,:]
    
    mcmc_c_spar = mcmc_sim(model_c_spar, prior_dist_hyper_param, data_disp_select_z_s)
    
    
    mcmc_c_spar.simulate(draws=20, burnin=0, initial_white_noise=None, initial_hyper_param=None, beta_white_noise=0.3, beta_hyper_param=0.3, accept_ratio_target=0.30, tune_interval=300)
    
    np.savetxt(os.path.join('.','data','post_samples_white_noise_20_exp_all.csv'), mcmc_c_spar.sample_mcmc_white_noise)
    np.savetxt(os.path.join('.','data','post_samples_hyper_param_20_exp_all.csv'), mcmc_c_spar.sample_mcmc_hyper_param)
    np.savetxt(os.path.join('.','data','post_model_output_20_exp_all.csv'), mcmc_c_spar.sample_mcmc_output)
    np.savetxt(os.path.join('.','data','log_likelihood_20_exp_all.csv'), mcmc_c_spar.log_likelihood_value)
    
    f = open(os.path.join('.', 'data', 'mcmc_stats_20_exp_all.csv'), 'w')
    f.write("beta_hyper_param: {}\n".format(mcmc_c_spar.beta_hyper_param))
    f.write("beta_white_noise: {}\n".format(mcmc_c_spar.beta_white_noise))
    f.write("num_accept_hyper_param: {}\n".format(mcmc_c_spar.num_accept_hyper_param))
    f.write("num_accept_white_noise: {}\n".format(mcmc_c_spar.num_accept_white_noise))
    f.write("num_accept_hyper_param_after_tune: {} \n".format(mcmc_c_spar.num_accept_hyper_param_after_tune))
    f.write("num_accept_white_noise_after_tune: {}\n".format(mcmc_c_spar.num_accept_white_noise_after_tune))
    f.close()
