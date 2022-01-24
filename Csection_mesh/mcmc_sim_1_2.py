"""
Created on Wed Dec 15 19:25:31 2021

@author: sx450
"""

import numpy as np
from scipy import stats
from node_element_extraction import cov_mat_eval
import time

class mcmc_sim():
    
    def __init__(self, model, prior_dist, data):
        
        self.model = model
        self.prior_dist = prior_dist
        self.data = data
        self.dim_in = len(prior_dist)  # number of parameters
        self.dim_out = len(model.output_node_indices)
        self.rng = np.random.default_rng()
        
    def simulate(self, draws=1000, burnin=1000, initial_sample=None, beta=0.3, accept_ratio_target=0.33, tune_interval=200):
        '''
        Get the posterior samples by running a Markov chain
        '''
        num_draws = draws + burnin
        self.beta = beta
        self.accept_ratio_target = accept_ratio_target
        self.tune_interval = tune_interval
        self.num_accept = 0
        self.num_accept_recently = 0
        self.num_accept_after_tune = 0
        self.sample_mcmc = np.ones((num_draws, self.dim_in))
        self.log_likelihood_value = np.ones(num_draws)
        self.sample_mcmc_output = np.ones((num_draws, self.dim_out))
        
        ## initialization
        current_sample_normal, current_sample_target, current_log_likelihood, current_logpdf, current_model_output = self.initalization(initial_sample)
        
        ## iterations
        start = time.time()
        for iteration in range(num_draws):
            
            current_sample_normal, current_sample_target, current_log_likelihood, current_logpdf, current_model_output, accepted = self.step(current_sample_normal, current_sample_target, current_log_likelihood, current_logpdf, current_model_output)
            
            # print(current_model_output)
            
            self.num_accept = self.num_accept + accepted
            self.num_accept_recently = self.num_accept_recently + accepted
            if iteration >= burnin: 
                self.num_accept_after_tune = self.num_accept_after_tune + accepted
            print("accept:", self.num_accept)
            
            # update beta value based on recent acceptance ratio
            if ((iteration+1) % self.tune_interval ==0) & (iteration<burnin) :
                # recent acceptance ratio
                alpha_recently = self.num_accept_recently/self.tune_interval
                # update beta value
                self.beta = self.tune(self.beta, alpha_recently)
                # reset recent number of accepted sample
                self.num_accept_recently = 0
            
            # save the current sample, likelihood and posterior value
            self.sample_mcmc[iteration,:]        = current_sample_target
            self.log_likelihood_value[iteration] = current_log_likelihood
            self.sample_mcmc_output[iteration,:] = current_model_output
            print('Iteration: %d in %d, Elapsed time: %g' % (iteration, num_draws, time.time()-start))
        
        # return self.sample_mcmc, self.log_likelihood_value
        
    def initalization(self, initial_sample=None):
        ''' 
        Initalization. Set the initial sample for the Markov chain
        If 'initial_sample' is None, use the random sample as initial sample
        '''
        
        if initial_sample is None:
            current_log_likelihood = None
            while current_log_likelihood is None:
                while True:
                    current_sample_normal = self.rng.standard_normal(self.dim_in)  # random initial sample in standard normal space
                    
                    current_sample_target = np.ones(self.dim_in)
                    for i in range(self.dim_in):
                        current_sample_target[i] = self.prior_dist[i].ppf(stats.norm.cdf(current_sample_normal[i])) # in target space
                    
                    if current_sample_target[0] > current_sample_target[1]: # constraints between the first two parameters
                        break
                current_log_likelihood, current_model_output = self.log_likelihood(current_sample_target)
        else:
            if len(initial_sample) != self.dim_in:
                raise ValueError("size of initial sample is not correct!")
            elif initial_sample[1] > initial_sample[0]:
                raise ValueError("first parameter should be larger than the second one!")
            else:
                current_sample_target = initial_sample  # in target space
                
                current_sample_normal = np.ones(self.dim_in) # in standard normal space
                for i in range(self.dim_in):
                    current_sample_normal[i] = stats.norm.ppf(self.prior_dist[i].cdf(current_sample_target[i]))
        
            # calculate the likelihood value of current sample
            current_log_likelihood, current_model_output = self.log_likelihood(current_sample_target)
        
        # calculate the logpdf value of current sample
        current_logpdf = 0
        for i in range(self.dim_in):
            current_logpdf = current_logpdf + stats.norm.logpdf(current_sample_normal[i])
            
        return current_sample_normal, current_sample_target, current_log_likelihood, current_logpdf, current_model_output
        
    def model_run(self, parameters):
        '''
        run c-spar model
        '''
        if self.model.gaussian_thickness == 1:
            return self.model.model_run(Ylength=parameters[0], ThicknessVar=parameters[1], Sigma=parameters[2], Length=parameters[3])
        elif self.model.gaussian_thickness == 0:
            return self.model.model_run(Ylength=parameters[0], ThicknessVar=parameters[1])
        
    def log_likelihood(self, parameters):
        '''
        Calculate the log likelihood value given a parameter value
        '''
        if self.model.gaussian_thickness == 1:
            # get the model output
            model_output = self.model_run(parameters[:-1])  # model parameters except the last one
            
            # calculate the likelihood
            if (model_output == 1e5).all():
                return None, None
            else:
                # log_like = stats.t.logpdf(model_output - self.data, df=3, scale=parameters[-1]).sum()
                log_like = stats.norm.logpdf(model_output - self.data, scale=parameters[-1]).sum()
                return log_like, model_output
        elif self.model.gaussian_thickness == 0:
            model_output = self.model_run(parameters[0:2])
            
            if (model_output == 1e5).all():
                return None, None
            else:
                cov_mat = cov_mat_eval(self.model.output_node_coords_2d, parameters[-2], parameters[-1])
                log_like = stats.multivariate_normal.logpdf(model_output - self.data, cov=cov_mat).sum()
                return log_like, model_output
    
    def get_sample_mh(self, current_sample_normal): # random walk
        '''
        Propose a candidate sample with the MH(random walk) proposal
        '''
        return current_sample_normal + self.beta*self.rng.standard_normal(self.dim_in)
        
    def get_sample_pcn(self, current_sample_normal): # pCN proposal
        '''
        Propose a candidate sample with the pCN proposal
        '''
        return  np.sqrt(1-self.beta**2)*current_sample_normal + self.beta*self.rng.standard_normal(self.dim_in)
        
    def step(self, current_sample_normal, current_sample_target, current_log_likelihood, current_logpdf, current_model_output):
        '''
        Get a candidate sample and determine accept it or not
        '''
        
        # get a candidate sample which can produce a reasonable model output (not None)
        candidate_log_likelihood = None
        while candidate_log_likelihood is None:
            # propose a candidate sample
            while True:
                # candidate_sample_normal = self.get_sample_mh(current_sample_normal)  # in standard normal space
                candidate_sample_normal = self.get_sample_pcn(current_sample_normal)
                
                candidate_sample_target = np.ones(self.dim_in)  # in target dpace
                for i in range(self.dim_in):
                    candidate_sample_target[i] = self.prior_dist[i].ppf(stats.norm.cdf(candidate_sample_normal[i]))
                    
                if candidate_sample_target[0] > candidate_sample_target[1]: # restriction between first two parameters
                    break
            
            # calculate the loglikelihood of candidate sample
            candidate_log_likelihood, candidate_model_output = self.log_likelihood(candidate_sample_target)
        
        # calculate the logpdf of candidate sample
        candidate_logpdf = 0
        for i in range(self.dim_in):
            candidate_logpdf = candidate_logpdf + stats.norm.logpdf(candidate_sample_normal[i])
            
        # calculate the ratio
        # ratio = np.exp(candidate_log_likelihood + candidate_logpdf - current_log_likelihood - current_logpdf)
        if candidate_log_likelihood > current_log_likelihood:  # avoid overflow of 'np.exp()'
            ratio = 1
        else:
            ratio = np.exp(candidate_log_likelihood - current_log_likelihood)
            
        print("ratio:", ratio)
        accepted = 0
        if ratio > self.rng.random(): # accept the candidate sample
            current_sample_normal = candidate_sample_normal
            current_sample_target = candidate_sample_target
            current_log_likelihood = candidate_log_likelihood
            current_logpdf = candidate_logpdf
            current_model_output = candidate_model_output
            accepted = 1
            
        return current_sample_normal, current_sample_target, current_log_likelihood, current_logpdf, current_model_output, accepted
        
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
    
    # c_spar_shell_model()  # get initial '.msh' file for node extraction
    # node_indices, node_coords = node_element_extraction.get_output_nodes()
    model_c_spar = c_spar_shell_model(gaussian_thickness=1, corner_thickness=1, output_node_indices=None, output_node_coords_2d=None)
    
    
    prior_dist = []
    prior_dist.append(stats.lognorm(s=0.101, scale=np.exp(1.844))) # Thickness
    prior_dist.append(stats.lognorm(s=0.236, scale=np.exp(0.549))) # ThicknessVar
    
    # model 1
    # prior_dist.append(stats.lognorm(s=0.644, scale=np.exp(-2.414))) # Sigma_Cov_Thickness
    prior_dist.append(stats.expon(scale=0.087)) # Sigma_Cov_Thickness
    prior_dist.append(stats.lognorm(s=0.236, scale=np.exp(2.159))) # Length_Cov_Thickness
    prior_dist.append(stats.lognorm(s=0.495, scale=np.exp(-4.658))) # Sigma_Error
    
    # # model 2
    # # prior_dist.append(stats.lognorm(s=0.644, scale=np.exp(-4.716))) # Sigma_GP_Error
    # prior_dist.append(stats.expon(scale=0.009)) # Sigma_GP_Error
    # prior_dist.append(stats.lognorm(s=0.2, scale=np.exp(2.0))) # Length_GP_Error
    
    
    data_disp_select_z = np.loadtxt(os.path.join('.','data','data_disp_select_z_full.csv'))
    data_disp_select_z_s = data_disp_select_z[0,:]
    
    mcmc_c_spar = mcmc_sim(model_c_spar, prior_dist, data_disp_select_z_s)
    
    # initial_sample = np.array([5.54412,4.49342,0.187711,3.66536,0.00939566])

    mcmc_c_spar.simulate(draws=7000, burnin=3000, initial_sample=None, beta=0.3, accept_ratio_target=0.30, tune_interval=300)
    
    np.savetxt(os.path.join('.','data','post_samples_full_exp.csv'), mcmc_c_spar.sample_mcmc)
    np.savetxt(os.path.join('.','data','log_likelihood_full_exp.csv'), mcmc_c_spar.log_likelihood_value)
    np.savetxt(os.path.join('.','data','post_samples_output_full_exp.csv'), mcmc_c_spar.sample_mcmc_output)
    f = open(os.path.join('.','data','beta_full_exp.txt'), 'w')
    f.write("beta: {}\n".format(mcmc_c_spar.beta))
    f.write("num_accept: {}\n".format(mcmc_c_spar.num_accept))
    f.write("num_accept_after_tune: {}\n".format(mcmc_c_spar.num_accept_after_tune))
    f.close()
