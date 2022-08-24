##This explores simple models of selection, simulating the allele frqeuency trajectory forward in time.

from __future__ import print_function

import numpy as np
import pylab as pl
import copy
import argparse

def selection_deterministic(initial_frequency = 0.01, population_size = 100, selection_coefficient = 0.01, time_steps = 100):
    """
    This function returns the deterministic frequency trace of an allele under selection given an initial frequency, selection coefficient, and number of time_steps.
    """
    assert initial_frequency >= 0.0 and initial_frequency <= 1.0
    assert population_size > 0
    assert selection_coefficient >= 0
    assert time_steps > 0

    frequencies_per_generation = [initial_frequency]
    for generation in range(time_steps - 1):
        #Get the frequency in the generation t
        frequency_in_current_generation = frequencies_per_generation[-1]
        #Calculate the number of individuals of both the selected (fitness 1 + s) and normal (fitness 1) alleles. Note that we assume population size is constant.
        number_after_reproduction_selected_allele = (frequency_in_current_generation) * population_size * (1.0 + selection_coefficient)
        number_after_reproduction_normal_allele = (1.0 - frequency_in_current_generation) * population_size * (1.0)
        #Convert this to a new frequency.
        frequency_after_population_regulation = number_after_reproduction_selected_allele / (number_after_reproduction_selected_allele + number_after_reproduction_normal_allele)
        #Record the new frequency in generation t + 1
        frequencies_per_generation.append(frequency_after_population_regulation)
    #(Convert to a numpy array for efficient manipulation)
    frequencies_per_generation = np.array(frequencies_per_generation)
    return frequencies_per_generation

def selection_stochastic_haploid(initial_frequency = 0.01, population_size = 100, selection_coefficient = 0.01, time_steps = 100):
    """
    This function returns a stochastic frequency trace of an allele under selection given an initial frequency, selection coefficient, and number of time_steps.
    """
    assert initial_frequency >= 0.0 and initial_frequency <= 1.0
    assert population_size > 0
    assert selection_coefficient >= 0
    assert time_steps > 0
    
    frequencies_per_generation = [initial_frequency]
    for generation in range(time_steps):
        #Get the frequency in the generation t
        frequency_in_current_generation = frequencies_per_generation[-1]
        #Draw a new population based on the cumulative fitness of selected vs normal alleles: fitness weighted by frequency
        relative_advantage_selected_allele = ((1.0 + selection_coefficient)*frequency_in_current_generation*population_size)/(((1.0 + selection_coefficient)*frequency_in_current_generation*population_size)+((1.0)*(1.0 - frequency_in_current_generation)*population_size))
        #A binomial random number simulated population_size independent reproduction events.
        #The probability of each reproduction event applying to a selected allele is relative_advantage_selected_allele.
        number_after_reproduction_selected_allele = np.random.binomial(population_size, relative_advantage_selected_allele)
        number_after_reproduction_normal_allele = population_size - number_after_reproduction_selected_allele
        #Convert this to a new frequency.
        frequency_after_population_regulation = number_after_reproduction_selected_allele / float(population_size)
        #Record the new frequency in generation t + 1
        frequencies_per_generation.append(frequency_after_population_regulation)
    #(Convert to a numpy array for efficient manipulation)
    frequencies_per_generation = np.array(frequencies_per_generation)
    return frequencies_per_generation

def plot_compare_deterministic_stochastic_selection(initial_frequency = 0.01, population_size = 100, selection_coefficient = 0.0, time_steps = 100, stochastic_replicates = 100, show_stochastic = True, show_average = True, show_deterministic = True, show_conditioned_averages = False, plot_freq_lim = 'None', outfile = None):
    """
    Visualisation function for comparing the deterministic and stochastic models of selection.
    """
    deterministic_result = selection_deterministic(initial_frequency = initial_frequency, population_size = population_size, selection_coefficient = selection_coefficient, time_steps = time_steps)
    if stochastic_replicates > 0:
        stochastic_results = []
        for i in range(stochastic_replicates):
            stochastic_results.append(selection_stochastic_haploid(initial_frequency = initial_frequency, population_size = population_size, selection_coefficient = selection_coefficient, time_steps = time_steps))
        stochastic_results = np.array(stochastic_results)
        stochastic_results_average = np.average(stochastic_results, 0)
        
        stochastic_results_not_fixed = copy.copy(stochastic_results)
        stochastic_results_not_fixed[(stochastic_results_not_fixed == 0) + (stochastic_results_not_fixed == 1)] = np.nan
        stochastic_results_not_fixed_average = np.nanmean(stochastic_results_not_fixed, 0)

        stochastic_results_not_extinct = stochastic_results[(stochastic_results[::,-1] > 0)]
        stochastic_results_not_extinct_average = np.nanmean(stochastic_results_not_extinct, 0)

        fixed = np.sum(stochastic_results[::,-1] == 1)
        extinct = np.sum(stochastic_results[::,-1] == 0)
        polymorphic = stochastic_replicates - fixed - extinct
        #
        print("In %d generations: %.1f%% of stochastic runs fixed, %.1f%% still polymorphic, %.1f%% extinct" %(time_steps, (fixed / float(stochastic_replicates))*100.0, (polymorphic / float(stochastic_replicates))*100.0, (extinct / float(stochastic_replicates))*100.0))
        print("s is %.3f; selection dominates drift if s > %.3f (%s). Probability of fixation if initial frequency is 1/N = %.5f is approximately 2s = %.3f" %(selection_coefficient, 1.0 / population_size, selection_coefficient > 1.0 / population_size, 1.0 / population_size, 2 * selection_coefficient))

    #Plotting
    fig, ax = pl.subplots()
    if show_stochastic is True:
        for stochastic_result in stochastic_results:
            if stochastic_result[-1] == 1.0:
                line = ax.plot(stochastic_result, color = 'tab:orange', alpha = 0.1)
            else:
                line = ax.plot(stochastic_result, color = 'tab:blue', alpha = 0.1)
    
        ax.plot([], color = 'tab:blue', alpha = 0.1, label = 'Stochastic runs')
        ax.plot([], color = 'tab:red', alpha = 0.1, label = 'Stochastic runs, fixed derived')
        if show_average is True:
            ax.plot(stochastic_results_average, color = 'tab:green', linestyle = '--', label = 'Average frequency trajectory (all runs)')
        if show_conditioned_averages is True:
            ax.plot(stochastic_results_not_fixed_average, color = 'tab:green', linestyle = '-.', label = 'Average frequency trajectory (polymorphic runs)')
            ax.plot(stochastic_results_not_extinct_average, color = 'tab:green', linestyle = ':', label = 'Average frequency trajectory ($f_{%d}(A_{1}) > 0$)' %(time_steps))
    
    if show_deterministic is True:
        ax.plot(deterministic_result, color = 'tab:red', label = 'Predicted frequency (determinsitic)')
       
    ax.legend()
    ax.set_xlim(0, time_steps)
    ax.set_xlabel(r'Time (generations, $t$)')
    if plot_freq_lim != 'None':
        ax.set_ylim(plot_freq_lim[0], plot_freq_lim[1])
    ax.set_ylabel(r'Allele frequency')
    title = r'$N = %d$, $f_{0}(A_{1}) = %g$, $s = %g$' %(population_size, initial_frequency, selection_coefficient) if selection_coefficient > 0.0 else r'$N = %d$, $f_{0}(A_{1}) = %g$' %(population_size, initial_frequency)
    ax.set_title(title)
    if type(outfile) is type(None):
        pl.show()
    else:
        pl.savefig(fname = outfile)
    return None


parser = argparse.ArgumentParser(description='Some functions to explore the impact of neutral drift and selection using forward time simulation.')

parser.add_argument('--outfile', metavar='outfile', type=str, nargs='?', default = 'example_run.png',
                    help='system location of the output visualisation')
parser.add_argument('--model', '-m', dest='model', action = 'store', nargs='?', default = 'neutral',
                    help='is the model neutral or does it include selection?')

parser.add_argument('--initial_frequency', '-freq', dest='initial_frequency', type = float, action = 'store', default = 0.5,
                    help="what is the initial frequency of the simulated allele?")
parser.add_argument('--selection_coefficient', '-s', dest='selection_coefficient', type = float, action = 'store', default = 0.01,
                    help="what is the selection coefficient of the allele?")
parser.add_argument('--time_steps', '-num_steps', dest='time_steps', type = int, action = 'store', default = 200,
                    help="how long to run the simulation for?")
parser.add_argument('--population_size', '-N', dest='population_size', type = int, action = 'store', nargs = '?', default = '100',
                    help="population size (haploid).")
parser.add_argument('--num_finite_runs', '-num_f', dest='num_finite_runs', type = int, action = 'store', nargs = '?', default = '100',
                    help="number of finite population size runs.")
parser.add_argument('--run_infinite_size', '-inf', dest='run_infinite_size', action = 'store_true', default = False,
                    help="should an infinite population size simulation be run? Flag.")


args = parser.parse_args()

if args.model == 'neutral':
    args.selection_coefficient = 0.0

if args.outfile == "None":
    print("Please specify an output file.")

elif args.model not in ['selection', 'neutral']:
    print("Please use one of 'selection' or 'neutral' for the model.")

elif args.population_size < 1:
    print("Please give a population size value > 1.")

elif args.initial_frequency < 0 or args.initial_frequency > 1.0:
    print("The initial frequency must be between 0 and 1.")

elif args.model == 'neutral' and args.selection_coefficient != 0:
    print("Please set the model to selection if the selection_coefficient is not 0.")

elif args.time_steps < 1:
    print("Please set the number of time steps to >= 1")

else:
    plot_compare_deterministic_stochastic_selection(initial_frequency = args.initial_frequency, population_size = args.population_size, selection_coefficient = 0.0 if args.model == 'neutral' else args.selection_coefficient, time_steps = args.time_steps, stochastic_replicates = args.num_finite_runs, show_stochastic = True if args.num_finite_runs > 0 else False, show_average = True if args.num_finite_runs > 0 else False, show_deterministic = args.run_infinite_size, show_conditioned_averages = False, plot_freq_lim = 'None', outfile = args.outfile)    
    


