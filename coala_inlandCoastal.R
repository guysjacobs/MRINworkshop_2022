#Code adapted from work by Evie Carter (2022)

library(coala)
library(ggplot2)

model1=coal_model(sample_size=c(0,417,427) ,loci_number=100, loci_length=16000)+ 
  
  #Sampling 417 and 427 men (one X copy) from populations 2 and 6
  
  feat_mutation(rate=1.25e-8*4*10000*16000) +
  feat_recombination(rate=1.0e-8*4*10000*16000)+ #Rate=locus_length * per_generation_per_base_pair_mutation_rate * 4 * Ne
  
  #Population sizes *0.75 as evolution is on the X chromosome (3/4 of Ne vs autosomes)                                                                
  
  feat_size_change(((5000)/10000)*0.75, p="all", time=0) +  #Set village pop sizes
  feat_size_change((800/10000)*0.75, p=1, time=0) +                #Pop. 1: 'Melanesian' ancestry
  feat_migration(rate=0.015*4*10000, pop_from=2, pop_to=3, time=0) +                 #Set migration rate
  feat_migration(rate=0.015*4*10000, pop_from=3, pop_to=2, time=0) +                 #Set migration rate
  
  feat_pop_merge(pop_source=3, pop_target=2, time=161/40000) +      #merge the coastal and inland village into pop 2 at 160 gens (~4000 years)
  feat_size_change(((1425/10000))*0.75, population=2, time=161/40000) + #Change pop 2 size at 160 gens 
  
  feat_migration(rate=0.26*4*10000, pop_from=1, pop_to=2, time=162/40000) +  #pop 2 receives 26% from 'Melanesian' pop 1 in a pulse over one generation
  feat_migration(rate=0.0, pop_from=1, pop_to=2, time = 163/40000) +            

  feat_size_change((2050/10000)*0.75, population=2, time=164/40000) +     #pop 2 population size becomes 2050
  
  feat_pop_merge(pop_source=2, pop_target=1, time=2000/40000) +                #Merge pop 2 into 1 at 2000 gens (50000 years) 
  feat_size_change((10500/10000)*0.75, population=1, time=2001/40000) +                #ancestral pop. size set to 10500
  sumstat_jsfs("jsfs_2_3", populations=c(2,3), per_locus=FALSE) +
  sumstat_nucleotide_div() +
  sumstat_file('./output/')

model1

stats=simulate(model1)

##Diversity

#Average number of segregating sites
mean(stats$pi)
ggplot() + aes(stats$pi) + geom_histogram(binwidth = 2) + geom_vline(aes(xintercept=mean(stats$pi)), color = 'darkred', size = 1, linetype ='dashed')


##Testing selection.

#We'll look for SNPs with alleles between 11 and 22 (ie 16 +- 5) freqeuncy in the inland population.
#And ask how high frequency these are in the coastal population.

freq_slice=as.data.frame(stats$jsfs_2_3[12:23,])
#ncol(freq_slice) #428 - coastal
sumofcols=colSums(freq_slice) #sum of columns, ie adding up obesrvations over all the simulations

#To use the ggplot histogram function we convert thesethe output into observations. Ignore this.
observations = vector("list", sum(sumofcols))
vec <- 1:length(sumofcols)
idx = 1
for (x in seq_along(vec)) {
  if (sumofcols[x] > 0) {
    for (y in seq_along(1:sumofcols[x])) {
      observations[idx] = x - 1
      idx = idx + 1
    }
  }
}
observations_df = do.call(rbind.data.frame, observations)
colnames(observations_df) = c("obs_freq")

#P-val. This is calculated as the proportion of SNPs that have higher freq difference than the G6PDd observation.

sum(sumofcols[46:428]) / sum(sumofcols) #col 46 will correspond to variants at freq. of 45
#proportion of overall SNPs [(sum(sumofcols))] that are above a freq. of 45 
#About 0.14 (for this run) - not significant, unable to detect non-neutrality

ggplot(data=observations_df, aes(obs_freq)) + geom_histogram(binwidth = 8) + geom_vline(aes(xintercept=45), color = 'darkred', size = 1, linetype ='dashed')
