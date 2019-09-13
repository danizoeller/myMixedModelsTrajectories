# myMixedModelsTrajectories: Trajectory fitting using mixed models regression
Mixed effect model toolbox for the analysis of longitudinal data

This toolbox allows to fit models of different orders (from constant to cubic models) to data with repeated measurements. The purpose is to estimate developmental curves over age, in a mixed sample, where subjects were recorded at different ages and at multiple time points. Aside from determining the best model (no age-relationship, linear age-relationship, quadratic age-relationship or cubic age-relationship), the toolbox allows to estimate the difference in development between multiple groups (cf. *Note on the grouping information*).

Main steps follow the algorithm fist proposed in Mutlu et al., Neuroimage 2015:
  1. Fit models of increasing order to the data and select the best model based on the Bayesian Information Criterion
  2. Estimate p-values of group differences in intercept and shape of the curves for multiple groups
  3. Correct for multiple comparisons using the False Discovery Rate
  4. Plot resulting model parameters and fitted curves
  
Please cite the following papers when using this code:
- Mutlu, A.K., Schneider, M., Debbané, M., Badoud, D., Eliez, S., Schaer, M., 2013. Sex differences in thickness, and folding developments throughout the cortex. Neuroimage 82, 200–207. doi:http://dx.doi.org/10.1016/j.neuroimage.2013.05.076
- Mancini, V., Sandini, C., Padula, M.C., Zöller, D., Schneider, M., Schaer, M., Eliez, S., 2019. Positive psychotic symptoms are associated with divergent developmental trajectories of hippocampal volume during late adolescence in patients with 22q11DS. Mol. Psychiatry d. doi:10.1038/s41380-019-0443-z
    
Example applications, which used this toolbox:
- Mancini, V., Sandini, C., Padula, M.C., Zöller, D., Schneider, M., Schaer, M., Eliez, S., 2019. Positive psychotic symptoms are associated with divergent developmental trajectories of hippocampal volume during late adolescence in patients with 22q11DS. Mol. Psychiatry d. doi:10.1038/s41380-019-0443-z
- Maeder, J., Sandini, C., Zöller, D., Schneider, M., Bostelmann, M., Pouillard, V., Caroni, P., Kliegel, M., Maeder, J., Sandini, C., Zöller, D., Schneider, M., Pouillard, V., Caroni, P., Kliegel, M., Long-term, S.E., 2019. Long-term verbal memory deficit and associated hippocampal alterations in 22q11 . 2 deletion syndrome. Child Neuropsychol. 00, 1–23. doi:10.1080/09297049.2019.1657392
- Franchini, M., Zöller, D., Gentaz, E., Glaser, B., De Wilde, H.W., Kojovic, N., Eliez, S., Schaer, M., 2018. Early adaptive functioning trajectories in preschoolers with autism spectrum disorders. J. Pediatr. Psychol. 43. doi:10.1093/jpepsy/jsy024
- Maeder, J., Schneider, M., Bostelmann, M., Debbané, M., Glaser, B., Menghetti, S., Schaer, M., Eliez, S., 2016. Developmental trajectories of executive functions in 22q11.2 deletion syndrome. J. Neurodev. Disord. 8, 10. doi:10.1186/s11689-016-9141-1

# Installation & Use
1. Clone or download this directory
2. Have a look at the three example main scripts, follow the comments and adapt the script for your purpose

### Note on the grouping information
The toolbox allows works on a single group, or multiple groups. Resulting p-values indicate whether a model with multiple groups is significantly better than assuming no grouping information. If you want to get a p-value for the developmental difference between two out of many groups, you need to run the algorithm separately for subjects of only these two groups.

### Correction for multiple comparisons
By default, the example scripts correct for multiple comparisons using the False Discovery Rate. If you do not want to correct within the toolbox, please remove the respective line from the script.

# Credits
This code was written by Daniela Zöller, based on original scripts by Kadir Mutlu. It uses the fdr_bh (https://ch.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh) function to correct for multiple comparisons and cbrewer (https://ch.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab), which includes color specifications and designs developed by Cynthia Brewer (http://colorbrewer.org/), and the function shadedErrorBar (https://ch.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar).

# Versions
v2.0 - 11.9.2019 Daniela Zöller
- update for publication on github

v1.4 - 11.4.2018 Daniela Zöller
- some bug fixes in the plotting function to make it work with only one group

v1.3 - 20.2.2018 Daniela Zöller
- added centering of covariates
- new function to correct for multiple comparisons (fdr_correct, example use in main_ma_example.m)
- added option to change the size of the figures
- changes default text size in startup to 12
- included the option to estimate only a GLM instead of mixed models (by using opts.mType=‘glm’)

v1.2  - 7.9.2016 Daniela Zöller
corrected bug for inclusion of multiple covariates

v1.0 - 7.9.2016 Daniela Zöller


# Known issues
- All data in the excel file should be numeric (also the grouping information)



