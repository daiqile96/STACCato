# STACCato
Identifying Condition-Related Cell-Cell Communication Events Using Supervised Tensor Analysis.

![STACCato Framework](Figure1.png)

- (1) perform tensor-based regression to estimate the effects of the biological condition of interest on CCC events while adjusting for other covariates; 
- (2) use the bootstrapping resampling method to assess the significance level of condition-related CCC events; 
- (3) conduct downstream analyses, including comparing significant CCC events across cell types and visualizing CCC events that are significantly associated with the condition of interest.

## Tutorials

- Prepare inputs for STACCato:
  - [use LIANA to build communication score tensor](Examples/0_construct_communication_score_tensor.ipynb)

- Run STACCato:
  - [Apply STACCato to identify condition-related CCC events](https://htmlpreview.github.io/?https://raw.githubusercontent.com/daiqile96/STACCato/main/Examples/STACCato.html)
  



