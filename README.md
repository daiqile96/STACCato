# STACCato
Supervised Tensor Analysis tool for studying Cell-cell Communication using scRNA-seq data across multiple samples and conditions

![STACCato Framework](Figure1.pdf)

- (1) construct 4-dimensional communication score tensor using multi-sample multi-condition scRNA-seq data; 
- (2) perform supervised tensor decomposition to estimate the effects of conditions on CCC events and infer activity patterns of cell types; 
- (3) use bootstrapping resampling to assess the significance level of the estimated effects; 
- (4) conduct downstream analyses including comparing significant CCC events across cell types and identify pathways significantly associated with conditions. 
