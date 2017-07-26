# ODonnelletal_2017_imbalances
Model code for O'Donnell et al, 2017.

Study on multidimensional imbalances in neural circuit properties in brain disorders.

Files include:

Brian circuit model of L2/3 somatosensory cortex
- layer23model.py; Brian simulation that sets up single L2/3 circuit model and simulates one pattern of input activity, plots summary results.
- layer23model_varyparams.py; Brian simulation of same circuit model that averages over many runs, input patterns, and input sparsities, then varies each parameter by 20%.
