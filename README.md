# LOBD
Codes that implement the Linear Orthogonal Basis Decomposition (LOBD) in MATLAB. The codes use the Tensorlab package which can be installed from https://www.tensorlab.net. 

The main function is contained in LOBD.m. See the examples folder for examples of using the function and the paper for more information about what is calculated. 

The required inputs are a cell of input data, all the same size, which must have length greater than or equal to two and the desired rank of the decompositon. 

The following are optional parameters that adapt the basic LOBD factorization, 
  
  * 'initialvalues' - the returned lobd.variables struct that can be used to warm start the optimization
  
  * 'tpoints' - the time points that the data are sampled at (only required if expt is true)
  
  * 'maxiters' - maximum number of iterations the optimization runs for (default: 100)
  
  * 'nonneg' - whether t is enforced to be nonnegative (default: false)
  
  * 'orthogonal' - whether X is enforced to have orthogonal columns (default: true)
  
  * 'useminf' - if true minf optimizer in Tensorlab is used rather nls (default: true)
  
  * 'cgiters' - number of CG iterations used by Tensorlab (default: 500)
  
  * 'sameT' - whether T is enforced to be the same in all factorizations (default: true)
  
  * 'tvreg' - whether a derivative regularization term on X and T is included (default: false)
  
  * 'projectcoeffs' - whether the cofficients are used as a fit parameter or enforced to be the projection of the basis against the initial condition.
  
  * 'tvweight' - relative weight of the regularization term (default: 1). Only used in 'tvreg' is true.
  
  * 'expt' - whether T is enforced to be exponential (default: false). 
