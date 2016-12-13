# CFODD: Cloud Microphysical Process Metrics
This is a model diagnostic metric that evaluates warm rain formation microphysical processes against satellite observations.
The methodology combines satellite cloud observables (e.g. radar reflectivity, cloud optical depth and cloud-top effective particle radius) 
to construct a particular statistics that "fingerprint" signatures of the warm rain formation process in the form of a 
contoured frequency diagram of radar reflectivity by in-cloud optical depth. 
The statistics are typically classified according to different ranges of the cloud-top particle size to depict 
how vertical microphysical structures tend to transition from non-precipitating clouds to precipitating clouds 
as a fairly monotonic function of the cloud-top particle size. 
The corresponding statistics are constructed from output of a given model for comparisons with the satellite-based statistics 
to expose how the model realistically represents the warm rain formation process in a statistical sense. 
The methodology is applicable to both global climate models and cloud-resolving models.

# Usage
Fortran 90 codes for constructing the metric are available at the directory "code/fortran" in this repository. 
There are two codes: one for satellite-based statistics and the other for model-based statistics using GFDL-CM3 sample output.
The satellite and model data required to run the codes are available upon request to Kentaroh Suzuki at ksuzuki@aori.u-tokyo.ac.jp.
The data will be shared through Google drive, to which the user will be given a permission to access.
The sample results obtained from the code and data are included in the directory "images".

# References
Suzuki, K., G. Stephens, A. Bodas-Salcedo, M. Wang, J.-C. Golaz, T. Yokohata, and T. Koshiro, 2015: Evaluation of the warm rain formation process in global models with satellite observations. J. Atmos. Sci., 72, 3996-4014, doi:10.1175/JAS-D-14-0251.1.

Suzuki, K., J.-C. Golaz, and G. L. Stephens, 2013: Evaluating cloud tuning in a climate model with satellite observations. Geophys. Res. Lett., 40, 4464-4468, doi:10.1002/grl.50874.

Suzuki, K., G. L. Stephens, S. C. van den Heever, and T. Y. Nakajima, 2011: Diagnosis of the warm rain process in cloud-resolving models using joint CloudSat and MODIS observations. J. Atmos. Sci., 68, 484-503.
