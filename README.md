# CFODD: Cloud Microphysical Process Metrics
CFODD (Contoured Frequency by Optical Depth Diagram) is a model diagnostic metric that evaluates warm rain formation microphysical processes against satellite observations. The methodology combines satellite cloud observables (e.g. radar reflectivity, cloud optical depth and cloud-top effective particle radius) to construct a particular statistics that "fingerprint" signatures of the warm rain formation process in the form of a contoured frequency diagram of radar reflectivity by in-cloud optical depth. The statistics are typically classified according to different ranges of the cloud-top particle size to depict how vertical microphysical structures tend to transition from non-precipitating clouds to precipitating clouds as a fairly monotonic function of the cloud-top particle size. The corresponding statistics are constructed from output of a given model for comparisons with the satellite-based statistics to expose how the model realistically represents the warm rain formation process in a statistical sense. The methodology is applicable to both global climate models and cloud-resolving models.

# Usage
Fortran 90 codes for constructing the metric are available at the directory "code/fortran" in this repository. There are two codes: one for satellite-based statistics using CloudSat and MODIS data product and the other for model-based statistics using GFDL-CM3 sample output. The satellite and model data required to run the codes are available upon request to Kentaroh Suzuki at ksuzuki@aori.u-tokyo.ac.jp. The data will be shared through Google Drive, to which the user will be given a permission to access. The sample results obtained from the code and data comparing the satellite-based and model-derived statistics are included in the directory "images".

# Related publications (selected)
[Suzuki, K., G. Stephens, A. Bodas-Salcedo, M. Wang, J.-C. Golaz, T. Yokohata, and T. Koshiro, 2015](http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-14-0265.1): Evaluation of the warm rain formation process in global models with satellite observations. J. Atmos. Sci., 72, 3996-4014, doi:10.1175/JAS-D-14-0265.1.

[Suzuki, K., J.-C. Golaz, and G. L. Stephens, 2013](http://onlinelibrary.wiley.com/doi/10.1002/grl.50874/abstract): Evaluating cloud tuning in a climate model with satellite observations. Geophys. Res. Lett., 40, 4464-4468, doi:10.1002/grl.50874.

[Suzuki, K., G. L. Stephens, S. C. van den Heever, and T. Y. Nakajima, 2011](http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-10-05026.1): Diagnosis of the warm rain process in cloud-resolving models using joint CloudSat and MODIS observations. J. Atmos. Sci., 68, 484-503.

# Input
| Frequency | Duration | Variables | Dimension | CMOR labels | Unit | File Format |
| --------- | -------- | --------- | --------- | ----------- | ---- | ----------- |
| 6 hourly  | 3 months | Cloud optical thickness (liquid) | 2D | tau | Unitless | nc |
|           |          | Cloud-top effective radius (liquid) | 2D | reffclwtop | micron | nc |
|           |          | Cloud-top temperature | 2D | N/A | K | nc |
|           |          | Liquid water path | 2D | N/A | kg/m2 | nc |
|           |          | In-cloud optical depth (liquid, St) | 3D | dtaus | Unitless | nc |
|           |          | Radar Reflectivity | 3D+Subcolumn | N/A | dBZ | nc |
|           |          | Fracout (Subcolumn scence index) | 3D+Subcolumn | N/A | Unitless | nc |

*CMOR labels denoted "N/A" indicate that the variable is not available in current archive of CMIP.

# Output
The output is the occurrence frequency of radar reflectivity normalized at each in-cloud optical depth in the form of the contoured frequency by optical depth diagram. The sample results obtained from GFDL CM3 model output and satellite observations are provided as "cfodd_gfdl_sg_std_jan_4class.txt" and "cfodd_r21_5class_djf.txt" in the "code/fortran" directory. These are visualized to be displayed as "cfodd_sat_gfdl.png" in the "image" directory.

# Program
The diagnostic code:

code/fortran/analysis_cfodd_gfdl.f90 (Model: GFDL CM3)

code/fortran/analysis_cfodd_sat.f90 (Satellite)

The visualization script:

code/fortran/draw_cfodd.exe for Gnuplot

# Data availability
The sample data from satellite and model required to create this diagnostic are available upon request to Kentaroh Suzuki at ksuzuki@aori.u-tokyo.ac.jp.
