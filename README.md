# SCEM

R scripts to perform the Splitting-Coalescence-Estimation Method to model birth seasonality in studies of herd animals.

## Authors

[Hannah Chazin](http://columbia.academia.edu/HannahChazin), [Soudeep Deb](http://soudeepd.github.io/), [Joshua Falk](http://home.uchicago.edu/~jsfalk/), Arun Srinivasan

## Methods

We introduce improved methods for statistically assessing birth seasonality and intra-annual variation in &delta;<sup>18</sup>O from faunal tooth enamel. 

The first method estimates input parameters for use with a previously-developed parametric approach (Tornero et al., 2013). The relevant code for this approach is `makeFits.R`, while `makeFitsWrong.R` is the code to implement the same method but with given initial conditions for two parameters. The latter can be used to show the disadvantage of the existing approach.

The second method we propose is a new idea that uses a nonparametric clustering procedure to group individuals with similar time series data and estimate birth seasonality based on the clusters. This method is more efficient across different scenarios, especially when less of the tooth row is preserved. The new approach offers a high level of statistical rigor and flexibility in dealing with the time series data produced through intra-individual sampling in isotopic analysis. One can use the function `SCEM.R` to implement this method. 

Example of implementing the above methods for real data can be found in `Usage-Example.R`. Other functions in this repository are used internally in the above-mentioned functions. 


## Contact

For any inquiries or questions related to this, please open an issue in this repository. You can also contact us at [hc2986@columbia.edu](hc2986@columbia.edu) or [sdeb@uchicago.edu](sdeb@uchicago.edu).
