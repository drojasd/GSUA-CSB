## Summary
This submission contains Matlab implementation of an iterative **multivariate outlier detection algorithm** 
described in [Hadi (1992)] [[1]]. In addition to flagging potential outliers, the main function
`DetectMultVarOutliers.m` also outputs robust estimates of the mean and covariance that it computes 
during execution. 

Deviating slightly from [Hadi (1992)], `DetectMultVarOutliers.m` initializes the sample mean with the [geometric median] 
of the dataset, instead of the coordinate-wise median. `GeometricMedian.m` is the function used compute this 
robust statistic; via the Weiszfeld's algorithm [[2]]. Note that this auxiliary function can be used on its own in 
any application that requires robust estimation of central tendency of multivariate data corrupted by sampling 
errors and/or noise.

For a quick demo on how to use the main function, see source code for `outliers_demo.m`	or simply enter `outliers_demo` 
into Matlab command window.

## References
[**[1]**] Hadi, A.S., 1992. Identifying multiple outliers in multivariate data. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 54(3), pp. 761-771.  

[**[2]**] Weiszfeld, E., 1937. Sur le point par lequel la somme des distances den points donnés est minimum. Tohoku Mathematics Journal, Vol. 43, pp. 355–386.

## License
[MIT] © 2019 Anton Semechko 
a.semechko@gmail.com

[Hadi (1992)]: https://www.researchgate.net/profile/Ali_Hadi/publication/243777821_Identifying_Multiple_Outliers_in_Multivariate_Data/links/5406dda50cf2c48563b2732e.pdf
[1]: https://www.researchgate.net/profile/Ali_Hadi/publication/243777821_Identifying_Multiple_Outliers_in_Multivariate_Data/links/5406dda50cf2c48563b2732e.pdf
[geometric median]: http://en.wikipedia.org/wiki/Geometric_median
[2]: http://en.wikipedia.org/wiki/Geometric_median 
[source code]: https://github.com/AntonSemechko/Multivariate-Outliers/blob/master/outliers_demo.m

[MIT]: https://github.com/AntonSemechko/Multivariate-Outliers/blob/master/LICENSE.md
