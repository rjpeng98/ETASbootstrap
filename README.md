<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ETASbootstrap)](https://cran.r-project.org/package=ETASbootstrap)
[![CRAN_Download_Count](http://cranlogs.r-pkg.org/badges/ETASbootstrap)](https://cran.r-project.org/package=ETASbootstrap)
[![R Build Status](https://github.com/rjpeng98/ETASbootstrap/workflows/R-CMD-check/badge.svg)](https://github.com/rjpeng98/ETASbootstrap/actions)
<!-- badges: end -->

# ETASbootstrap

For detailed documentation, please visit [the user manual](https://cran.r-project.org/web/packages/ETASbootstrap/ETASbootstrap.pdf).

*ETASbootstrap* is developed to compute the bootstrap confidence intervals based on empirical quantiles for the parameters of the 2-D spatial and temporal 'ETAS' model. 
See Dutilleul, P., Genest, C., Peng, R., 2024. [Bootstrapping for parameter uncertainty in the space–time epidemic-type aftershock sequence model](https://academic.oup.com/gji/article/236/3/1601/7511107).
Geophysical Journal International 236, 1601–1608, for details of the algorithm.
This version improves on Version 0.1.0 of the package by enabling the study space window to be polygonal rather than merely rectangular. 

## Installation

To install the package from [CRAN](https://cran.r-project.org/package=ETASbootstrap), run the following in R:
```R
install.packages("ETASbootstrap")
```

You can also install the current version of the package on GitHub by running:
```R
require(remotes)
install_github('rjpeng98/ETASbootstrap')
```

If [remotes](https://github.com/r-lib/remotes) is not installed, you should first run:

```R
install.packages('remotes')
```

## Test

A quick-test file has been included in the repository and can be run directly after downloading.
The online checking results for the package are available at [https://mirror.its.dal.ca/cran/web/checks/check_results_ETASbootstrap.html](https://mirror.its.dal.ca/cran/web/checks/check_results_ETASbootstrap.html).

## License

This project is licensed under

- [MIT license](https://opensource.org/licenses/MIT) ([`LICENSE-MIT`](https://github.com/rjpeng98/ETASbootstrap/blob/main/LICENSE)).


