<!-- README.md is generated from README.Rmd. Please edit that file -->

# metid <img src="man/figures/logo.png" align="right" alt="" width="120" />

[![](https://www.r-pkg.org/badges/version/massdatabase?color=green)](https://cran.r-project.org/package=massdatabase)
[![](https://img.shields.io/github/languages/code-size/tidymass/massdatabase.svg)](https://github.com/tidymass/massdatabase)
[![Dependencies](https://tinyverse.netlify.com/badge/massdatabase)](https://cran.r-project.org/package=massdatabase)
[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)


`massdatabase` is a part of [tidymass](https://tidymass.github.io/tidymass/)

-------

# About

`massdatabase` is a R package which is used for metabolite identification based
on in-house database and public database based on accurate mass (m/z),
retention time (RT) and/or MS2 spectra.

<img src="man/figures/Figure_1.png" align="middle" alt="" width = "100%"/>

# Installation

You can install `massdatabase` from [GitLab](https://gitlab.com/jaspershen/massdatabase)

``` r
if(!require(remotes)){
install.packages("remotes")
}
remotes::install_gitlab("jaspershen/massdatabase")
```
or [Github](https://github.com/tidymass/massdatabase)

``` r
remotes::install_github("tidymass/massdatabase")
```

`massdatabase` is a part of `tidymass`, so you can also install it by installing [`tidymass`](https://www.tidymass.org/).

# Usage

Please see the `Help documents` page to get the instruction of `massdatabase`.

## Need help?

If you have any questions about `massdatabase`, please don’t hesitate to email me (<shenxt@stanford.edu>).

<i class="fa fa-weixin"></i>
[shenzutao1990](https://www.shenxt.info/files/wechat_QR.jpg)

<i class="fa fa-envelope"></i> <shenxt@stanford.edu>

<i class="fa fa-twitter"></i>
[Twitter](https://twitter.com/JasperShen1990)

<i class="fa fa-map-marker-alt"></i> 
[M339, Alway Buidling, Cooper Lane, Palo Alto, CA 94304](https://www.google.com/maps/place/Alway+Building/@37.4322345,-122.1770883,17z/data=!3m1!4b1!4m5!3m4!1s0x808fa4d335c3be37:0x9057931f3b312c29!8m2!3d37.4322345!4d-122.1748996)

# Citation

If you use `massdatabase` in your publications, please cite this paper:

TidyMass: An Object-oriented Reproducible Analysis Framework for LC-MS Data.

Xiaotao Shen, Hong Yan, Chuchu Wang, Peng Gao, Caroline H. Johnson, Michael P. Snyder.

[Web Link](https://www.biorxiv.org/content/10.1101/2022.03.15.484499v1).

Thanks very much!
