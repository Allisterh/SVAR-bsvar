# bsvar 

Bayesian estimation of non-Gaussian structural vector autoregressions (SVARs) via Stan. More specifically, this R-package aims at making Bayesian estimation of statistically identified non-Gaussian SVARs as easy as possible by providing an easy-to-use interface to [cmdstanr](https://mc-stan.org/cmdstanr/) along with various ready-made tools for the analysis of SVARs. 

## Getting Started

### Prerequisites 

You must have [cmdstanr](https://mc-stan.org/cmdstanr/) installed (and properly set-up) on your computer before installing the package. However...

**IMPORTANT: For now, cmdstan versions 2.32.1 and up do not work with the package, so you must use older release (but not older than 2.29.0)** 

Thus, in setting up [cmdstanr](https://mc-stan.org/cmdstanr/) you may follow the instructions in [here](https://mc-stan.org/cmdstanr/articles/cmdstanr.html), but with the exception of **NOT** installing the latest version of cmdstan. If you're getting confused and just came here for an easy-to-use R-package, don't worry, just follow the instructions below and soon all these unpleasant technicalities will be over and done with. 


First, if you don't have [cmdstanr](https://mc-stan.org/cmdstanr/) R-package already installed, open a fresh R session and install [cmdstanr](https://mc-stan.org/cmdstanr/) by running the following command

```
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```


Second, as required by [cmdstanr](https://mc-stan.org/cmdstanr/), make sure you have a suitable C++ toolchain installed on your computer. What this means in practice? 

**On Windows:** This means that you have to have Rtools installed on your computer. This is not an R-package and you may download Rtools installer for instance from [here](https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html).

**On Mac:** This means that you must probably have at least Xcode command line tools installed on your computer. If you are not sure, the prompt in the next step should tell you if you need to do something.

**On Linux:** You probably already have what you need.


Third, make sure your toolchain is set up properly by running the following commands:

```
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE)
```
This should tell you if your toolchain is set up properly and it might also prompt you to accept any automatic fixes if this isn't the case. If you get stuck here, see [this](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) or [this](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for more thorough instructions (or you can just carry on and hope for the best).


Fourth, install a suitable [release](https://github.com/stan-dev/cmdstan/releases) (< 2.32.1) of cmdstan on your computer (not the [cmdstanr](https://mc-stan.org/cmdstanr/) R-package we have already installed!). For now it is recommended to use version 2.29.2 since this package has been most thoroughly tested with that version. You may install this version of cmdstan on your computer by running the following command (this might take a few minutes):

```
install_cmdstan(version = "2.29.2")
```

**NOTE: If the above command fails when using Rstudio, try running the command without Rstudio. You can get back to Rstudio after succesfully installing a suitable version of cmdstan.**

Now [cmdstanr](https://mc-stan.org/cmdstanr/) should be all set up. Even if you have an installation of the wrong cmdstan version (<2.29.0 or >2.32.0) on your computer (e.g. you have installed it earlier or by accident), `bsvar()` should find the right version and use that. However, if problems do persist, you may navigate to the location of cmdstan installations as given by `cmdstan_path()` and manually remove the folder(s) cmdstan-2.32.1, cmdstan-2.32.2... (or something similar) and see if that helps.


### Installing the package 

When you have [cmdstanr](https://mc-stan.org/cmdstanr/) all set up you are ready to install the package. Simply run the following commands: 

```
if(!("devtools" %in% installed.packages())) install.packages("devtools")
devtools::install_github("jetroant/bsvar")
```

The installation of the package itself should be a breeze and only take a few seconds, however if there are a lot of dependencies to be installed (that is, if you do not already have other R-packages that are needed for this package to work installed on your computer) this might take a few minutes.

### Test drive 

Everything should now be ready for you to use [bsvar](https://github.com/jetroant/bsvar). First time you call `bsvar()`, a cmdstan program is compiled, which might take a few minutes. This however only needs to be done once, after which such speed bumps cease to exist. To test whether everything works as intended and to get that first time model compilation out of your way, you may run the following commands: 

```
library(bsvar) # Loads the package
fit <- bsvar(y = matrix(rnorm(200), ncol = 2), lags = 2) # Estimates a bivariate model with two lags using random Gaussian data
prnt(fit) # Summarizes the parameter estimates
irfs(fit) # Plots impulse response functions
dist_plot(fit) # Depicts the estimated shock distribution parameters
```

These are by no means all the available tools included in the package, but hopefully give some idea of what can be done with the package. More elaborate examples will most probably be added here in the near future as the documentation of the package progresses. 

### Marginal likelihood estimation 

Finally, a few words regarding the marginal likelihood estimation via the excellent [bridgesampling](https://github.com/quentingronau/bridgesampling) R-package. `bridgesampling::bridge_sampler` can be called by runnning the following command: 

```
marginal_likelihood(fit)
```
`marginal_likelihood()` is a wrapper for `bridgesampling::bridge_sampler()` that takes care of some pesky practicalities. Providing the `stanfit` object outputted by `bsvar()` directly to `bridgesampling::bridge_sampler()` will not work. This is mainly because (for now) [bridgesampling](https://github.com/quentingronau/bridgesampling) only supports Stan models built with [rstan](https://mc-stan.org/users/interfaces/rstan), not with [cmdstanr](https://mc-stan.org/cmdstanr/). To this end, `marginal_likelihood()` recompiles an [rstan](https://mc-stan.org/users/interfaces/rstan) program based on the same code as the original [cmdstanr](https://mc-stan.org/cmdstanr/) program used for sampling from the posterior. For now, this needs to be once every session (annoying, I know!), so calling `marginal_likelihood()` for the first time in a session may take a few minutes as opposed to seconds later on. Also, very importantly, calling `marginal_likelihood()` requires an installation of `rstan 2.26.0` or higher. You can check the version of your installation by running `rstan::stan_version()`. Check the instructions for updating your [rstan](https://mc-stan.org/users/interfaces/rstan) installation [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).














