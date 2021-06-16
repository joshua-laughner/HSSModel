# HSSModel.jl Documentation

HSSModel implements a basic NOx-HOx steady state model that
can be used to compute NOx lifetime, ozone production efficency,
and other quantities for theoretical conditions.

## Installation

Currently, this package is only available to install from GitHub. From the Julia
REPL:

```
julia> ]add https://github.com/joshua-laughner/HSSModel.git
```

Since this will track the `master` branch, you may also wish to pin this package 
so that it does not update until you want it to:

```
julia> ]pin HSSModel
```

## Basic Walkthrough

Importing `HSSModel` exposes the main driver function, `hox_ss_solver` from
the `SteadyStateSolvers` module. This function requires 5 inputs:

1. NO concentration (in molec. cm``^{-3}``)
1. NO2 concentration (in molec. cm``^{-3}``)
1. HOx production rate (in molec. cm``^{-3}`` s``^{-1}``)
1. Total VOC reactivity (in s``^{-1}``)
1. The NO + RO2 branching ratio, ``\alpha`` (unitless).

For this example, we'll be using values from Table 4 of 
[Murphy et al. 2006](https://acp.copernicus.org/preprints/6/11971/2006/acpd-6-11971-2006.pdf)
for P(HOx), VOC reactivity, and ``\alpha``, plus we will assume
5 ppb of NOx. Since our inputs need to be in number density,
rather than mixing ratio, we'll need the number density of
air to convert.

```@example 1
# gas constant in Pa cm^3 K^-1 molec^-1
R = 8.314 * 1e6 / 6.022e23

# number density of air at STP
nair = 101325 / (R * 298)

# Convert from NOx in ppb to NO and NO2 in molec. cm^-3
# assuming a 4:1 NO2:NO ratio
nox_ppb = 5;
no_nd = nox_ppb*1e-9 * 0.2 * nair;
no2_nd = nox_ppb*1e-9 * 0.8 * nair;

# Other default values from Murphy et al. 2006
phox = 2.5e-9 * nair / 3600;  # convert ppb hr^-1 -> molec. cm^-3 s^-1
vocr = 5.8;
alpha = 0.04;

nothing # hide
```

With these conditions, we call the steady state solver, which
will compute the steady state concentrations of OH, HO2, and RO2
given these fixed values.

```@example 1
using HSSModel;
result = hox_ss_solver(no_nd, no2_nd, phox, vocr, alpha);

println("OH = $(result.oh), HO2 = $(result.ho2), RO2 = $(result.ro2)");
```

This result structure contains information about the
steady state concentrations, as well as the model configuration.
Generally, most information you would need to plot the model
results is contained in this structure. 

These result structure can also be use to calculate secondary
quantities. The `DerivedQuantities` module exports two functions
to do so: `nox_lifetime` and `ozone_production_efficiency`. Calling
these with your results produces:

```@example 1
nox_lifetime(result)
```

This returns a dictionary with three different lifetimes, in hours:

* `total` is the overall NOx lifetime.
* `hno3` is the lifetime with respect to loss to HNO3
* `ans` is the lifetime with respect to loss to alkyl nitrates

```@example 1
ozone_prod_efficiency(result)
```

Unlike `nox_lifetime`, this returns a single value. This represents
the ratio of ozone production to NOx loss.

## Running an ensemble

The steady state model function works well with Julia's built in 
broadcasting behavior. For a simple example, let's look at running
a range of NOx concentrations:

```@example 1
# Compute NO and NO2 number densities for logarithmically-spaced
# NOx mixing ratios between 0.01 and 10 ppb, still assuming a
# 4:1 NO2:NO ratio.
nox_nd = exp10.(range(-2., stop=1., length=20)) .* 1e-9 .* nair;
no2_nd = 0.8 .* nox_nd;
no_nd = 0.2 .* nox_nd;

# Note the broadcasting operator (the period following `hox_ss_solver`)
results = hox_ss_solver.(no_nd, no2_nd, phox, vocr, alpha);
size(results)
```

If you want to test an ensemble with multiple variables, you can 
use broadcasting in multiple dimensions. In this example, we'll
use the same NOx vector (a 1D, 75 element vector) along with a
vector of VOC reactivity values. By making the VOCR values a
row (1-by-``N``) vector, the broadcasting will result in a 
2D array of results.

```@example 1
# Create a 1-by-10 vector of VOC reactivities
vocr_vec = collect(transpose(range(1., stop=10., length=10)));
results = hox_ss_solver.(no_nd, no2_nd, phox, vocr_vec, alpha);
size(results)
```

A useful trick for quickly extracting values from these results
arrays (for plotting or downstread analysis) is to broadcast the
`getfield` function:

```@example 1
getfield.(results, :oh)
```
