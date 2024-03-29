module SteadyStateSolvers

using NLsolve;
import ..Rates;

export nonlin_nox_analytic_model,
       hox_ss_solver,
       SteadyStateOptions;


struct SteadyStateOptions
    T::Real
    M::Real
    rh::Real

    k_RO2NO::Real
    k_RO2HO2::Real
    k_RO2RO2::Real

    k4::Real
    k2eff::Real
    k5eff::Real
end

"""
    SteadyStateOptions(;T::Real=298, M::Real=2e19, rh::Real=0.01, 
                        k_RO2NO::Real=8e-12, k_RO2HO2::Real=8e-12, k_RO2RO2::Real=6.8e-14,
                        k4::Real=1.1e-11, k2eff::Real=8e-12, k5eff::Real=5e-12)

Construct a `SteadyStateOptions` instance, that is a structure that controls the
options for the HOx steady state solver. All parameters are keyword arguments. Available
parameters and their default values are:

* `T` (default: 298) - air temperature in Kelvin
* `M` (default: ``2 × 10^{19}``) - number density of air in molec. cm``^{-3}``.
* `rh` (default: 0.01) - mole fraction of H2O in the air, i.e. `rh * M` 
  is the number density of H2O used in the model.
* `k_RO2NO` (default: ``8 × 10^{-12}``) - fixed rate constant for the NO + RO2 → NO2 + RO 
  reaction. Since this needs to be an effective rate for the mixing of RO2 radicals, it is just
  used as a fixed value and doesn't change with T and M. Units of molec``^{-1}`` cm``^{3}`` s``^{-1}``.
* `k_RO2HO2` (default: ``8 × 10^{-12}``) - fixed rate constant for the RO2 + HO2 loss pathway.
  Like `k_RO2NO`, this doesn't vary with T and M. Units of molec``^{-1}`` cm``^{3}`` s``^{-1}``.
* `k_RO2RO2` (default: ``6.8 × 10^{-14}``) - fixed rate constant for the RO2 + RO2 loss pathway.
  Like `k_RO2NO`, this doesn't vary with T and M. Units of molec``^{-1}`` cm``^{3}`` s``^{-1}``.
* `k4` (default: ``1.1 × 10^{-11}``) - rate constant for NO2 + OH → HNO3 reaction *only* in the 
  initial guess for OH. Units of molec``^{-1}`` cm``^{3}`` s``^{-1}``.
* `k2eff` (default: ``8 × 10^{-12}``) - effective rate constant for the NO + RO2 → NO2 + RO reaction 
  *only* in the initial guess for OH. Units of molec``^{-1}`` cm``^{3}`` s``^{-1}``.
* `k5eff` (default: ``5 × 10^{-12}``) - effective rate constant for the RO2 and HO2 self-reaction 
  *only* in the initial guess for OH. Units of molec``^{-1}`` cm``^{3}`` s``^{-1}``.
"""
function SteadyStateOptions(;T::Real=298, M::Real=2e19, rh::Real=0.01, 
                            k_RO2NO::Real=8e-12, k_RO2HO2::Real=8e-12, k_RO2RO2::Real=6.8e-14,
                            k4::Real=1.1e-11, k2eff::Real=8e-12, k5eff::Real=5e-12)

    return SteadyStateOptions(T, M, rh, k_RO2NO, k_RO2HO2, k_RO2RO2, k4, k2eff, k5eff);
end

"""
Struct to hold the final state of the NOx-HOx steady state model. All concentrations
are in molec. cm``^{-3}`` s``^{-1}``. The following fields have the steady state
solution:

* `oh` - OH steady-state concentration
* `ho2` - HO2 steady-state concentration
* `ro2` - RO2 steady-state concentration

Additional fields store input or internal variables to aid in plotting,
calculating derived quantities, or debugging:

* `no` - NO concentration (from input)
* `no2` - NO2 concentration (from input)
* `vocr` - VOC reactivity in s``^{-1}``
* `phox` - HOx production rate in molec. cm``^{-3}`` s``^{-1}``
* `alpha` - NO + RO2 branching ratio
* `rates` - a dictionary of all rate constants (both specified and calculated) used in the solution
* `options` - the input `SteadyStateOptions` structure.
* `solver_results` - the `NLSolve.SolverResults` structure returned by the solver that found the
  optimum solution to the NOx-HOx steady state equations.
"""
struct SteadyStateResult
    no::Real
    no2::Real
    oh::Real
    ho2::Real
    ro2::Real
    vocr::Real
    phox::Real
    alpha::Real
    rates::Dict
    options::SteadyStateOptions
    solver_results::NLsolve.SolverResults
end


function nox_to_no_and_no2(nox, no2_no)
    no = 1/(no2_no + 1) * nox;
    no2 = no2_no/(no2_no + 1) * nox;
    return no, no2;
end


"""
    nonlin_nox_analytic_model(no::Real, no2::Real; kwargs...)::Real

Compute OH concentrations given NO and NO2 concentrations along
with additional kinetics parameters defined by the keyword argument.
All concentrations are in molec. cm``^{-3}``.

This uses the steady state model described in Murphy et al., ACP, 2006: 
"The weekend effect within and downwind of Sacramento: Part 2. 
Observational evidence for chemical and dynamical contributions"
(doi: 10.5194/acpd-6-11971-2006). See Eqs. 3, 4, and 5.

The keyword arguments, with defaults in  are:

* phox (6.25e6) - HOx production rate in molec. cm``^{-3}`` s``^{-1}``.
* vocr (5.8) - Total VOC OH reactivity in s``^{-1}``
* alpha (0.04) - RO2 + NO branching ratio, unitless.
* k4 (1.1e-11) - Rate constant for OH + NO2 --> HNO3 in cm``^3`` molec``^{-1}`` s``^{-1}``
* k2eff (8e-12) - Effective reaction rate of NO with RO2 in cm``^3`` molec``^{-1}`` s``^{-1}``
* k5eff (5e-12) - Effective reaction rate of RO2 and HO2 self reaction in cm``^3`` molec``^{-1}`` s``^{-1}``
"""
function nonlin_nox_analytic_model(no::Real, no2::Real; phox::Real=6.25e6, vocr::Real=5.8, alpha::Real=0.04, options::SteadyStateOptions=SteadyStateOptions())::Real
    k4 = options.k4;
    k2eff = options.k2eff;
    k5eff = options.k5eff;

    # The equation from Murphy et al. is basically the quadratic equation
    a = 6 * k5eff * (vocr / (k2eff * no))^2;
    b = k4 * no2 + alpha * vocr;
    c = -phox;

    return (-b + sqrt( b^2 - 4*a*c)) / (2 * a);
end 

"""
    nonlin_nox_analytic_model(nox::Real; no2_no::Real=4, kwargs...)::Real

Compute OH concentration from NOx concentration and an NO2:NO ratio, which defaults to 4:1.
"""
function nonlin_nox_analytic_model(nox::Real; no2_no::Real=4, kwargs...)::Real
    no, no2 = nox_to_no_and_no2(nox, no2_no)
    return nonlin_nox_analytic_model(no, no2; kwargs...);
end


"""
    hox_ss_solver(no::Real, no2::Real, phox::Real, vocr::Real, alpha::Real; no2_no::Real=4, options::SteadyStateOptions=SteadyStateOptions())::SteadyStateResult

Solve the NOx-HOx steady state system for a given set of NOx, P(HOx), VOC reactivity, and
RO2 branching ratio. The positional inputs are:

* `no` - NO concentration in molec. cm``^{-3}``
* `no2` - NO2 concentration in molec. cm``^{-3}``
* `phox` - HOx production rate in molec. cm``^{-3}`` s``^{-1}``
* `alpha` - NO + RO2 branching ratio

`options` must be a `SteadyStateOptions` instance, it can be used to change more detailed options
in the model. 

This returns a `SteadyStateResult` instance which contains the final state of the model.
"""
function hox_ss_solver(no::Real, no2::Real, phox::Real, vocr::Real, alpha::Real; no2_no::Real=4, options::SteadyStateOptions=SteadyStateOptions())::SteadyStateResult
    T = options.T;
    M = options.M;
    h2o = options.rh * M;

    k_RO2NO = options.k_RO2NO;
    k_RO2HO2 = options.k_RO2HO2;
    k_RO2RO2 = options.k_RO2RO2;
    k_HO2NO = Rates.kNOHO2(T, M);
    k_HO2HO2 = Rates.kHO2self(T, M, h2o);
    k_OHNO2 = Rates.KOHNO2a(T, M);

    rates = Dict("RO2+NO"=>k_RO2NO, "RO2+HO2"=>k_RO2HO2, "RO2+RO2"=>k_RO2RO2, "HO2+NO"=>k_HO2NO, "HO2+HO2"=>k_HO2HO2, "NO2+OH"=>k_OHNO2);

    x_initial = zeros(3);
    x_initial[1] = nonlin_nox_analytic_model(no, no2; phox=phox, vocr=vocr, alpha=alpha, options=options);
    x_initial[2] = x_initial[1] * vocr / (k_RO2NO * no);
    x_initial[3] = x_initial[2];

    # Convert from molec/cm3 to ppt to try to improve numeric stability
    # Was useful in the Matlab version, seems less so in Julia
    mcc_per_ppt = M / 1e12;

    no /= mcc_per_ppt;
    no2 /= mcc_per_ppt;
    phox /= mcc_per_ppt;
    k_RO2NO *= mcc_per_ppt;
    k_HO2NO *= mcc_per_ppt;
    k_HO2HO2 *= mcc_per_ppt;
    k_RO2HO2 *= mcc_per_ppt;
    k_RO2RO2 *= mcc_per_ppt;
    k_OHNO2 *= mcc_per_ppt;
    x_initial /= mcc_per_ppt

    # Solve the system of equations:
    #   [HO2] = (k_RO2NO * [RO2] * [NO] * (1-alpha) ) / (k_HO2NO * [NO] + 2 * k_HO2HO2 * [HO2] + k_RO2HO2 * [RO2])
    #   [RO2] = (VOCR * [OH]) / (k_RO2NO * [NO] + k_RO2HO2 * [HO2] + 2 * k_RO2RO2 * [RO2])
    #   PHOx  = k_OHNO2 * [OH] * [NO2] + alpha * k_RO2NO * [RO2] * [NO] + 2 * k_RO2HO2 * [RO2] * [HO2] + 2 * k_RO2RO2 * [RO2]^2 + 2 * k_HO2HO2 * [HO2]^2
    #
    # These solve the same system of reactions in Murphy 2006, but do not assume that [RO2] == [HO2]
    # For the solver function, x = ([HO2], [RO2], [OH])
    function f!(F, x)
        ho2 = x[1];
        ro2 = x[2];
        oh  = x[3];

        F[1] = (k_RO2NO * ro2 * no * (1-alpha) ) / (k_HO2NO * no + 2 * k_HO2HO2 * ho2 + k_RO2HO2 * ro2) - ho2;
        F[2] = (vocr * oh) / (k_RO2NO * no + k_RO2HO2 * ho2 + 2 * k_RO2RO2 * ro2) - ro2;
        F[3] = k_OHNO2 * oh * no2 + alpha * k_RO2NO * ro2 * no + 2 * k_RO2HO2 * ro2 * ho2 + 2 * k_RO2RO2 * ro2.^2 + 2 * k_HO2HO2 * ho2.^2 - phox;
    end

    result = nlsolve(f!, x_initial, autodiff=:forward, iterations=100000, method=:newton);
    ho2, ro2, oh = result.zero .* mcc_per_ppt;
    #ho2, ro2, oh = result.zero;

    
    
    SteadyStateResult(no*mcc_per_ppt, no2*mcc_per_ppt, oh, ho2, ro2, vocr, phox*mcc_per_ppt, alpha, rates, options, result)
    #SteadyStateResult(no, no2, oh, ho2, ro2, vocr, phox, alpha, rates, options, result)
end


end # module
