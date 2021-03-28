module SteadyStateSolvers

using NLsolve;
import ..Rates;

export nonlin_nox_analytic_model,
       hox_ss_solver;


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

function SteadyStateOptions(;T::Real=298, M::Real=2e19, rh::Real=0.01, 
                            k_RO2NO::Real=8e-12, k_RO2HO2::Real=8e-12, k_RO2RO2::Real=6.8e-14,
                            k4::Real=1.1e-11, k2eff::Real=8e-12, k5eff::Real=5e-12)

    return SteadyStateOptions(T, M, rh, k_RO2NO, k_RO2HO2, k_RO2RO2, k4, k2eff, k5eff);
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
All concentrations are in molec./cm^3

This uses the steady state model described in Murphy et al., ACP, 2006: 
"The weekend effect within and downwind of Sacramento: Part 2. 
Observational evidence for chemical and dynamical contributions"
(doi: 10.5194/acpd-6-11971-2006). See Eqs. 3, 4, and 5.

The keyword arguments, with defaults in  are:

* phox (6.25e6) - HOx production rate in molec. cm^-3 s^-1.
* vocr (5.8) - Total VOC OH reactivity in s^-1
* alpha (0.04) - RO2 + NO branching ratio, unitless.
* k4 (1.1e-11) - Rate constant for OH + NO2 --> HNO3 in cm^3 molec.^-1 s^-1
* k2eff (8e-12) - Effective reaction rate of NO with RO2 in cm^3 molec.^-1 s^-1
* k5eff (5e-12) - Effective reaction rate of RO2 and HO2 self reaction in cm^3 molec.^-1 s^-1
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


function hox_ss_solver(no::Real, no2::Real, phox::Real, vocr::Real, alpha::Real; no2_no::Real=4, options::SteadyStateOptions=SteadyStateOptions())::Array{<:Real,1}
    T = options.T;
    M = options.M;
    h2o = options.rh * M;

    k_RO2NO = options.k_RO2NO;
    k_RO2HO2 = options.k_RO2HO2;
    k_RO2RO2 = options.k_RO2RO2;
    k_HO2NO = Rates.kNOHO2(T, M);
    k_HO2HO2 = Rates.kHO2self(T, M, h2o);
    k_OHNO2 = Rates.KOHNO2a(T, M);

    # Convert from molec/cm3 to ppt to try to improve numeric stability
    mcc_per_ppt = M / 1e12;

    no /= mcc_per_ppt;
    no2 /= mcc_per_ppt;
    phox /= mcc_per_ppt;
    k_RO2HO2 *= mcc_per_ppt;
    k_HO2NO *= mcc_per_ppt;
    k_HO2HO2 *= mcc_per_ppt;
    k_RO2HO2 *= mcc_per_ppt;
    k_RO2RO2 *= mcc_per_ppt;
    k_OHNO2 *= mcc_per_ppt;

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

    x_initial = zeros(3);
    x_initial[1] = nonlin_nox_analytic_model(no, no2; phox=phox, vocr=vocr, alpha=alpha, options=options);
    x_initial[2] = x_initial[1] * vocr / (k_RO2NO * no);
    x_initial[3] = x_initial[2];

    result = nlsolve(f!, x_initial, autodiff=:forward);
    concentrations = result.zero .* mcc_per_ppt;
    return concentrations;
end


end # module