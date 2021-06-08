module DerivedQuantities

import ..Rates;
import ..SteadyStateSolvers;

export nox_lifetime

"""
    nox_lifetime(result::SteadyStateSolvers.SteadyStateResult)::Dict

Compute NOx lifetimes from a steady state model result. The return
value is a dictionary with keys "total", "hno3", and "ans" which correspond
to the overall lifetime, lifetime w.r.t. loss to HNO3, and lifetime w.r.t.
loss to alkyl nitrates, respectively. Lifetimes are in hours.
"""
function nox_lifetime(result::SteadyStateSolvers.SteadyStateResult)::Dict
    # Get quantities out of the result 
    nox = result.no + result.no2;
    T = result.options.T;
    M = result.options.M;

    k2eff = result.options.k2eff;
    k_OHNO2 = result.options.k4;
    k_HO2NO = Rates.kNOHO2(T, M);
    k_RO2HO2 = result.options.k_RO2HO2;
    k_RO2RO2 = result.options.k_RO2RO2;

    tau_hno3 = nox / (k_OHNO2 * result.oh * result.no2 * 3600); # 3600 converts seconds to hours
    # From Murphy 2006, ACPD, we assume some effective k_RO2_NO (k2eff) that is
    # the weighted average of k's for various RO2+NO reactions. [RO2] is
    # calculated by assuming:
    #  1) RO2 + HO2 is in steady state
    #  2) All RO2's go on to produce HO2, thus [RO2] = [HO2]
    #
    # This means P(RO2) = VOCR*[OH] == L(HO2) = k_HO2+NO * [HO2] * [NO]
    # => [HO2] = [RO2] = VOCR * [OH] / (k_HO2+NO * [NO])
    #
    # This gives us an expression for [RO2], however, we need k_RO2+NO. In
    # Murphy, k2eff is the effective rate of NO oxidation by peroxy radicals,
    # so
    #
    #   k2_eff = ( [HO2]*k_HO2+NO + \sum_i (k_RO2_i * [RO2]_i) )/( [HO2 + \sum_i [RO2]_i )
    #
    # But [HO2] = \sum_i [RO2]_i (let's call this C) so that this becomes:
    #
    #   k2_eff = ( C*k_HO2+NO + C*k_RO2+NO ) / 2C 
    # => 2*k2_eff - k_HO2+NO = k_RO2+NO
    k_RO2NO = 2*k2eff - k_HO2NO;
    tau_ans = nox ./ (result.alpha .* k_RO2NO .* result.ro2 .* result.no .* 3600);

    tau = (1/tau_hno3 + 1/tau_ans)^(-1)
    return Dict("total" => tau, "hno3" => tau_hno3, "ans" => tau_ans)
end


end
