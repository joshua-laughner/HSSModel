"""
Kinetic rate constants needed for the steady state models.
For all rate functions, `T` is the temperature in K and 
`M` the number density of air in molec. cm^-3. For rate functions
that include `h2o` is the number density of water, also in
molec. cm^-3.
"""
module Rates

"Reaction rate for HO2 + NO -> NO2 + OH"
kNOHO2(T::Real,M::Real) = 3.5e-12 * exp(250.0/T)  # Based on JPL data evaluation #15

"Reaction rate for HO2 + HO2"
kHO2self(T::Real,M::Real,h2o::Real) = (3.5e-13*exp(430.0/T) + 1.7e-33*(M-h2o)*exp(1000.0/T))*(1 + 1.4e-21*h2o*exp(2200.0/T));  # Based on JPL data evaluation #15

"Reaction rate for OH + NO2 -> HNO3"
function KOHNO2a(T::Real,M::Real)
    # From Mollner et al. 2010 (Science, doi: 10.1126/science.1193030 )
    # Yes that is T raised to the 0 to remove the T dependence.
    k9o= 1.51e-30*(T/300.0)^0 * M;
    k9oo=2.58e-11*(T./300)^0;
    k9=(k9o/(1 + (k9o/k9oo)))*0.6^((1+(log10(k9o/k9oo))^2)^(-1.0));
    return k9; 
end

end
