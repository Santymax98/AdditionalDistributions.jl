using TestItems
using AdditionalDistributions
using Distributions

@testitem "Generic – BetaNegBinomial" tags=[:generic, :discrete, :betanegbinomial] setup=[M] begin
    M.check(BetaNegBinomial(2, 1.5, 0.8))
    M.check(BetaNegBinomial(5, 2.0, 3.0))
    M.check(BetaNegBinomial(10, 1.5, 0.5))
    #M.check(BetaNegBinomial(22, 0.5, 1.1))
end

@testitem "Generic – Borel" tags=[:generic, :discrete, :borel] setup=[M] begin
    M.check(Borel(0.1))
    M.check(Borel(0.5))
    M.check(Borel(0.9))
end

@testitem "Generic – Conway" tags=[:generic, :discrete, :conway] setup=[M] begin
    M.check(Conway(0.1, 0.0))
    M.check(Conway(1.0, 0.5))
    M.check(Conway(2.5, 1.0))
end

@testitem "Generic – Delaporte" tags=[:generic, :discrete, :delaporte] setup=[M] begin
    M.check(Delaporte(0.1, 0.5, 0.5))
    M.check(Delaporte(2.5, 2.0, 1.0))
    M.check(Delaporte(5.0, 3.0, 2.0))
end

@testitem "Generic – FlorySchulz" tags=[:generic, :discrete, :floryschulz] setup=[M] begin
    M.check(FlorySchulz(0.1))
    M.check(FlorySchulz(0.7))
    M.check(FlorySchulz(0.9))
end

@testitem "Generic – GaussKuzmin" tags=[:generic, :discrete, :gausskuzmin] setup=[M] begin
    M.check(GaussKuzmin())
end

@testitem "Generic – Logarithmic" tags=[:generic, :discrete, :logarithmic] setup=[M] begin
    M.check(Logarithmic(0.1))
    M.check(Logarithmic(0.5))
    M.check(Logarithmic(0.9))
end

@testitem "Generic – PoissonInvGaussian" tags=[:generic, :discrete, :poissoninvgaussian] setup=[M] begin
    M.check(PoissonInvGaussian(2.0, 50.0))  # low overdispersion, close to Poisson(2)
    M.check(PoissonInvGaussian(10.0, 5.0))  # moderate overdispersion
    M.check(PoissonInvGaussian(25.0, 1.5))  # strong overdispersion and heavier tail
end

@testitem "Generic – Rademacher" tags=[:generic, :discrete, :rademacher] setup=[M] begin
    M.check(Rademacher())
end

@testitem "Generic – Weibull_Type1" tags=[:generic, :discrete, :weibull_type1] setup=[M] begin
    M.check(Weibull_Type1(0.5, 1.0)) # Geometric distribution
    M.check(Weibull_Type1(0.8, 1.5)) # PMF decreasing
    M.check(Weibull_Type1(0.9, 3.0)) # non-trivial case
end


@testitem "Generic – Yule" tags=[:generic, :discrete, :yule] setup=[M] begin
    M.check(Yule(0.8))   # caso con media infinita
    M.check(Yule(2.0))   # caso con media finita
end

@testitem "Generic – Zeta" tags=[:generic, :discrete, :zeta] setup=[M] begin
    M.check(Zeta(1.5))   # media infinita
    M.check(Zeta(2.5))   # media finita
    M.check(Zeta(3.5))
end

@testitem "Generic – ZIB" tags=[:generic, :discrete, :zib] setup=[M] begin
    M.check(ZIB(10, 0.3, 0.7))
    M.check(ZIB(20, 0.2, 0.9))
end

@testitem "Generic – ZINB" tags=[:generic, :discrete, :zinb] setup=[M] begin
    M.check(ZINB(10, 0.3, 0.7))
    M.check(ZINB(20, 0.4, 0.5))
end

@testitem "Generic – ZIP" tags=[:generic, :discrete, :zip] setup=[M] begin
    M.check(ZIP(2.0, 0.3))
    M.check(ZIP(5.0, 0.7))
end

@testitem "Generic – Zipf" tags=[:generic, :discrete, :zipf] setup=[M] begin
    M.check(Zipf(10, 1.0))
    M.check(Zipf(50, 2.0))
end

