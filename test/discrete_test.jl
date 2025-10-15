using TestItems
using ProbabilityDistributions
using Distributions

@testitem "Generic – BetaNegBinomial" tags=[:Generic, :Discrete, :BetaNegBinomial] setup=[M] begin
    M.check(BetaNegBinomial(2, 1.5, 0.8))
    M.check(BetaNegBinomial(5, 2.0, 3.0))
    M.check(BetaNegBinomial(10, 1.5, 0.5))
    M.check(BetaNegBinomial(22, 0.5, 1.1))
end

@testitem "Generic – Borel" tags=[:Generic, :Discrete, :Borel] setup=[M] begin
    M.check(Borel(0.1))
    M.check(Borel(0.5))
    M.check(Borel(0.9))
end

@testitem "Generic – Conway" tags=[:Generic, :Discrete, :Conway] setup=[M] begin
    M.check(Conway(0.1, 0.0))
    M.check(Conway(1.0, 0.5))
    M.check(Conway(2.5, 1.0))
end

@testitem "Generic – Delaporte" tags=[:Generic, :Discrete, :Delaporte] setup=[M] begin
    M.check(Delaporte(0.1, 0.5, 0.5))
    M.check(Delaporte(2.5, 2.0, 1.0))
    M.check(Delaporte(5.0, 3.0, 2.0))
end

@testitem "Generic – FlorySchulz" tags=[:Generic, :Discrete, :FlorySchulz] setup=[M] begin
    M.check(FlorySchulz(0.1))
    M.check(FlorySchulz(0.7))
    M.check(FlorySchulz(0.9))
end

@testitem "Generic – GaussKuzmin" tags=[:Generic, :Discrete, :GaussKuzmin] setup=[M] begin
    M.check(GaussKuzmin())
end

@testitem "Generic – Logarithmic" tags=[:Generic, :Discrete, :Logarithmic] setup=[M] begin
    M.check(Logarithmic(0.1))
    M.check(Logarithmic(0.5))
    M.check(Logarithmic(0.9))
end

@testitem "Generic – Rademacher" tags=[:Generic, :Discrete, :Rademacher] setup=[M] begin
    M.check(Rademacher())
end

@testitem "Generic – Yule" tags=[:Generic, :Discrete, :Yule] setup=[M] begin
    M.check(Yule(0.8))   # caso con media infinita
    M.check(Yule(2.0))   # caso con media finita
end

@testitem "Generic – Zeta" tags=[:Generic, :Discrete, :Zeta] setup=[M] begin
    M.check(Zeta(1.5))   # media infinita
    M.check(Zeta(2.5))   # media finita
    M.check(Zeta(3.5))
end

@testitem "Generic – ZIB" tags=[:Generic, :Discrete, :ZIB] setup=[M] begin
    M.check(ZIB(10, 0.3, 0.7))
    M.check(ZIB(20, 0.2, 0.9))
end

@testitem "Generic – ZINB" tags=[:Generic, :Discrete, :ZINB] setup=[M] begin
    M.check(ZINB(10, 0.3, 0.7))
    M.check(ZINB(20, 0.4, 0.5))
end

@testitem "Generic – ZIP" tags=[:Generic, :Discrete, :ZIP] setup=[M] begin
    M.check(ZIP(2.0, 0.3))
    M.check(ZIP(5.0, 0.7))
end

@testitem "Generic – Zipf" tags=[:Generic, :Discrete, :Zipf] setup=[M] begin
    M.check(Zipf(10, 1.0))
    M.check(Zipf(50, 2.0))
end

