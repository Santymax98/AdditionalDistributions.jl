@testitem "Generic – Alpha" tags=[:generic, :continuous, :alpha] setup=[M] begin
    M.check(Alpha(1.0, 1.0))
    M.check(Alpha(3.0, 4.0))
    M.check(Alpha(5.0, 1.0))
end

@testitem "Generic – Argus" tags=[:generic, :continuous, :argus] setup=[M] begin
    M.check(Argus(1.0, 1.0))
    M.check(Argus(2.0, 2.0))
    M.check(Argus(3.0, 4.0))
end

@testitem "Generic – Asymmetric Laplace" tags=[:generic, :continuous, :asymlaplace] setup=[M] begin
    M.check(AsymLaplace(0.0, 1.0, 1.0))
    M.check(AsymLaplace(2.0, 3.0, 0.5))
    M.check(AsymLaplace(-1.0, 2.0, 2.0))
end


@testitem "Generic – Benini" tags=[:generic, :continuous, :benini] setup=[M] begin
    M.check(Benini(1.0, 1.0, 1.0))
    M.check(Benini(3.0, 4.0, 5.0))
    M.check(Benini(5.0, 1.0, 2.0))
end

@testitem "Generic – Benktander_Type1" tags=[:generic, :continuous, :benktander_type1] setup=[M] begin
    M.check(Benktander_Type1(1.0, 0.5))
    M.check(Benktander_Type1(3.0, 3.0))
    M.check(Benktander_Type1(5.0, 7.5))
end

@testitem "Generic – Benktander_Type2" tags=[:generic, :continuous, :benktander_type2] setup=[M] begin
    M.check(Benktander_Type2(1.0, 0.5))
    M.check(Benktander_Type2(3.0, 0.8))
    M.check(Benktander_Type2(10.0, 0.2))
end

@testitem "Generic – Bhattacharjee" tags=[:generic, :continuous, :bhattacharjee] setup=[M] begin
    M.check(Bhattacharjee(1.0, 2.0, 0.5))
    M.check(Bhattacharjee(0.0, 1.0, 0.1))
    M.check(Bhattacharjee(-2.0, -1.0, 0.3))
end

@testitem "Generic – BirnbaumSaunders" tags=[:generic, :continuous, :birnbaumsaunders] setup=[M] begin
    M.check(BirnbaumSaunders(0.0, 1.0, 1.0))
    M.check(BirnbaumSaunders(1.0, 0.5, 2.0))
    M.check(BirnbaumSaunders(-2.0, 0.1, 0.1))
end

@testitem "Generic – Bradford" tags=[:generic, :continuous, :bradford] setup=[M] begin
    M.check(Bradford(1.0))
    M.check(Bradford(5.5))
    M.check(Bradford(12.1))
end

@testitem "Generic – Burr" tags=[:generic, :continuous, :burr] setup=[M] begin
    M.check(Burr(1.0, 1.0, 1.0))
    M.check(Burr(2.0, 3.0, 4.0))
    M.check(Burr(0.2, 0.3, 0.4))
end

@testitem "Generic – CrystalBall" tags=[:generic, :continuous, :crystalball] setup=[M] begin
    M.check(CrystalBall(1.5, 3.5, 0.0, 1.0))
    M.check(CrystalBall(2.0, 5.0, 1.0, 0.5))
end

@testitem "Generic – Dagum" tags=[:generic, :continuous, :dagum] setup=[M] begin
    M.check(Dagum(0.5, 1.5, 2.5))
    M.check(Dagum(2.0, 3.0, 4.0))
    M.check(Dagum(0.2, 0.3, 0.4))
end

@testitem "Generic – Gompertz" tags=[:generic, :continuous, :gompertz] setup=[M] begin
    M.check(Gompertz(0.5, 1.5))
    M.check(Gompertz(2.0, 3.0))
    M.check(Gompertz(0.1, 0.2))
end

@testitem "Generic – HalfCauchy" tags=[:generic, :continuous, :halfcauchy] setup=[M] begin
    M.check(HalfCauchy())
    M.check(HalfCauchy(0.5))
    M.check(HalfCauchy(3.0))
end

@testitem "Generic – HalfNormal" tags=[:generic, :continuous, :halfnormal] setup=[M] begin
    M.check(HalfNormal())
    M.check(HalfNormal(0.5))
    M.check(HalfNormal(3.0))
end

@testitem "Generic – HalfTDist" tags=[:generic, :continuous, :halftdist] setup=[M] begin
    M.check(HalfTDist())
    M.check(HalfTDist(0.5))
    M.check(HalfTDist(3.0))
end

@testitem "Generic – IrwinHall" tags=[:generic, :continuous, :irwinhall] setup=[M] begin
    M.check(IrwinHall(2))
    M.check(IrwinHall(4, -1.0, 1.0))
    M.check(IrwinHall(10, 2.0, 5.0))
end

@testitem "Generic – Lomax" tags=[:generic, :continuous, :lomax] setup=[M] begin
    M.check(Lomax(0.5, 1.5))
    M.check(Lomax(2.0, 3.0))
    M.check(Lomax(1.0, 0.2))
end

@testitem "Generic – Maxwell" tags=[:generic, :continuous, :maxwell] setup=[M] begin
    M.check(Maxwell(1.0))
    M.check(Maxwell(5.5))
    M.check(Maxwell(12.1))
end

@testitem "Generic – Nakagami" tags=[:generic, :continuous, :nakagami] setup=[M] begin
    M.check(Nakagami(0.6, 1.5))
    M.check(Nakagami(2.0, 3.0))
    M.check(Nakagami(3.1, 0.2))
end

@testitem "Generic – PERT" tags=[:generic, :continuous, :pert] setup=[M] begin
    M.check(PERT(0.0, 1.0, 2.0))
    M.check(PERT(-1.0, 0.0, 1.0))
    M.check(PERT(5.0, 6.0, 7.0))
end