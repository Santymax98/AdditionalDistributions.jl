@testitem "Generic – Alpha" tags=[:Generic, :Continuous, :Alpha] setup=[M] begin
    M.check(Alpha(1.0, 1.0))
    M.check(Alpha(3.0, 4.0))
    M.check(Alpha(5.0, 1.0))
end

@testitem "Generic – Argus" tags=[:Generic, :Continuous, :Argus] setup=[M] begin
    M.check(Argus(1.0, 1.0))
    M.check(Argus(2.0, 2.0))
    M.check(Argus(3.0, 4.0))
end

@testitem "Generic – Benini" tags=[:Generic, :Continuous, :Benini] setup=[M] begin
    M.check(Benini(1.0, 1.0, 1.0))
    M.check(Benini(3.0, 4.0, 5.0))
    M.check(Benini(5.0, 1.0, 2.0))
end

@testitem "Generic – Benktander_Type1" tags=[:Generic, :Continuous, :Benktander_Type1] setup=[M] begin
    M.check(Benktander_Type1(1.0, 0.5))
    M.check(Benktander_Type1(3.0, 3.0))
    M.check(Benktander_Type1(5.0, 7.5))
end

@testitem "Generic – Benktander_Type2" tags=[:Generic, :Continuous, :Benktander_Type2] setup=[M] begin
    M.check(Benktander_Type2(1.0, 0.5))
    M.check(Benktander_Type2(3.0, 0.8))
    M.check(Benktander_Type2(10.0, 0.2))
end

@testitem "Generic – Bhattacharjee" tags=[:Generic, :Continuous, :Bhattacharjee] setup=[M] begin
    M.check(Bhattacharjee(1.0, 2.0, 0.5))
    M.check(Bhattacharjee(0.0, 1.0, 0.1))
    M.check(Bhattacharjee(-2.0, -1.0, 0.3))
end

@testitem "Generic – BirnbaumSaunders" tags=[:Generic, :Continuous, :BirnbaumSaunders] setup=[M] begin
    M.check(BirnbaumSaunders(0.0, 1.0, 1.0))
    M.check(BirnbaumSaunders(1.0, 0.5, 2.0))
    M.check(BirnbaumSaunders(-2.0, 0.1, 0.1))
end

@testitem "Generic – Bradford" tags=[:Generic, :Continuous, :Bradford] setup=[M] begin
    M.check(Bradford(1.0))
    M.check(Bradford(5.5))
    M.check(Bradford(12.1))
end

@testitem "Generic – Burr" tags=[:Generic, :Continuous, :Burr] setup=[M] begin
    M.check(Burr(1.0, 1.0, 1.0))
    M.check(Burr(2.0, 3.0, 4.0))
    M.check(Burr(0.1, 0.2, 0.3))
end

@testitem "Generic – Dagum" tags=[:Generic, :Continuous, :Dagum] setup=[M] begin
    M.check(Dagum(0.5, 1.5, 2.5))
    M.check(Dagum(2.0, 3.0, 4.0))
    M.check(Dagum(0.2, 0.3, 0.4))
end

@testitem "Generic – Gompertz" tags=[:Generic, :Continuous, :Gompertz] setup=[M] begin
    M.check(Gompertz(0.5, 1.5))
    M.check(Gompertz(2.0, 3.0))
    M.check(Gompertz(0.1, 0.2))
end

@testitem "Generic – Lomax" tags=[:Generic, :Continuous, :Lomax] setup=[M] begin
    M.check(Lomax(0.5, 1.5))
    M.check(Lomax(2.0, 3.0))
    M.check(Lomax(1.0, 0.2))
end

@testitem "Generic – Maxwell" tags=[:Generic, :Continuous, :Maxwell] setup=[M] begin
    M.check(Maxwell(1.0))
    M.check(Maxwell(5.5))
    M.check(Maxwell(12.1))
end

@testitem "Generic – Nakagami" tags=[:Generic, :Continuous, :Nakagami] setup=[M] begin
    M.check(Nakagami(0.6, 1.5))
    M.check(Nakagami(2.0, 3.0))
    M.check(Nakagami(3.1, 0.2))
end

@testitem "Generic – PERT" tags=[:Generic, :Continuous, :PERT] setup=[M] begin
    M.check(PERT(0.0, 1.0, 2.0))
    M.check(PERT(-1.0, 0.0, 1.0))
    M.check(PERT(5.0, 6.0, 7.0))
end