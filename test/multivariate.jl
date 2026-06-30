# =============================================================
# Validation test of MvTStudent against mvtnorm::pmvt(R)
# ==============================================================
#
# This test suite uses the same scenarios by Genz (2002)
# used in the `MvNormalCDF.jl` package to verify the accuracy
# of multivariate normal integrators.
#
# In this case, the same `GenzTestData.td` configuration matrix has been taken
# (identical to the one used in `MvNormalCDF.jl`) and adapted to be
# compatible with `MvGaussian` and `MvTStudent`, which are equivalent in structure
# but use their own QMC engine for CDF evaluation.
#
# The numerical reference comes from `mvtnorm::pmvt()` in R, using
# `sigma = Σ` (full covariance matrix) and `delta = 0` (central t).
# The comparison is performed on probabilities of rectangles defined
# by the Genz limits (a,b), for different degrees of freedom ν.
#
# The tables `P_REF` and `E_REF` were generated directly in R with:
# pmvt(lower=a, upper=b, df=ν, sigma=Σ, delta=rep(0, d),
# algorithm=GenzBretz(maxpts=m, abseps=1e-8, releps=1e-8))
#
# ----------------------------------------------------------
# Author: Santiago Jiménez (2025)
# Based on: Genz & Bretz (2002), mvtnorm R
# ===========================================================

@testmodule GenzTestData begin

td = Array{Any}(undef,(14,9))
#  1-cov mtx 2-a 3-b 4-m 5-p 6-ptol 7-e 8-etol

# from MATLAB documenation 4 dim
td[1,1] = [4 3 2 1;3 5 -1 1;2 -1 4 2;1 1 2 5]  # Σ cov Matrix
td[1,2] = [-Inf; -Inf; -Inf; -Inf]             # a lower integration limit
td[1,3] = [1; 2; 3; 4 ]                        # b upper integration limit
td[1,4] = 5000
td[1,5] = 0.605653                              # expected p value
td[1,6] = 0.001557374                           # ± p tolerance
td[1,7] = 0.001394971                           # expected e (error) value
td[1,8] = 0.0009277058                          # ± e tolerance

# 3 dim
td[2,1] =[1  3/5  1/3; 3/5 1  11/15; 1/3 11/15 1]
td[2,2] = [-Inf;-Inf;-Inf]
td[2,3] = [1;4;2]
td[2,4] = 3000
td[2,5] = 0.827985
td[2,6] = 3.529068e-5
td[2,7] = 2.8608229e-5
td[2,8] = 3.72604931e-5

# 3 dim
td[3,1] = [1 0.25 0.2; 0.25 1 0.333333333; 0.2 0.333333333 1]
td[3,2] = [-1;-4;-2]
td[3,3] = [1;4;2]
td[3,4] = 4000
td[3,5] = 0.65368
td[3,6] = 0.000002089699
td[3,7] = 0.000001799225
td[3,8] = 0.000001161755

# Genz book eq. 1.5 p. 4-5 & p. 63
# Genz gives wrong answer on p. 4-5
td[4,1] = [1/3 3/5 1/3; 3/5 1.0 11/15; 1/3 11/15 1.0]
td[4,2] = [-Inf; -Inf; -Inf]
td[4,3] = [1; 4; 2]
td[4,4] = 4000
td[4,5] = 0.943174
td[4,6] = 0.00006724871
td[4,7] = 0.0000509768
td[4,8] = 0.00003337435

# Genz p. 63 uses different r matrix
td[5,1] = [1 0 0; 3/5 1 0; 1/3 11/15 1]
td[5,2] = [-Inf; -Inf; -Inf]
td[5,3] = [1; 4; 2]
td[5,4] = 4000
td[5,5] = 0.827985
td[5,6] = 0.00001024322
td[5,7] = 0.000008933135
td[5,8] = 0.000005194712

# singular example
# problem reduces to univariate problem with
# p = cdf.(Normal(),1) = 0.84134476068543
td[6,1] = [1 1 1; 1 1 1; 1 1 1]
td[6,2] = [-Inf, -Inf, -Inf]
td[6,3] = [1, 1, 1]
td[6,4] = 3000
td[6,5] = 0.841345
td[6,6] = 4.440892E-16
td[6,7] = 0.0
td[6,8] = 2.220446E-16

# 5 dim example
td[7,1] = [1 1 1 1 1;
         1 2 2 2 2;
         1 2 3 3 3;
         1 2 3 4 4;
         1 2 3 4 5]
td[7,2] = [-1,-2,-3,-4,-5]
td[7,3] = [2,3,4,5,6]
td[7,4] = 6000
td[7,5] = 0.761243
td[7,6] = 0.0004715559
td[7,7] = 0.0003493157
td[7,8] = 0.0002319606

# Genz used wrong integration limits when computing
# above. see p. 63
td[8,1] = td[7,1]
td[8,2] = sort(td[7,2])
td[8,3] = 1 .- td[8,2]
td[8,4] = 6000
td[8,5] = 0.474128
td[8,6] = 0.000037335
td[8,7] = 0.00003169846
td[8,8] = 0.0000178166

# positive orthant probability of above
td[9,1] = td[7,1]
td[9,2] = [0,0,0,0,0]
td[9,3] = td[8,3]
td[9,4] = 6000
td[9,5] = 0.113537
td[9,6] = 0.0001740044
td[9,7] = 0.0001584365
td[9,8] = 0.0001053655

# test 7 Cov matrix, but now -Inf lower limit
td[10,1] = td[7,1]
td[10,2] = [-Inf,-Inf,-Inf,-Inf,-Inf]
td[10,3] = td[8,3]
td[10,4] = 6000
td[10,5] = 0.810315
td[10,6] = 0.00003136998
td[10,7] = 0.0000287006
td[10,8] = 0.00001567681

# eight dimensional test
td[11,1] = [1 1 1 1 1 1 1 1;
            1 2 2 2 2 2 2 2;
            1 2 3 3 3 3 3 3;
            1 2 3 4 4 4 4 4;
            1 2 3 4 5 5 5 5;
            1 2 3 4 5 6 6 6;
            1 2 3 4 5 6 7 7;
            1 2 3 4 5 6 7 8]
td[11,2] =  -1*[1,2,3,4,5,6,7,8]
td[11,3] = [2,3,4,5,6,7,8,9]
td[11,4] = 9000
td[11,5] = 0.759474
td[11,6] = 0.0006248131
td[11,7] = 0.0005161476
td[11,8] = 0.0003407714

# orthant probability of above
td[12,1] = td[11,1]
td[12,2] = [0,0,0,0,0,0,0,0]
td[12,3] = [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf]
td[12,4] = 9000
td[12,5] = 0.19638
td[12,6] = 0.0002654318
td[12,7] = 0.0002384606
td[12,8] = 0.0001693536

# with dim with lower limit -Inf (exact result)
td[13,1] = td[11,1]
td[13,2] = -Inf*[1,1,1,1,1,1,1,1]
td[13,3] = td[12,3]
td[13,4] = 9000
td[13,5] = 1.0
td[13,6] = 2.220446E-16
td[13,7] = 0.0
td[13,8] = 2.220446E-16

# 25 dimensions
td[14,1] = [59.227 2.601 3.38 8.303 -0.334 11.029 10.908 0.739 4.703 7.075 8.049 1.403 9.838 5.46 11.949 2.272 7.234 15.215 -9.091 12.265 3.01 -3.199 10.608 8.464 -8.685;
 2.601 77.213 0.882 2.99 -2.536 -4.55 -3.874 -3.607 6.023 3.129 15.7 -7.271 -10.655 8.456 15.387 -5.764 6.617 -6.331 -2.244 -0.925 -7.516 13.836 9.243 -0.84 -3.781;
 3.38 0.882 79.72 -1.465 3.179 -1.799 9.842 9.165 -1.54 -8.03 0.778 9.053 -2.598 -8.844 15.857 13.613 -1.878 8.18 11.806 -0.242 4.711 2.258 9.554 1.184 -8.047;
 8.303 2.99 -1.465 65.041 3.803 2.13 -0.936 -5.996 2.719 -4.648 4.611 4.486 13.38 -0.376 0.179 14.654 -7.089 -1.194 9.357 5.12 4.943 -0.475 4.764 -8.56 2.337;
 -0.334 -2.536 3.179 3.803 64.128 4.85 -14.767 -10.044 12.437 5.065 8.191 0.391 0.068 9.754 -0.062 0.429 9.265 5.502 4.227 0.559 0.811 3.169 4.558 -1.878 -3.885;
 11.029 -4.55 -1.799 2.13 4.85 64.912 12.989 2.675 1.227 8.205 3.3 3.545 6.225 11.936 -2.956 6.188 -2.206 1.184 10.546 2.492 12.035 -15.789 4.296 1.086 7.225;
 10.908 -3.874 9.842 -0.936 -14.767 12.989 78.424 2.03 5.595 4.921 1.573 -7.355 -9.425 -4.024 4.912 -12.05 -2.034 -2.435 -10.355 0.985 8.23 2.806 6.254 -1.494 5.531;
 0.739 -3.607 9.165 -5.996 -10.044 2.675 2.03 93.184 8.84 3.843 -14.968 16.386 -0.223 4.398 -4.786 1.731 4.025 0.479 3.12 -15.591 12.614 -8.279 -3.582 6.597 -1.915;
 4.703 6.023 -1.54 2.719 12.437 1.227 5.595 8.84 73.981 4.019 -6.404 5.869 -4.305 5.936 2.016 5.987 10.051 -0.705 13.229 -1.715 7.102 12.89 10.967 5.262 15.954;
 7.075 3.129 -8.03 -4.648 5.065 8.205 4.921 3.843 4.019 68.048 7.346 7.412 9.956 6.743 2.547 0.177 0.844 2.147 -4.072 11.832 -3.55 -0.096 -1.96 -1.381 -3.249;
 8.049 15.7 0.778 4.611 8.191 3.3 1.573 -14.968 -6.404 7.346 66.708 -6.22 -8.952 0.647 3.039 -12.078 7.618 10.398 -5.055 1.577 -12.77 11.477 8.272 2.071 7.728;
 1.403 -7.271 9.053 4.486 0.391 3.545 -7.355 16.386 5.869 7.412 -6.22 55.025 0.273 -12.049 -2.71 11.531 12.604 1.294 2.791 -5.698 -2.231 15.025 12.229 -5.876 -3.374;
 9.838 -10.655 -2.598 13.38 0.068 6.225 -9.425 -0.223 -4.305 9.956 -8.952 0.273 82.184 3.466 -9.297 -4.347 12.586 4.372 13.705 -6.795 -5.818 -6.78 5.11 10.099 5.05;
 5.46 8.456 -8.844 -0.376 9.754 11.936 -4.024 4.398 5.936 6.743 0.647 -12.049 3.466 69.352 10.855 7.282 6.615 4.58 0.306 6.482 14.589 5.081 -9.141 -4.657 10.763;
 11.949 15.387 15.857 0.179 -0.062 -2.956 4.912 -4.786 2.016 2.547 3.039 -2.71 -9.297 10.855 87.988 2.504 4.226 10.461 9.703 -3.112 -13.348 0.944 -2.824 -4.498 10.551;
 2.272 -5.764 13.613 14.654 0.429 6.188 -12.05 1.731 5.987 0.177 -12.078 11.531 -4.347 7.282 2.504 63.611 3.024 2.35 -5.191 -6.101 -6.324 -0.483 9.899 5.768 2.382;
 7.234 6.617 -1.878 -7.089 9.265 -2.206 -2.034 4.025 10.051 0.844 7.618 12.604 12.586 6.615 4.226 3.024 78.3 1.54 -6.868 2.613 6.006 5.49 9.06 -4.229 -4.395;
 15.215 -6.331 8.18 -1.194 5.502 1.184 -2.435 0.479 -0.705 2.147 10.398 1.294 4.372 4.58 10.461 2.35 1.54 61.499 11.083 15.428 1.771 2.517 5.181 13.476 6.829;
 -9.091 -2.244 11.806 9.357 4.227 10.546 -10.355 3.12 13.229 -4.072 -5.055 2.791 13.705 0.306 9.703 -5.191 -6.868 11.083 76.538 4.287 6.564 -2.49 15.558 11.202 -16.964;
 12.265 -0.925 -0.242 5.12 0.559 2.492 0.985 -15.591 -1.715 11.832 1.577 -5.698 -6.795 6.482 -3.112 -6.101 2.613 15.428 4.287 59.057 -5.619 4.374 8.934 -6.203 14.687;
 3.01 -7.516 4.711 4.943 0.811 12.035 8.23 12.614 7.102 -3.55 -12.77 -2.231 -5.818 14.589 -13.348 -6.324 6.006 1.771 6.564 -5.619 69.833 -1.233 -7.922 -7.027 -2.315;
 -3.199 13.836 2.258 -0.475 3.169 -15.789 2.806 -8.279 12.89 -0.096 11.477 15.025 -6.78 5.081 0.944 -0.483 5.49 2.517 -2.49 4.374 -1.233 65.458 5.786 4.613 2.223;
 10.608 9.243 9.554 4.764 4.558 4.296 6.254 -3.582 10.967 -1.96 8.272 12.229 5.11 -9.141 -2.824 9.899 9.06 5.181 15.558 8.934 -7.922 5.786 69.574 3.527 -2.802;
 8.464 -0.84 1.184 -8.56 -1.878 1.086 -1.494 6.597 5.262 -1.381 2.071 -5.876 10.099 -4.657 -4.498 5.768 -4.229 13.476 11.202 -6.203 -7.027 4.613 3.527 65.98 3.495;
 -8.685 -3.781 -8.047 2.337 -3.885 7.225 5.531 -1.915 15.954 -3.249 7.728 -3.374 5.05 10.763 10.551 2.382 -4.395 6.829 -16.964 14.687 -2.315 2.223 -2.802 3.495 58.327]
 td[14,2] = vec(-Inf*fill(1,(25,1)))
 td[14,3] = [6.0; 9.0; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf]
 td[14,4] = 35000
 td[14,5] = 0.665342
 td[14,6] = 0.000001975944
 td[14,7] = 0.00000181746
 td[14,8] = 0.000001083068


 td[1,9] = 0.6054064156442294
 td[2,9] = 0.8279766144831429
 td[3,9] = 0.6536793504453967
 td[4,9] = 0.9431576518489296
 td[5,9] = 0.8279829100356028
 td[6,9] = 0.8413447460685427
 td[7,9] = 0.7612423324911164
 td[8,9] = 0.4741170481094123
 td[9,9] = 0.11351703757272677
 td[10,9] = 0.8103089978617527
 td[11,9] = 0.7593693258920418
 td[12,9] = 0.19640599118630003
 td[13,9] = 1.0
 td[14,9] = 0.6653426686040154
end

@testitem "MvGaussian test" tags=[:Multivariate, :MvGaussian] setup=[GenzTestData] begin
    using AdditionalDistributions, Test, Distributions
    using StableRNGs: StableRNG
    using LinearAlgebra
    td = GenzTestData.td

    # Lower triangular array detector with diag ~ 1 (Genz)
    is_genz_r(r; atol=1e-12) =
        istril(r) &&
        all(abs.(diag(r) .- 1) .<= atol) &&
        all(abs.(triu(r,1)) .<= atol)

    # Robust construction of Σ usable in MvNormal
    function to_cov(r::AbstractMatrix; atol=1e-12)
        if all(abs.(r .- r') .<= atol)
            Σ = r
        elseif is_genz_r(r; atol=atol)
            n = size(r,1)
            Σ = Matrix{eltype(r)}(I, n, n)
            @inbounds for j in 1:n-1, i in j+1:n
                Σ[i,j] = Σ[j,i] = r[i,j]
            end
        else
            Σ = r * r'
        end
        return Symmetric((Σ + Σ')/2)
    end

    for i in 1:14
        r_raw = Float64.(td[i,1])
        a = Float64.(td[i,2])
        b = Float64.(td[i,3])
        m = td[i,4]*10
        pexpected = td[i,5]
        ptol = td[i,6]
        eexpected = td[i,7]
        etol = td[i,8]

        Σ = to_cov(r_raw)

        if !isposdef(Matrix(Σ))
            if i == 6
                @test round(cdf(Normal(), 1.0), digits=6) ≈ pexpected atol=ptol
            else
                @test_skip "Σ no es ≻0 (i=$i)"
            end
            continue
        end

        μ = zeros(size(Σ,1))
        d = AdditionalDistributions.MvGaussian(μ, Matrix(Σ))
        p, e, _ = cdf(d, a, b; m=m, rng=StableRNG(1234), full=true)

        @test abs(p - pexpected) ≤ max(ptol, 5e-6)
        @test isfinite(e) && e ≥ 0
    end

    # sanity check
    r_nowarn = Float64.(td[3,1])
    d_nowarn = AdditionalDistributions.MvGaussian(zeros(3), Matrix(Symmetric(r_nowarn)))
    @test_nowarn cdf(d_nowarn, [0.0,0.0,0.0], [Inf,Inf,Inf]; m=td[3,4], rng=StableRNG(1234))
end


 @testitem "Special cases test" tags=[:Multivariate, :MvGaussian] setup=[GenzTestData] begin
    using AdditionalDistributions 
    using Test, Distributions, StableRNGs, ForwardDiff
    mu    = [1.0, 1.0]
    sigma = [1.0 0.5; 0.5 1.0]
    a     = [0.0, 0.0]
    b     = [Inf, Inf]

    # --- ADAPTATION ---
    # We create our MvGaussian
    d_sp_case1 = AdditionalDistributions.MvGaussian(mu, sigma)
    # We call cdf (full=false by default, returns only the value)
    @test Distributions.cdf(d_sp_case1, a, b; rng = StableRNG(1234)) ≈ 0.7450418725220342 rtol=1e-3
    # --- END ADAPTATION ---

    # --- ADAPTATION ---
    # We create our MvGaussian with zero mean
    d_sp_case2 = AdditionalDistributions.MvGaussian([0,0], sigma)
    @test Distributions.cdf(d_sp_case2, a, b) ≈ 0.33333333333333337 rtol=1e-3
    # --- END ADAPTATION ---
end

@testitem "ForwardDiff test" tags=[:Multivariate, :MvGaussian] setup=[GenzTestData] begin
    using AdditionalDistributions 
    using Test, Distributions, StableRNGs, ForwardDiff, LinearAlgebra
    # NOTE: This is expected to fail if _qf_std is not
    # differentiable, but we'll test it anyway!

    μ = [1., 2., 3.] 
    Σ = [1 0.25 0.2; 0.25 1 0.333333333; 0.2 0.333333333 1]
    ag = float.([-1; -4; -2])
    bg = float.([1; 4; 2])
    
    # --- ADAPTATION ---
    #gf difference from the mean (μ)
    gf(x) = Distributions.cdf(AdditionalDistributions.MvGaussian(x, Σ), ag, bg)
    # --- END ADAPTATION ---
    @test_nowarn ForwardDiff.gradient(gf, [1, 2, 3])

    μ_f2 = [1., 2., 3.] 
    Σ_f2 = [1 0.25 0.2; 0.25 1 0.333333333; 0.2 0.333333333 1]
    ag_f2 = float.([-1; -4; -2])
    bg_f2 = float.([1; 4; 2])

    # --- ADAPTATION ---
    # gf2 difference with respect to the covariance matrix (Σ)
    gf2(x) = begin
    # x has 9 free parameters (3x3). We construct lower triangular L.
    Lraw = LowerTriangular(reshape(x, 3, 3))
    # Enforce positivity on the diagonal:
    L = copy(Lraw)
    @inbounds for i in 1:3
        L[i,i] = exp(Lraw[i,i])   # diag positive
    end
    Σ = Matrix(L*L')              # SPD for construction
    Distributions.cdf(AdditionalDistributions.MvGaussian(μ_f2, Σ), ag_f2, bg_f2)
    end
    # --- END ADAPTATION ---
    @test_nowarn ForwardDiff.gradient(gf2, [1, 0.25, 0.2, 0.25, 1, 0.333333333, 0.2, 0.333333333, 1])

    # --- ADAPTATION ---
    # f difference with respect to the limits of integration (a)
    f(x) = Distributions.cdf(AdditionalDistributions.MvGaussian([0.;0.], reshape([0.1; 0.; 0.; 0.1], 2, 2)), [x; -1.], [1.; 1.])
    # --- END ADAPTATION ---
    @test_nowarn ForwardDiff.derivative(f, -1.)
end


@testitem "MvTStudent - accuracy vs. mvtnorm R (sigma=Σ, delta=0)" tags=[:Multivariate, :MvTStudent] setup=[GenzTestData] begin
    using AdditionalDistributions, Test, Distributions
    using LinearAlgebra, StableRNGs

    # --- References from R: pmvt(sigma = Σ, delta = 0) ---
    const P_REF = Dict{Tuple{Int,Int},Float64}(
        (1,2)=>0.532260496339,(1,5)=>0.572645743266,(1,10)=>0.588975102021,(1,15)=>0.594045002330,(1,20)=>0.597095820258,
        (2,2)=>0.743181905473,(2,5)=>0.791437175658,(2,10)=>0.809332705795,(2,15)=>0.815487454646,(2,20)=>0.818550938537,
        (3,2)=>0.514239757460,(3,5)=>0.588309390085,(3,10)=>0.619278561723,(3,15)=>0.630360540604,(3,20)=>0.636071039301,
        (4,2)=>0.839504817198,(4,5)=>0.895263839008,(4,10)=>0.916362772467,(4,15)=>0.923707841454,(4,20)=>0.927377414052,
        (5,2)=>0.712539450198,(5,5)=>0.761831282499,(5,10)=>0.780381208720,(5,15)=>0.786850690308,(5,20)=>0.790057475252,
        (6,2)=>0.788675381326,(6,5)=>0.818402825939,(6,10)=>0.829564248221,(6,15)=>0.833399182458,(6,20)=>0.835385683306,
        (7,2)=>0.578377449893,(7,5)=>0.674554800823,(7,10)=>0.715447314183,(7,15)=>0.730553924669,(7,20)=>0.737802525677,
        (8,2)=>0.382405869190,(8,5)=>0.432845704968,(8,10)=>0.453071754909,(8,15)=>0.460038226172,(8,20)=>0.463540504453,
        (9,2)=>0.091143335450,(9,5)=>0.102902450219,(9,10)=>0.107880962688,(9,15)=>0.109672304862,(9,20)=>0.110684549005,
        (10,2)=>0.743303951663,(10,5)=>0.783310103926,(10,10)=>0.796965841983,(10,15)=>0.801551879401,(10,20)=>0.803692584793,
        (11,2)=>0.569321301840,(11,5)=>0.667614157446,(11,10)=>0.711210902547,(11,15)=>0.726685019835,(11,20)=>0.734769061854,
        (12,2)=>0.196366820978,(12,5)=>0.196418748911,(12,10)=>0.196400088376,(12,15)=>0.196326430960,(12,20)=>0.196371314153,
        (13,2)=>1.0,(13,5)=>1.0,(13,10)=>1.0,(13,15)=>1.0,(13,20)=>1.0,
        (14,2)=>0.603345346425,(14,5)=>0.638056361651,(14,10)=>0.651253223039,(14,15)=>0.655848256465,(14,20)=>0.658183450650,
    )

    const E_REF = Dict{Tuple{Int,Int},Float64}(
        (1,2)=>0.001730677598,(1,5)=>0.001205320260,(1,10)=>0.001867720246,(1,15)=>0.001604069723,(1,20)=>0.001589323985,
        (2,2)=>0.000145336664,(2,5)=>0.000151409087,(2,10)=>0.000075117045,(2,15)=>0.000055166279,(2,20)=>0.000049306075,
        (3,2)=>0.000014789282,(3,5)=>0.000025868617,(3,10)=>0.000070196320,(3,15)=>0.000041753968,(3,20)=>0.000032147181,
        (4,2)=>0.000056631937,(4,5)=>0.000125912091,(4,10)=>0.000122884671,(4,15)=>0.000122164031,(4,20)=>0.000121694095,
        (5,2)=>0.000094209100,(5,5)=>0.000045971497,(5,10)=>0.000045387064,(5,15)=>0.000029292859,(5,20)=>0.000077222247,
        (6,2)=>0.000015166101,(6,5)=>0.000024939568,(6,10)=>0.000023340550,(6,15)=>0.000055830537,(6,20)=>0.000011352798,
        (7,2)=>0.001282978081,(7,5)=>0.001156497595,(7,10)=>0.000931498746,(7,15)=>0.001139424127,(7,20)=>0.000785414943,
        (8,2)=>0.000061738358,(8,5)=>0.000085466266,(8,10)=>0.000044724015,(8,15)=>0.000073748685,(8,20)=>0.000042237886,
        (9,2)=>0.000112900173,(9,5)=>0.000163249821,(9,10)=>0.000109978431,(9,15)=>0.000235075316,(9,20)=>0.000126493906,
        (10,2)=>0.000552667194,(10,5)=>0.000323907400,(10,10)=>0.000142662800,(10,15)=>0.000068791688,(10,20)=>0.000105562310,
        (11,2)=>0.000432868597,(11,5)=>0.001079054630,(11,10)=>0.000914218573,(11,15)=>0.000834198228,(11,20)=>0.000347980050,
        (12,2)=>0.000097152147,(12,5)=>0.000210888005,(12,10)=>0.000184283910,(12,15)=>0.000132046751,(12,20)=>0.000147748150,
        (13,2)=>0.0,(13,5)=>0.0,(13,10)=>0.0,(13,15)=>0.0,(13,20)=>0.0,
        (14,2)=>1e-15,(14,5)=>1e-15,(14,10)=>1e-15,(14,15)=>1e-15,(14,20)=>1e-15,
    )

    # --- Building Σ consistent with R script (sigma = Σ) ---
    _is_unit_lower_tri(M; atol=1e-12) = istril(M) && all(abs.(diag(M) .- 1) .<= atol) && all(abs.(triu(M,1)) .<= atol)

    function to_cov(M::AbstractMatrix; atol=1e-12)
        if all(abs.(M .- M') .<= atol)
            Σ = Symmetric((M + M')/2)
        elseif _is_unit_lower_tri(M; atol=atol)
            # r: factor (Genz type). Use Σ = r*r'
            F = LowerTriangular(M)
            Σ = Symmetric(F*F')
        else
            # non-symmetric generic case: use as factor
            Σ = Symmetric(M*M')
        end
        return Σ
    end

    rng = StableRNG(123)
    td = GenzTestData.td

    for key in keys(P_REF)
        i, ν = key
        Σ = to_cov(Float64.(td[i,1]))
        if !isposdef(Matrix(Σ))
            @test_skip "Σ not posdef (i=$i) → skip"
            continue
        end

        a = Float64.(td[i,2])
        b = Float64.(td[i,3])
        m = td[i,4]

        # center t with scatter Σ; AdditionalDistributions uses MvTDist(ν, μ, Σ)
        d = AdditionalDistributions.MvTStudent(ν, zeros(size(Σ,1)), Matrix(Σ))
        p̂, ê, _ = cdf(d, a, b; m=m, rng=rng, full=true)

        pref = P_REF[(i,ν)]
        eref = get(E_REF, (i,ν), 0.0)
        atol = max(5e-4, 4e0*ê, 4e0*eref)  # reasonable margin vs pmvt
        @test abs(p̂ - pref) ≤ atol
    end
end
# ===========================================================
# Reproducible reference in R (to regenerate the tables)
# ===========================================================
#
# library(mvtnorm)
# source("genz_testdata.R") # td array identical to that of MvNormalCDF.jl
#
# pmvt_wrap_sigma <- function(lower, upper, df, sigma, maxpts,
# delta = NULL, abseps = 1e-8, releps = 1e-8) {
#if (is.null(delta)) delta <- rep(0, length(upper))
# alg <- GenzBretz(maxpts = maxpts, abseps = abseps, releps = releps)
# res <- pmvt(lower = lower, upper = upper, df = df,
# delta = delta, sigma = sigma, algorithm = alg)
# p <- as.numeric(res)
# err <- attr(res, "error"); if (is.null(err)) err <- NA_real_
# msg <- attr(res, "msg"); if (is.null(msg)) msg ​​<- NA_character_
# list(p = p, error = err, msg = msg)
# }
#
# dfs <- c(2L,5L,10L,15L,20L)
# rows <- list()
# for (i in seq_along(td)) {
# S <- td[[i]]$Sigma; a <- td[[i]]$a; b <- td[[i]]$b; m <- td[[i]]$m
# d <- length(b)
# for (nu in dfs) ​​{
# out <- pmvt_wrap_sigma(lower=a, upper=b, df=nu, sigma=S, maxpts=m)
# rows[[length(rows)+1L]] <- data.frame(i=i, d=d, nu=nu,
# p=out$p, error=out$error, stringsAsFactors=FALSE)
# }
# }
# res <- do.call(rbind, rows)
# res <- res[order(res$i, res$nu), ]
#
# emit_julia_dict(res, name="p_ref_t", key=c("i","nu"), val="p", digits=12)
# emit_julia_dict(res, name="err_ref_t",key=c("i","nu"), val="error", digits=12)
# ===========================================================

@testitem "Multivariate CDFResult API" tags=[:multivariate] begin
    using LinearAlgebra
    using Distributions
    using AdditionalDistributions

    @testset "MvGaussian cdf_result" begin
        d = MvGaussian(zeros(2), Matrix(I, 2, 2))
        a = [-1.0, -1.0]
        b = [1.0, 1.0]

        res = cdf_result(d, a, b; m=20_000)

        @test res isa CDFResult
        @test isfinite(res.value)
        @test isfinite(res.error)
        @test 0.0 <= res.value <= 1.0
        @test res.inform in (0, 1, 2, 3)
        @test res.neval == 20_000
        @test res.algorithm == :mvsort_rqmc

        value, error, inform = res
        @test value == res.value
        @test error == res.error
        @test inform == res.inform
    end

    @testset "MvTStudent cdf_result" begin
        d = MvTStudent(4.5, zeros(2), Matrix(I, 2, 2))
        a = [-1.0, -1.0]
        b = [1.0, 1.0]

        res = cdf_result(d, a, b; m=20_000)

        @test res isa CDFResult
        @test isfinite(res.value)
        @test isfinite(res.error)
        @test 0.0 <= res.value <= 1.0
        @test res.inform in (0, 1, 2, 3)
        @test res.neval == 20_000
        @test res.algorithm == :mvsort_rqmc_t
    end

    @testset "cdf remains scalar-valued" begin
        d = MvGaussian(zeros(2), Matrix(I, 2, 2))
        a = [-1.0, -1.0]
        b = [1.0, 1.0]

        val = cdf(d, a, b; m=20_000)
        res = cdf_result(d, a, b; m=20_000)

        @test val isa Real
        @test 0.0 <= val <= 1.0
        @test 0.0 <= res.value <= 1.0
    end

    @testset "legacy full=true remains supported" begin
        d = MvGaussian(zeros(2), Matrix(I, 2, 2))
        a = [-1.0, -1.0]
        b = [1.0, 1.0]

        value, error, inform = cdf(d, a, b; m=20_000, full=true)

        @test isfinite(value)
        @test isfinite(error)
        @test inform in (0, 1, 2, 3)
    end
end