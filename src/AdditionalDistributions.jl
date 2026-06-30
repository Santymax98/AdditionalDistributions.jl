module AdditionalDistributions
    
    using Distributions
    import Distributions: @check_args, @distr_support, @inline, params, partype
    import Roots
    import SpecialFunctions
    import Random
    import Base: @propagate_inbounds
    import Statistics
    import StatsBase
    import HypergeometricFunctions
    import QuadGK
    import LambertW
    import LogExpFunctions
    import LinearAlgebra
    import Primes

    const normal_dist = Distributions.Normal()
    const sqrt2π = sqrt(2.0 * π)


    include("utils.jl")
    
    #discrete Distributions
    include("discrete/BetaNegBinomial.jl")
    include("discrete/Borel.jl")
    include("discrete/Conway.jl")
    include("discrete/Delaporte.jl")
    include("discrete/FlorySchulz.jl")
    include("discrete/GaussKuzmin.jl")
    include("discrete/Logarithmic.jl")
    include("discrete/Rademacher.jl")
    include("discrete/Yule.jl")
    include("discrete/Zeta.jl")
    include("discrete/ZIB.jl")
    include("discrete/ZINB.jl")
    include("discrete/ZIP.jl")
    include("discrete/Zipf.jl")


    export 
        BetaNegBinomial,
        Borel,
        Conway,
        Delaporte,
        FlorySchulz,
        GaussKuzmin,
        Logarithmic,
        Rademacher,
        Yule,
        Zeta,
        ZIB,
        ZINB,
        ZIP,
        Zipf

    #Continuous Distributions
    include("continuous/Alpha.jl")
    include("continuous/Argus.jl")
    include("continuous/Benini.jl")
    include("continuous/Benktander_Type1.jl")
    include("continuous/Benktander_Type2.jl")
    include("continuous/Bhattacharjee.jl")
    include("continuous/BirnbaumSaunders.jl")
    include("continuous/Bradford.jl")
    include("continuous/Burr.jl")
    include("continuous/CrystalBall.jl")
    include("continuous/Dagum.jl")
    include("continuous/Gompertz.jl")
    include("continuous/IrwinHall.jl")
    include("continuous/Lomax.jl")
    include("continuous/Maxwell.jl")
    include("continuous/Nakagami.jl")
    include("continuous/PERT.jl")

    export
        Alpha,
        Argus,
        Benini,
        Benktander_Type1,
        Benktander_Type2,
        Bhattacharjee,
        BirnbaumSaunders,
        Bradford,
        Burr,
        CrystalBall,
        Dagum,
        Gompertz,
        IrwinHall,
        Lomax,
        Maxwell,
        Nakagami,
        PERT
    
        include("multivariate/CDFResult.jl")
        include("multivariate/MvGaussian.jl")
        include("multivariate/MvStudent.jl")

        export 
            MvGaussian,
            MvTStudent,
            CDFResult,
            cdf_result
        
end
