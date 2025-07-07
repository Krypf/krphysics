module RelativityTools

include("Metrics.jl")
export SchwarzschildMetric, polar, myreplace, metric_value, metric_test
include("Accels.jl")
export AccelSchwarzschild, AccelKerr
include("RKS.jl")
export Rungekutta, newk, SCprodsum, RK4, RKF45, RK8

end
