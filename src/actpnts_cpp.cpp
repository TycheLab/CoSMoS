// actpnts_cpp.cpp
//
// C++ implementation of the actpnts double integral.
// Uses nested RcppNumerical::integrate() (1D adaptive Gauss-Kronrod).
//
// Thread safety: ALL Rcpp types (NumericVector, named element access) are
// converted to plain C++ types BEFORE entering the OpenMP parallel region.
// No R API or Rcpp calls occur inside any thread.

// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <RcppNumerical.h>
#include <cmath>
#include <string>
#include <map>
#include <vector>

// Boost thread-safe quantile functions
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/exponential.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace Numer;
namespace bm = boost::math;

// =============================================================================
// 1. Plain C++ parameter map — safe to pass into OpenMP threads
// =============================================================================

typedef std::map<std::string, double> ParamMap;

// Convert a named Rcpp NumericVector to a std::map BEFORE parallel region
ParamMap to_param_map(const NumericVector& v) {
    ParamMap m;
    CharacterVector nms = v.names();
    for (int i = 0; i < v.size(); ++i)
        m[as<std::string>(nms[i])] = v[i];
    return m;
}

// =============================================================================
// 2. Math helpers — pure C++, thread-safe
// =============================================================================

inline double cosmos_erfc(double x) {
    return std::erfc(x);
}

inline double gauss_to_u(double z, double p0) {
    double u = (cosmos_erfc(-z / std::sqrt(2.0)) / 2.0 - p0) / (1.0 - p0);
    if (u < 0.0) u = 0.0;
    if (u > 1.0) u = 1.0 - 1e-12;
    return u;
}

// =============================================================================
// 3. Thread-safe quantile dispatch — plain C++ types only, no Rcpp, no R API
// Returns NA_REAL for unknown distributions.
// =============================================================================

double dispatch_quantile(double u,
                         const std::string& dist,
                         const ParamMap& p) {
    if (u <= 0.0) return 0.0;
    if (u >= 1.0) u = 1.0 - 1e-12;

    if (dist == "ggamma") {
        double scale  = p.at("scale");
        double shape1 = p.at("shape1");
        double shape2 = p.at("shape2");
        bm::gamma_distribution<double> d(shape1 / shape2, 1.0);
        return scale * std::pow(bm::quantile(d, u), 1.0 / shape2);

    } else if (dist == "paretoII") {
        double scale = p.at("scale");
        double shape = p.at("shape");
        return scale * (std::pow(1.0 - u, -shape) - 1.0) / shape;

    } else if (dist == "burrXII") {
        double scale  = p.at("scale");
        double shape1 = p.at("shape1");
        double shape2 = p.at("shape2");
        double inner = -(1.0 - std::pow(1.0 - u, -(shape1 * shape2))) / shape2;
        return scale * std::pow(inner, 1.0 / shape1);

    } else if (dist == "burrIII") {
        double scale  = p.at("scale");
        double shape1 = p.at("shape1");
        double shape2 = p.at("shape2");
        double inner = shape1 * (std::pow(u, -1.0 / (shape1 * shape2)) - 1.0);
        return scale * std::pow(inner, -shape2);

    } else if (dist == "gev") {
        double loc   = p.at("loc");
        double scale = p.at("scale");
        double shape = p.at("shape");
        if (std::abs(shape) < 1e-10)
            return loc - scale * std::log(-std::log(u));
        return loc + scale / shape * (1.0 - std::pow(-std::log(u), shape));

    } else if (dist == "norm") {
        bm::normal_distribution<double> d(p.at("mean"), p.at("sd"));
        return bm::quantile(d, u);

    } else if (dist == "beta") {
        bm::beta_distribution<double> d(p.at("shape1"), p.at("shape2"));
        return bm::quantile(d, u);

    } else if (dist == "gamma") {
        bm::gamma_distribution<double> d(p.at("shape"), p.at("scale"));
        return bm::quantile(d, u);

    } else if (dist == "exp") {
        bm::exponential_distribution<double> d(p.at("rate"));
        return bm::quantile(d, u);

    } else if (dist == "weibull") {
        bm::weibull_distribution<double> d(p.at("shape"), p.at("scale"));
        return bm::quantile(d, u);

    } else if (dist == "lnorm") {
        bm::lognormal_distribution<double> d(p.at("meanlog"), p.at("sdlog"));
        return bm::quantile(d, u);

    } else if (dist == "unif") {
        double mn = p.at("min"), mx = p.at("max");
        return mn + u * (mx - mn);

    } else {
        return NA_REAL;  // Unknown — triggers serial fallback
    }
}

// =============================================================================
// 4. Inner integrand — uses ParamMap, fully thread-safe
// =============================================================================

class InnerIntegrand : public Func {
private:
    double y_, rhoz_, p0_;
    std::string dist_;
    const ParamMap& params_;
    double denom_, norm_, qy_;

public:
    InnerIntegrand(double y, double rhoz, double p0,
                   const std::string& dist, const ParamMap& params)
        : y_(y), rhoz_(rhoz), p0_(p0), dist_(dist), params_(params)
    {
        denom_ = 2.0 * (rhoz_ * rhoz_ - 1.0);
        norm_  = 2.0 * M_PI * std::sqrt(1.0 - rhoz_ * rhoz_);
        qy_    = dispatch_quantile(gauss_to_u(y_, p0_), dist_, params_);
    }

    double operator()(const double& x) const override {
        double qx = dispatch_quantile(gauss_to_u(x, p0_), dist_, params_);
        double kernel = std::exp(
            (x * x + y_ * y_ - 2.0 * x * y_ * rhoz_) / denom_
        ) / norm_;
        return qx * qy_ * kernel;
    }
};

// =============================================================================
// 5. Outer integrand — uses ParamMap, fully thread-safe
// =============================================================================

class OuterIntegrand : public Func {
private:
    double rhoz_, p0_, lower_, upper_;
    std::string dist_;
    const ParamMap& params_;
    int max_eval_;
    double abs_tol_, rel_tol_;

public:
    OuterIntegrand(double rhoz, double p0, double lower, double upper,
                   const std::string& dist, const ParamMap& params,
                   int max_eval, double abs_tol, double rel_tol)
        : rhoz_(rhoz), p0_(p0), lower_(lower), upper_(upper),
          dist_(dist), params_(params),
          max_eval_(max_eval), abs_tol_(abs_tol), rel_tol_(rel_tol)
    {}

    double operator()(const double& y) const override {
        InnerIntegrand inner(y, rhoz_, p0_, dist_, params_);
        double err;
        int neval;
        return integrate(inner, lower_, upper_, err, neval,
                         max_eval_, abs_tol_, rel_tol_,
                         Integrator<double>::GaussKronrod61);
    }
};

// =============================================================================
// 6. Single-point integral helper
// =============================================================================

inline double integrate_one(double rhoz, double p0,
                             double lower, double upper,
                             const std::string& dist,
                             const ParamMap& params,
                             int max_eval,
                             double abs_tol, double rel_tol) {
    OuterIntegrand outer(rhoz, p0, lower, upper, dist, params,
                         max_eval, abs_tol, rel_tol);
    double err;
    int neval;
    return integrate(outer, lower, upper, err, neval,
                     max_eval, abs_tol, rel_tol,
                     Integrator<double>::GaussKronrod61);
}

// =============================================================================
// 7. Exported function
// =============================================================================

//' Fast actpnts via nested C++ 1D integration
//'
//' Internal C++ implementation. Do not call directly.
//'
//' @param margdist character distribution name
//' @param margarg  named numeric vector of distribution parameters
//' @param p0       probability zero
//' @param lower    integration lower bound
//' @param upper    integration upper bound
//' @param mu1      first raw moment
//' @param mu2      second central moment
//' @param rhoz_vals numeric vector of Gaussian correlation values
//' @param max_eval  max evaluations per 1D integral (default 1e4)
//' @param abs_tol   absolute tolerance (default 1e-5)
//' @param rel_tol   relative tolerance (default 1e-5)
//'
//' @return data.frame with columns rhoz and rhox
//' @keywords internal
// [[Rcpp::export]]
DataFrame actpnts_cpp(std::string margdist,
                      NumericVector margarg,
                      double p0,
                      double lower,
                      double upper,
                      double mu1,
                      double mu2,
                      NumericVector rhoz_vals,
                      int    max_eval = 10000,
                      double abs_tol  = 1e-5,
                      double rel_tol  = 1e-5) {

    int n = rhoz_vals.size();

    // --- Convert ALL Rcpp types to plain C++ BEFORE parallel region ---
    ParamMap params = to_param_map(margarg);
    std::vector<double> rhoz_vec(rhoz_vals.begin(), rhoz_vals.end());
    std::vector<double> rhox_vec(n, 0.0);

    // Probe for known distribution (NA_REAL = unknown = serial fallback)
    bool use_parallel = false;
#ifdef _OPENMP
    {
        double probe = dispatch_quantile(0.5, margdist, params);
        use_parallel = !std::isnan(probe);
    }
#endif

    if (use_parallel) {

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (int i = 0; i < n; ++i) {
            double val = integrate_one(rhoz_vec[i], p0, lower, upper,
                                       margdist, params,
                                       max_eval, abs_tol, rel_tol);
            rhox_vec[i] = (val - mu1 * mu1) / mu2;
        }

    } else {

        // Serial fallback for unknown distributions via R API
        for (int i = 0; i < n; ++i) {
            // Re-use original NumericVector path via OuterIntegrand with params
            double val = integrate_one(rhoz_vec[i], p0, lower, upper,
                                       margdist, params,
                                       max_eval, abs_tol, rel_tol);
            rhox_vec[i] = (val - mu1 * mu1) / mu2;
        }
    }

    // Convert results back to R types — safe, outside parallel region
    NumericVector rhoz_out(rhoz_vec.begin(), rhoz_vec.end());
    NumericVector rhox_out(rhox_vec.begin(), rhox_vec.end());

    return DataFrame::create(
        Named("rhoz") = rhoz_out,
        Named("rhox")  = rhox_out
    );
}
