#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// Revealed preference relations for a finite price/quantity data set.
//
// e[t]      = p[t, ] . q[t, ]                 expenditure at t
// c[t, s]   = p[t, ] . q[s, ]                 cost of bundle s at prices t
// R0[t, s]  = (efficiency * e[t] >= c[t, s])  direct revealed preference
// P0[t, s]  = (efficiency * e[t] >  c[t, s])  strict direct revealed preference
// R         = transitive closure of R0        (Warshall 1962)
//
// [[Rcpp::export]]
List rp_relations(NumericMatrix p, NumericMatrix q, double efficiency) {
  int tt = p.nrow(), n = p.ncol();

  std::vector<double> e(tt, 0.0);
  for (int t = 0; t < tt; t++) {
    double s = 0.0;
    for (int j = 0; j < n; j++) s += p(t, j) * q(t, j);
    e[t] = s;
  }

  LogicalMatrix R0(tt, tt), P0(tt, tt), R(tt, tt);
  for (int t = 0; t < tt; t++) {
    double lhs = efficiency * e[t];
    for (int s = 0; s < tt; s++) {
      double cost = 0.0;
      for (int j = 0; j < n; j++) cost += p(t, j) * q(s, j);
      R0(t, s) = (lhs >= cost);
      P0(t, s) = (lhs >  cost);
      R(t, s)  = R0(t, s);
    }
  }

  // Warshall transitive closure.
  for (int k = 0; k < tt; k++) {
    for (int t = 0; t < tt; t++) {
      if (!R(t, k)) continue;
      for (int s = 0; s < tt; s++) {
        if (R(k, s)) R(t, s) = true;
      }
    }
  }

  return List::create(_["R0"] = R0, _["P0"] = P0, _["R"] = R);
}

// Scan for pairs (t, s) with A[t, s] and B[s, t]. Note the transposed index on
// B: an axiom violation pairs a preference from t to s against a strict
// preference running back from s to t.
//
// [[Rcpp::export]]
IntegerMatrix scan_violations(LogicalMatrix A, LogicalMatrix B,
                              bool skip_diagonal) {
  int tt = A.nrow();
  std::vector<int> vt, vs;
  for (int t = 0; t < tt; t++) {
    for (int s = 0; s < tt; s++) {
      if (skip_diagonal && t == s) continue;
      if (A(t, s) && B(s, t)) {
        vt.push_back(t + 1);
        vs.push_back(s + 1);
      }
    }
  }
  IntegerMatrix out(vt.size(), 2);
  for (size_t i = 0; i < vt.size(); i++) {
    out(i, 0) = vt[i];
    out(i, 1) = vs[i];
  }
  colnames(out) = CharacterVector::create("t", "s");
  return out;
}
