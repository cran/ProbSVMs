#ifndef TUNESAUX_H
#define TUNESAUX_H

void cnvVct2armamat(const arma::vec& v, arma::mat& M);

struct srchres {
  arma::mat alpha;
  double fval;
  unsigned iter;
  double hitrate;
};

#endif
