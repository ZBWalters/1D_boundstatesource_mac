//#include "./basis.h"

//basis of laguerre polynomials times exponentials: L_n(x)*exp(kappa*x)
//if kappa>0, basis extends from -infinity <y<= 0,
//if kappa <0, basis extends from 0<=y<infinity
class laguerrebasis : public basis{
 public:
  laguerrebasis(int order_in, rl kappa_in);
  ~laguerrebasis();
  rl kappa;

  cmplx* funcmatrix(cmplx (*func)(rl x), rl (*ytox)(rl y, basis* bas1, basis* bas2),basis* dbasis);
  cmplx* evalbasisfuncs(rl y);
  cmplx* evallaguerrefuncs(rl y);
  rl** integrationtable(int glorder);
  rl mapy2range(rl y, rl x1, rl x2);
  rl mapfromchildrange(rl childx, rl childxmin, rl childxmax);
  rl maptochildrange(rl y, rl childxmin, rl childxmax);
  rl dy2dx(rl x1, rl x2);
};
