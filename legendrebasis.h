//#include "./basis.h"

class legendrebasis : public basis{
 public:
  //methods
  //legendrebasis(legendrebasis* inpbasis);
  legendrebasis(int order_in);
  ~legendrebasis();
  cmplx* funcmatrix(cmplx (*func)(rl x), rl (*ytox)(rl y, basis* bas1, basis* bas2),basis* dbasis);
  cmplx* evalbasisfuncs(rl y);
  rl** integrationtable(int glorder);
  rl mapy2range(rl y, rl x1, rl x2);
  rl mapfromchildrange(rl childx, rl childxmin, rl childxmax);
  rl maptochildrange(rl y, rl childxmin, rl childxmax);
  rl dy2dx(rl x1, rl x2);
};
