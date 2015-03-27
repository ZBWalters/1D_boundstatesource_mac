class basis{
 public:
  int order;
  cmplx* omat;
  cmplx* nablamat;
  cmplx* nablasqmat;

//basis is a virtual class, so it has no constructors or destructors of
//its own.
//The motivation for defining an abstract basis class is that it will often
//prove useful to define a wavefunction using different primitive bases in
//different regions -- for example, in infinite range ecs, to use laguerre
//polynomials for the end regions, but legendre polynomials for the interior
//regions.  Another possibility would be to use orthogonal polynomials with
//different weighting functions in regions which contain a singularity, so that
//the calculated wavefunction remains well behaved.  For this reason, it is
//useful to define sets of primitive bases -- such as orthogonal polynomials --
//and sets of derived bases, which consist of linear combinations of primitive
//basis functions. (If desired, derived bases can be derived from linear
//combinations of other derived basis functions).


//  basis();
//  basis(int order_in);
//  basis(int order_in, cmplx* omat_in, cmplx* nablamat_in,
//	cmplx* nablasqmat_in);
//  basis(basis* inpbasis);
//  
//  ~basis();

  virtual cmplx* evalbasisfuncs(rl x)=0;
  virtual cmplx* funcmatrix(cmplx (*func)(rl x), rl (*ytox)(rl y, basis* bas1, basis* bas2),basis* dbasis)=0;
  virtual rl** integrationtable(int glorder)=0;
  virtual rl mapy2range(rl y, rl x1, rl x2)=0;
  virtual rl dy2dx(rl x1, rl x2)=0;

  virtual rl mapfromchildrange(rl childx, rl childxmin, rl childxmax)=0;
  virtual rl maptochildrange(rl y, rl childxmin, rl childxmax)=0;

  void printomat();
  void printnablamat();
  void printnablasqmat();

  cmplx* oinvpsi(cmplx* psi);
};