class potential{
 public:
  potential();
  potential(rl q, rl rc); 
  rl V(rl x);
  cmplx Vecs(rl x, rl recs, rl thetaecs);
  rl rc;
  rl q;
 private:
};
