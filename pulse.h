class pulse{
  public:
  rl E0, A0, omega, pulselength;
  rl Ez(rl t);
  rl Az(rl t);

  pulse();
  pulse(rl in_E0, rl in_omega);
  pulse(rl in_E0, rl in_omega, rl in_pulselength);
 
// private:
//  static const rl defaultE0=0.1;
//  static const rl defaultomega=0.0565;
//  static const rl defaultpulselength=1000;
};
