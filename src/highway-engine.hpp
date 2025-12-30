// This file is intended to be included inside a namespace where 'hn' is defined
// as an alias to hwy::HWY_NAMESPACE (or similar).

template <typename T_base>
struct HighwayEngine {
  // Highway handles vectors of T_base (double or float)
  using T = T_base;
  using D = hn::ScalableTag<T>;
  using V = hn::Vec<D>;
  using M = hn::Mask<D>;
  using RealV = V;

  static size_t Lanes() { return hn::Lanes(D()); }

  // Represents a batch of complex numbers (Real vectors, Imag vectors)
  struct ComplexV { V re; V im; };

  static inline V sin(V v) { return hn::Sin(D(), v); }
  static inline V cos(V v) { return hn::Cos(D(), v); }
  static inline V exp(V v) { return hn::Exp(D(), v); }
  static inline V tan(V v) { 
      return hn::Div(hn::Sin(D(), v), hn::Cos(D(), v));
  }
  static inline V sqrt(V v) { return hn::Sqrt(v); }
  static inline V abs(V v) { return hn::Abs(v); }
  static inline V abs(ComplexV z) {
    return hn::Sqrt(hn::Add(hn::Mul(z.re, z.re), hn::Mul(z.im, z.im)));
  }
  static inline V ceil(V v) { return hn::Ceil(v); }
  static inline V floor(V v) { return hn::Floor(v); }
  static inline V max(V a, V b) { return hn::Max(a, b); }
  static inline V pow(V b, V e) { 
      return hn::Exp(D(), hn::Mul(e, hn::Log(D(), b)));
  }
  static inline V log(V v) { return hn::Log(D(), v); }
  static inline ComplexV select(M mask, ComplexV a, ComplexV b) {
    return {hn::IfThenElse(mask, a.re, b.re), hn::IfThenElse(mask, a.im, b.im)};
  }

  static inline V select(M mask, V a, V b) {
    return hn::IfThenElse(mask, a, b);
  }

  static inline M gt(V a, V b) { return hn::Gt(a, b); }
  static inline M lt(V a, V b) { return hn::Lt(a, b); }
  static inline M ge(V a, V b) { return hn::Ge(a, b); }
  static inline M le(V a, V b) { return hn::Le(a, b); }
  static inline M neq(V a, V b) { return hn::Ne(a, b); }
  static inline M lor(M a, M b) { return hn::Or(a, b); }
  static inline M land(M a, M b) { return hn::And(a, b); }
  static inline V add(V a, V b) { return hn::Add(a, b); }
  static inline V sub(V a, V b) { return hn::Sub(a, b); }
  static inline V mul(V a, V b) { return hn::Mul(a, b); }
  static inline V div(V a, V b) { return hn::Div(a, b); }
  static inline V set(T val) { return hn::Set(D(), val); }
  
  static inline V sign(V v) {
    auto zero = hn::Zero(D());
    auto one = hn::Set(D(), 1.0);
    auto minus_one = hn::Set(D(), -1.0);
    auto mask = hn::Gt(v, zero);
    return hn::IfThenElse(mask, one, minus_one);
  }

  // Complex arithmetic
  static inline ComplexV add(ComplexV a, ComplexV b) {
    return {hn::Add(a.re, b.re), hn::Add(a.im, b.im)};
  }

  static inline ComplexV sub(ComplexV a, ComplexV b) {
    return {hn::Sub(a.re, b.re), hn::Sub(a.im, b.im)};
  }

  static inline ComplexV mul(ComplexV a, ComplexV b) {
    // (a.re*b.re - a.im*b.im) + i(a.re*b.im + a.im*b.re)
    return {hn::Sub(hn::Mul(a.re, b.re), hn::Mul(a.im, b.im)),
            hn::Add(hn::Mul(a.re, b.im), hn::Mul(a.im, b.re))};
  }
  
  static inline ComplexV div(ComplexV a, ComplexV b) {
    V mag_sq = hn::Add(hn::Mul(b.re, b.re), hn::Mul(b.im, b.im));
    V re = hn::Div(hn::Add(hn::Mul(a.re, b.re), hn::Mul(a.im, b.im)), mag_sq);
    V im = hn::Div(hn::Sub(hn::Mul(a.im, b.re), hn::Mul(a.re, b.im)), mag_sq);
    return {re, im};
  }

  static inline V get_real(ComplexV z) { return z.re; }
  static inline V get_imag(ComplexV z) { return z.im; }
  static inline ComplexV make_complex(V re, V im) { return {re, im}; }

  static inline ComplexV load_interleaved(const T* data) {
    D d;
    V re, im;
    hn::LoadInterleaved2(d, data, re, im);
    return {re, im};
  }

  static inline void store_interleaved(ComplexV z, T* data) {
    D d;
    hn::StoreInterleaved2(z.re, z.im, d, data);
  }

  static inline void store(ComplexV z, std::complex<T>* ptr) {
    store_interleaved(z, reinterpret_cast<T*>(ptr));
  }

  static inline ComplexV load(const std::complex<T>* ptr) {
    return load_interleaved(reinterpret_cast<const T*>(ptr));
  }

  static inline V load(const T* ptr) {
    return hn::LoadU(D(), ptr);
  }

  static inline void store(V v, T* ptr) {
    hn::StoreU(v, D(), ptr);
  }

  static inline T reduce_max(V v) {
    size_t lanes = hn::Lanes(D());
    std::vector<T> buf(lanes);
    hn::StoreU(v, D(), buf.data());
    T max_val = buf[0];
    for(size_t i=1; i<lanes; ++i) {
      if(buf[i] > max_val) max_val = buf[i];
    }
    return max_val;
  }

  static inline ComplexV log(ComplexV z) {
    V r = hn::Sqrt(hn::Add(hn::Mul(z.re, z.re), hn::Mul(z.im, z.im)));
    V phi = hn::Atan2(D(), z.im, z.re);
    return {hn::Log(D(), r), phi};
  }

  static inline ComplexV sqrt(ComplexV z) {
    V r = hn::Sqrt(hn::Add(hn::Mul(z.re, z.re), hn::Mul(z.im, z.im)));
    V u = hn::Sqrt(hn::Mul(hn::Add(r, z.re), hn::Set(D(), 0.5)));
    V v_abs = hn::Sqrt(hn::Mul(hn::Sub(r, z.re), hn::Set(D(), 0.5)));
    V zero = hn::Zero(D());
    auto mask = hn::Ge(z.im, zero);
    V v = hn::IfThenElse(mask, v_abs, hn::Neg(v_abs));
    return {u, v};
  }
};
