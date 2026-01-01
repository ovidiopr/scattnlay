// This file is intended to be included inside a namespace where 'hn' is defined
// as an alias to hwy::HWY_NAMESPACE (or similar).

template <typename T_base>
struct HighwayEngine {
  // Highway handles vectors of T_base (double or float)
  using T = T_base;
  using D = hn::ScalableTag<T>;
  using V = hn::Vec<D>;
  using M = hn::Mask<D>;
  using MaskV = M;
  
  // Wrapper for V to support operator overloading
  struct RealV {
      V v;
      RealV(V _v) : v(_v) {}
      RealV() = default;
      operator V() const { return v; }
      
      RealV operator+(const RealV& other) const { return RealV(hn::Add(v, other.v)); }
      RealV operator-(const RealV& other) const { return RealV(hn::Sub(v, other.v)); }
      RealV operator*(const RealV& other) const { return RealV(hn::Mul(v, other.v)); }
      RealV operator/(const RealV& other) const { return RealV(hn::Div(v, other.v)); }
      RealV operator-() const { return RealV(hn::Neg(v)); }
      
      RealV operator*(T s) const { return RealV(hn::Mul(v, hn::Set(D(), s))); }
      friend RealV operator*(T s, const RealV& rv) { return RealV(hn::Mul(hn::Set(D(), s), rv.v)); }
      RealV operator/(T s) const { return RealV(hn::Div(v, hn::Set(D(), s))); }
      
      RealV& operator=(const V& other) { v = other; return *this; }
  };

  template <typename U> using Vector = std::vector<U>;

  static size_t Lanes() { return hn::Lanes(D()); }

  // Represents a batch of complex numbers (Real vectors, Imag vectors)
  struct ComplexV { 
    RealV re; RealV im; 

    friend inline ComplexV operator+(ComplexV a, ComplexV b) {
      return {RealV(hn::Add(a.re.v, b.re.v)), RealV(hn::Add(a.im.v, b.im.v))};
    }
    friend inline ComplexV operator-(ComplexV a, ComplexV b) {
      return {RealV(hn::Sub(a.re.v, b.re.v)), RealV(hn::Sub(a.im.v, b.im.v))};
    }
    friend inline ComplexV operator*(ComplexV a, ComplexV b) {
      return {RealV(hn::Sub(hn::Mul(a.re.v, b.re.v), hn::Mul(a.im.v, b.im.v))),
              RealV(hn::Add(hn::Mul(a.re.v, b.im.v), hn::Mul(a.im.v, b.re.v)))};
    }
    friend inline ComplexV operator/(ComplexV a, ComplexV b) {
      V mag_sq = hn::Add(hn::Mul(b.re.v, b.re.v), hn::Mul(b.im.v, b.im.v));
      V re = hn::Div(hn::Add(hn::Mul(a.re.v, b.re.v), hn::Mul(a.im.v, b.im.v)), mag_sq);
      V im = hn::Div(hn::Sub(hn::Mul(a.im.v, b.re.v), hn::Mul(a.re.v, b.im.v)), mag_sq);
      return {RealV(re), RealV(im)};
    }

    friend inline ComplexV operator+(ComplexV a, RealV b) {
      return {RealV(hn::Add(a.re.v, b.v)), a.im};
    }
    friend inline ComplexV operator+(RealV a, ComplexV b) {
      return {RealV(hn::Add(a.v, b.re.v)), b.im};
    }
    friend inline ComplexV operator-(ComplexV a, RealV b) {
      return {RealV(hn::Sub(a.re.v, b.v)), a.im};
    }
    friend inline ComplexV operator-(RealV a, ComplexV b) {
      return {RealV(hn::Sub(a.v, b.re.v)), RealV(hn::Neg(b.im.v))};
    }
    friend inline ComplexV operator*(ComplexV a, RealV b) {
      return {RealV(hn::Mul(a.re.v, b.v)), RealV(hn::Mul(a.im.v, b.v))};
    }
    friend inline ComplexV operator*(RealV a, ComplexV b) {
      return {RealV(hn::Mul(a.v, b.re.v)), RealV(hn::Mul(a.v, b.im.v))};
    }
    friend inline ComplexV operator/(ComplexV a, RealV b) {
      return {RealV(hn::Div(a.re.v, b.v)), RealV(hn::Div(a.im.v, b.v))};
    }
    friend inline ComplexV operator/(RealV a, ComplexV b) {
      V mag_sq = hn::Add(hn::Mul(b.re.v, b.re.v), hn::Mul(b.im.v, b.im.v));
      V re = hn::Div(hn::Mul(a.v, b.re.v), mag_sq);
      V im = hn::Div(hn::Mul(a.v, hn::Neg(b.im.v)), mag_sq);
      return {RealV(re), RealV(im)};
    }
    friend inline ComplexV operator-(ComplexV a) {
      return {RealV(hn::Neg(a.re.v)), RealV(hn::Neg(a.im.v))};
    }
  };

  static inline RealV sin(RealV v) { return RealV(hn::Sin(D(), v.v)); }
  static inline RealV cos(RealV v) { return RealV(hn::Cos(D(), v.v)); }
  static inline RealV exp(RealV v) { return RealV(hn::Exp(D(), v.v)); }
  static inline RealV tan(RealV v) { 
      return RealV(hn::Div(hn::Sin(D(), v.v), hn::Cos(D(), v.v)));
  }
  static inline RealV sqrt(RealV v) { return RealV(hn::Sqrt(v.v)); }
  static inline RealV abs(RealV v) { return RealV(hn::Abs(v.v)); }
  static inline RealV abs(ComplexV z) {
    return RealV(hn::Sqrt(hn::Add(hn::Mul(z.re.v, z.re.v), hn::Mul(z.im.v, z.im.v))));
  }
  static inline RealV ceil(RealV v) { return RealV(hn::Ceil(v.v)); }
  static inline RealV floor(RealV v) { return RealV(hn::Floor(v.v)); }
  static inline RealV max(RealV a, RealV b) { return RealV(hn::Max(a.v, b.v)); }
  static inline RealV pow(RealV b, RealV e) { 
      return RealV(hn::Exp(D(), hn::Mul(e.v, hn::Log(D(), b.v))));
  }
  static inline RealV log(RealV v) { return RealV(hn::Log(D(), v.v)); }
  static inline ComplexV select(M mask, ComplexV a, ComplexV b) {
    return {RealV(hn::IfThenElse(mask, a.re.v, b.re.v)), RealV(hn::IfThenElse(mask, a.im.v, b.im.v))};
  }

  static inline RealV select(M mask, RealV a, RealV b) {
    return RealV(hn::IfThenElse(mask, a.v, b.v));
  }

  static inline M gt(RealV a, RealV b) { return hn::Gt(a.v, b.v); }
  static inline M lt(RealV a, RealV b) { return hn::Lt(a.v, b.v); }
  static inline M ge(RealV a, RealV b) { return hn::Ge(a.v, b.v); }
  static inline M le(RealV a, RealV b) { return hn::Le(a.v, b.v); }
  static inline M eq(RealV a, RealV b) { return hn::Eq(a.v, b.v); }
  static inline M neq(RealV a, RealV b) { return hn::Ne(a.v, b.v); }
  static inline M lor(M a, M b) { return hn::Or(a, b); }
  static inline M land(M a, M b) { return hn::And(a, b); }

  static inline RealV set(T val) { return RealV(hn::Set(D(), val)); }
  
  static inline RealV sign(RealV v) {
    auto zero = hn::Zero(D());
    auto one = hn::Set(D(), 1.0);
    auto minus_one = hn::Set(D(), -1.0);
    auto mask = hn::Gt(v.v, zero);
    return RealV(hn::IfThenElse(mask, one, minus_one));
  }



  static inline RealV get_real(ComplexV z) { return z.re; }
  static inline RealV get_imag(ComplexV z) { return z.im; }
  static inline ComplexV make_complex(RealV re, RealV im) { return {re, im}; }

  static inline ComplexV load_interleaved(const T* data) {
    D d;
    V re, im;
    hn::LoadInterleaved2(d, data, re, im);
    return {RealV(re), RealV(im)};
  }

  static inline void store_interleaved(ComplexV z, T* data) {
    D d;
    hn::StoreInterleaved2(z.re.v, z.im.v, d, data);
  }

  static inline void store(ComplexV z, std::complex<T>* ptr) {
    store_interleaved(z, reinterpret_cast<T*>(ptr));
  }

  static inline ComplexV load(const std::complex<T>* ptr) {
    return load_interleaved(reinterpret_cast<const T*>(ptr));
  }

  static inline RealV load(const T* ptr) {
    return RealV(hn::LoadU(D(), ptr));
  }

  static inline void store(RealV v, T* ptr) {
    hn::StoreU(v.v, D(), ptr);
  }

  static inline T reduce_max(RealV v) {
    size_t lanes = hn::Lanes(D());
    std::vector<T> buf(lanes);
    hn::StoreU(v.v, D(), buf.data());
    T max_val = buf[0];
    for(size_t i=1; i<lanes; ++i) {
      if(buf[i] > max_val) max_val = buf[i];
    }
    return max_val;
  }

  static inline ComplexV log(ComplexV z) {
    V r = hn::Sqrt(hn::Add(hn::Mul(z.re.v, z.re.v), hn::Mul(z.im.v, z.im.v)));
    V phi = hn::Atan2(D(), z.im.v, z.re.v);
    return {RealV(hn::Log(D(), r)), RealV(phi)};
  }

  static inline ComplexV sqrt(ComplexV z) {
    V r = hn::Sqrt(hn::Add(hn::Mul(z.re.v, z.re.v), hn::Mul(z.im.v, z.im.v)));
    V u = hn::Sqrt(hn::Mul(hn::Add(r, z.re.v), hn::Set(D(), 0.5)));
    V v_abs = hn::Sqrt(hn::Mul(hn::Sub(r, z.re.v), hn::Set(D(), 0.5)));
    V zero = hn::Zero(D());
    auto mask = hn::Ge(z.im.v, zero);
    V v = hn::IfThenElse(mask, v_abs, hn::Neg(v_abs));
    return {RealV(u), RealV(v)};
  }
};
