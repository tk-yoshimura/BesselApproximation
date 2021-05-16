using MultiPrecision;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BesselApproximation {

    public static class MultiPrecisionSandbox<N> where N : struct, IConstant {

        public static MultiPrecision<N> BesselJ(MultiPrecision<N> nu, MultiPrecision<N> x) {
            if (MultiPrecision<N>.Abs(nu) > 64) {
                throw new ArgumentOutOfRangeException(nameof(nu));
            }
            if (nu.IsNaN || x.IsNaN) {
                return MultiPrecision<N>.NaN;
            }

            if (x.Sign == Sign.Minus) {
                if (nu != MultiPrecision<N>.Truncate(nu)) {
                    return MultiPrecision<N>.NaN;
                }

                long n = (long)nu;
                return ((n & 1L) == 0) ? BesselJ(nu, MultiPrecision<N>.Abs(x)) : -BesselJ(nu, MultiPrecision<N>.Abs(x));
            }
            if (x.IsZero && nu.IsZero) {
                return 1;
            }
            if (nu.Sign == Sign.Minus && nu == MultiPrecision<N>.Truncate(nu)) {
                long n = (long)nu;
                return ((n & 1L) == 0) ? BesselJ(MultiPrecision<N>.Abs(nu), x) : -BesselJ(MultiPrecision<N>.Abs(nu), x);
            }

            if (nu - MultiPrecision<N>.Point5 == MultiPrecision<N>.Truncate(nu)) {
                long n = (long)MultiPrecision<N>.Truncate(nu);

                if (n >= -2 || n < 2) {
                    MultiPrecision<Plus1<N>> x_ex = x.Convert<Plus1<N>>();
                    MultiPrecision<Plus1<N>> envelope = MultiPrecision<Plus1<N>>.Sqrt(2 / (MultiPrecision<Plus1<N>>.PI * x_ex));

                    if (n == -2) {
                        return -(envelope * (MultiPrecision<Plus1<N>>.Cos(x_ex) / x_ex + MultiPrecision<Plus1<N>>.Sin(x_ex))).Convert<N>();
                    }
                    if (n == -1) {
                        return (envelope * MultiPrecision<Plus1<N>>.Cos(x_ex)).Convert<N>();
                    }
                    if (n == 0) {
                        return (envelope * MultiPrecision<Plus1<N>>.Sin(x_ex)).Convert<N>();
                    }
                    if (n == 1) {
                        return (envelope * (MultiPrecision<Plus1<N>>.Sin(x_ex) / x_ex - MultiPrecision<Plus1<N>>.Cos(x_ex))).Convert<N>();
                    }
                }
            }

            if (MultiPrecision<N>.Length < 8) {
                return MultiPrecisionSandbox<Pow2.N8>.BesselJ(nu.Convert<Pow2.N8>(), x.Convert<Pow2.N8>()).Convert<N>();
            }

            if (x < Consts.BesselJY.ApproxThreshold) {
                return BesselJNearZero(nu, x).Convert<N>();
            }
            else {
                return BesselJLimit(nu, x).Convert<N>();
            }
        }

        public static MultiPrecision<N> BesselY(MultiPrecision<N> nu, MultiPrecision<N> x) {
            if (MultiPrecision<N>.Abs(nu) > 64) {
                throw new ArgumentOutOfRangeException(nameof(nu));
            }
            if (nu.IsNaN || x.IsNaN) {
                return MultiPrecision<N>.NaN;
            }

            if (x.Sign == Sign.Minus) {
                if (nu != MultiPrecision<N>.Truncate(nu)) {
                    return MultiPrecision<N>.NaN;
                }

                long n = (long)nu;
                return ((n & 1L) == 0) ? BesselY(nu, MultiPrecision<N>.Abs(x)) : -BesselY(nu, MultiPrecision<N>.Abs(x));
            }
            if (x.IsZero) {
                return MultiPrecision<N>.NegativeInfinity;
            }
            if (nu.Sign == Sign.Minus && nu == MultiPrecision<N>.Truncate(nu)) {
                long n = (long)nu;
                return ((n & 1L) == 0) ? BesselY(MultiPrecision<N>.Abs(nu), x) : -BesselY(MultiPrecision<N>.Abs(nu), x);
            }

            if (nu - MultiPrecision<N>.Point5 == MultiPrecision<N>.Truncate(nu)) {
                long n = (long)MultiPrecision<N>.Truncate(nu);

                if (n >= -2 || n < 2) {
                    MultiPrecision<Plus1<N>> x_ex = x.Convert<Plus1<N>>();
                    MultiPrecision<Plus1<N>> envelope = MultiPrecision<Plus1<N>>.Sqrt(2 / (MultiPrecision<Plus1<N>>.PI * x_ex));

                    if (n == -2) {
                        return -(envelope * (MultiPrecision<Plus1<N>>.Sin(x_ex) / x_ex - MultiPrecision<Plus1<N>>.Cos(x_ex))).Convert<N>();
                    }
                    if (n == -1) {
                        return (envelope * MultiPrecision<Plus1<N>>.Sin(x_ex)).Convert<N>();
                    }
                    if (n == 0) {
                        return -(envelope * MultiPrecision<Plus1<N>>.Cos(x_ex)).Convert<N>();
                    }
                    if (n == 1) {
                        return -(envelope * (MultiPrecision<Plus1<N>>.Cos(x_ex) / x_ex + MultiPrecision<Plus1<N>>.Sin(x_ex))).Convert<N>();
                    }
                }
            }

            if (MultiPrecision<N>.Length < 8) {
                return MultiPrecisionSandbox<Pow2.N8>.BesselY(nu.Convert<Pow2.N8>(), x.Convert<Pow2.N8>()).Convert<N>();
            }

            if (x < Consts.BesselJY.ApproxThreshold) {
                return BesselYNearZero(nu, x).Convert<N>();
            }
            else {
                return BesselYLimit(nu, x).Convert<N>();
            }
        }

        private static MultiPrecision<Plus1<N>> BesselJNearZero(MultiPrecision<N> nu, MultiPrecision<N> z) {
            Consts.BesselNearZeroCoef table = Consts.Bessel.NearZeroCoef(nu);

            MultiPrecision<Double<N>> z_ex = z.Convert<Double<N>>();
            MultiPrecision<Double<N>> u = 1;
            MultiPrecision<Double<N>> w = z_ex * z_ex;

            MultiPrecision<Double<N>> x = 0;

            Sign sign = Sign.Plus;

            for (int k = 0; k < int.MaxValue; k++) {
                MultiPrecision<Double<N>> c = u * table.Value(k);

                if (sign == Sign.Plus) {
                    x += c;
                    sign = Sign.Minus;
                }
                else {
                    x -= c;
                    sign = Sign.Plus;
                }
                u *= w;

                if (c.IsZero || x.Exponent - c.Exponent > MultiPrecision<Plus1<N>>.Bits) {
                    break;
                }
            }

            MultiPrecision<Double<N>> p;
            if (nu == 0) {
                p = 1;
            }
            else if (nu == 1) {
                p = z_ex / 2;
            }
            else if (nu == 2) {
                p = z_ex * z_ex / 4;
            }
            else {
                p = MultiPrecision<Double<N>>.Pow(z_ex / 2, nu.Convert<Double<N>>());
            }

            MultiPrecision<Double<N>> y = x * p;

            return y.Convert<Plus1<N>>();
        }

        private static MultiPrecision<Plus1<N>> BesselJLimit(MultiPrecision<N> nu, MultiPrecision<N> z) {
            Consts.BesselLimitCoef table = Consts.Bessel.LimitCoef(nu);

            MultiPrecision<Plus4<N>> z_ex = z.Convert<Plus4<N>>();
            MultiPrecision<Plus4<N>> v = 1 / z_ex;
            MultiPrecision<Plus4<N>> w = v * v;

            MultiPrecision<Plus4<N>> x = 0, y = 0, p = 1, q = v;

            Sign sign = Sign.Plus;

            for (int k = 0; k <= Consts.BesselJY.LimitApproxTerms; k++) {
                MultiPrecision<Plus4<N>> c = p * table.Value(k * 2);
                MultiPrecision<Plus4<N>> s = q * table.Value(k * 2 + 1);

                if (sign == Sign.Plus) {
                    x += c;
                    y += s;
                    sign = Sign.Minus;
                }
                else {
                    x -= c;
                    y -= s;
                    sign = Sign.Plus;
                }

                p *= w;
                q *= w;

                if (!c.IsZero && x.Exponent - c.Exponent <= MultiPrecision<Plus1<N>>.Bits) {
                    continue;
                }
                if (!s.IsZero && y.Exponent - s.Exponent <= MultiPrecision<Plus1<N>>.Bits) {
                    continue;
                }

                break;
            }

            MultiPrecision<Plus4<N>> omega = z_ex - (2 * nu.Convert<Plus4<N>>() + 1) * MultiPrecision<Plus4<N>>.PI / 4;
            MultiPrecision<Plus4<N>> m = x * MultiPrecision<Plus4<N>>.Cos(omega) - y * MultiPrecision<Plus4<N>>.Sin(omega);
            MultiPrecision<Plus4<N>> t = m * MultiPrecision<Plus4<N>>.Sqrt(2 / (MultiPrecision<Plus4<N>>.PI * z_ex));

            return t.Convert<Plus1<N>>();
        }

        private static MultiPrecision<Plus1<N>> BesselYNearZero(MultiPrecision<N> nu, MultiPrecision<N> z) {
            if (nu != MultiPrecision<N>.Truncate(nu)) {
                MultiPrecision<Plus4<N>> nu_ex = nu.Convert<Plus4<N>>(), z_ex = z.Convert<Plus4<N>>();

                MultiPrecision<Plus4<N>> bessel_j_pos = MultiPrecisionSandbox<Plus4<N>>.BesselJ(nu_ex, z_ex);
                MultiPrecision<Plus4<N>> bessel_j_neg = MultiPrecisionSandbox<Plus4<N>>.BesselJ(-nu_ex, z_ex);

                (MultiPrecision<Plus4<N>> sin, MultiPrecision<Plus4<N>> cos) = Consts.Bessel.SinCos(nu);

                MultiPrecision<Plus4<N>> y = (bessel_j_pos * cos - bessel_j_neg) / sin;

                return y.Convert<Plus1<N>>();
            }
            else {
                int n = (int)MultiPrecision<N>.Truncate(nu);

                Consts.BesselIntegerFiniteTermCoef finite_table = Consts.Bessel.IntegerFiniteTermCoef(n);
                Consts.BesselIntegerConvergenceTermCoef convergence_table = Consts.Bessel.IntegerConvergenceTermCoef(n);

                MultiPrecision<Double<N>> z_ex = z.Convert<Double<N>>();
                MultiPrecision<Double<N>> u = 1;
                MultiPrecision<Double<N>> w = z_ex * z_ex;

                MultiPrecision<Double<N>> r
                    = 2 * BesselJNearZero(n, z).Convert<Double<N>>() * MultiPrecision<Double<N>>.Log(z.Convert<Double<N>>() / 2);

                long r_exponent = MultiPrecision<Pow2.N4>.Sqrt(2 / (MultiPrecision<Pow2.N4>.PI * (z.Convert<Pow2.N4>() + MultiPrecision<Pow2.N4>.Point5))).Exponent;
                
                MultiPrecision<Double<N>> m = MultiPrecision<Double<N>>.Pow(z_ex / 2, n);

                MultiPrecision<Double<N>> x = 0, y = 0;

                Sign sign = Sign.Plus;

                for (int k = 0; k < int.MaxValue; k++, u *= w) {
                    MultiPrecision<Double<N>> c = u * convergence_table.Value(k);

                    if (sign == Sign.Plus) {
                        x += c;
                        sign = Sign.Minus;
                    }
                    else {
                        x -= c;
                        sign = Sign.Plus;
                    }

                    if (k < n) {
                        y += u * finite_table.Value(k);                        
                        continue;
                    }

                    if (c.IsZero || Math.Min(x.Exponent - c.Exponent, r_exponent - c.Exponent - m.Exponent) > MultiPrecision<Plus1<N>>.Bits) {
                        break;
                    }
                }

                MultiPrecision<Double<N>> d = (r - y / m - x * m) / MultiPrecision<Double<N>>.PI;

                return d.Convert<Plus1<N>>();
            }
        }

        private static MultiPrecision<Plus1<N>> BesselYLimit(MultiPrecision<N> nu, MultiPrecision<N> z) {
            Consts.BesselLimitCoef table = Consts.Bessel.LimitCoef(nu);

            MultiPrecision<Plus4<N>> z_ex = z.Convert<Plus4<N>>();
            MultiPrecision<Plus4<N>> v = 1 / z_ex;
            MultiPrecision<Plus4<N>> w = v * v;

            MultiPrecision<Plus4<N>> x = 0, y = 0, p = 1, q = v;

            Sign sign = Sign.Plus;

            for (int k = 0; k <= Consts.BesselJY.LimitApproxTerms; k++) {
                MultiPrecision<Plus4<N>> c = p * table.Value(k * 2);
                MultiPrecision<Plus4<N>> s = q * table.Value(k * 2 + 1);

                if (sign == Sign.Plus) {
                    x += c;
                    y += s;
                    sign = Sign.Minus;
                }
                else {
                    x -= c;
                    y -= s;
                    sign = Sign.Plus;
                }

                p *= w;
                q *= w;

                if (!c.IsZero && x.Exponent - c.Exponent <= MultiPrecision<Plus1<N>>.Bits) {
                    continue;
                }
                if (!s.IsZero && y.Exponent - s.Exponent <= MultiPrecision<Plus1<N>>.Bits) {
                    continue;
                }

                break;
            }

            MultiPrecision<Plus4<N>> omega = z_ex - (2 * nu.Convert<Plus4<N>>() + 1) * MultiPrecision<Plus4<N>>.PI / 4;
            MultiPrecision<Plus4<N>> m = x * MultiPrecision<Plus4<N>>.Sin(omega) + y * MultiPrecision<Plus4<N>>.Cos(omega);
            MultiPrecision<Plus4<N>> t = m * MultiPrecision<Plus4<N>>.Sqrt(2 / (MultiPrecision<Plus4<N>>.PI * z_ex));

            return t.Convert<Plus1<N>>();
        }

        private static partial class Consts {
            public static class BesselJY {
                public static MultiPrecision<N> ApproxThreshold { private set; get; }

                public static int LimitApproxTerms { private set; get; }

                static BesselJY() {
                    if (MultiPrecision<N>.Length > 128) {
                        throw new ArgumentOutOfRangeException(nameof(MultiPrecision<N>.Length));
                    }

                    ApproxThreshold = Math.Ceiling(50 + 11.0965 * MultiPrecision<N>.Length);

                    LimitApproxTerms = (int)Math.Ceiling(41 + 11.0190 * MultiPrecision<N>.Length);

#if DEBUG
                    Trace.WriteLine($"BesselJY<{MultiPrecision<N>.Length}> initialized.");
#endif
                }
            }

            public static class Bessel {
                private readonly static Dictionary<MultiPrecision<N>, BesselNearZeroCoef> nearzero_table = new();
                private readonly static Dictionary<MultiPrecision<N>, BesselLimitCoef> limit_table = new();

                private readonly static Dictionary<int, BesselIntegerFiniteTermCoef> integer_finite_table = new();
                private readonly static Dictionary<int, BesselIntegerConvergenceTermCoef> interger_convergence_table = new();

                private readonly static Dictionary<MultiPrecision<N>, (MultiPrecision<Plus4<N>> sin, MultiPrecision<Plus4<N>> cos)> sincos_table = new();

                public static BesselNearZeroCoef NearZeroCoef(MultiPrecision<N> nu) {
                    BesselNearZeroCoef table;
                    if (nearzero_table.ContainsKey(nu)) {
                        table = nearzero_table[nu];
                    }
                    else {
                        table = new BesselNearZeroCoef(nu);
                        nearzero_table.Add(nu, table);
                    }

                    return table;
                }

                public static BesselLimitCoef LimitCoef(MultiPrecision<N> nu) {
                    BesselLimitCoef table;
                    if (limit_table.ContainsKey(nu)) {
                        table = limit_table[nu];
                    }
                    else {
                        table = new BesselLimitCoef(nu);
                        limit_table.Add(nu, table);
                    }

                    return table;
                }

                public static BesselIntegerFiniteTermCoef IntegerFiniteTermCoef(int n) {
                    BesselIntegerFiniteTermCoef table;
                    if (integer_finite_table.ContainsKey(n)) {
                        table = integer_finite_table[n];
                    }
                    else {
                        table = new BesselIntegerFiniteTermCoef(n);
                        integer_finite_table.Add(n, table);
                    }

                    return table;
                }

                public static BesselIntegerConvergenceTermCoef IntegerConvergenceTermCoef(int n) {
                    BesselIntegerConvergenceTermCoef table;
                    if (interger_convergence_table.ContainsKey(n)) {
                        table = interger_convergence_table[n];
                    }
                    else {
                        table = new BesselIntegerConvergenceTermCoef(n);
                        interger_convergence_table.Add(n, table);
                    }

                    return table;
                }

                public static (MultiPrecision<Plus4<N>> sin, MultiPrecision<Plus4<N>> cos) SinCos(MultiPrecision<N> nu) {
                    if (!sincos_table.ContainsKey(nu)) {
                        sincos_table.Add(nu, (MultiPrecision<Plus4<N>>.SinPI(nu.Convert<Plus4<N>>()), MultiPrecision<Plus4<N>>.CosPI(nu.Convert<Plus4<N>>())));
                    }

                    return sincos_table[nu];
                }
            }

            public class BesselNearZeroCoef {
                private readonly MultiPrecision<Double<N>> nu;
                private readonly List<MultiPrecision<Double<N>>> a_table = new();
                private readonly List<MultiPrecision<Double<N>>> c_table = new();

                public BesselNearZeroCoef(MultiPrecision<N> nu) {
                    this.nu = nu.Convert<Double<N>>();

                    MultiPrecision<Double<N>> a0;

                    if (nu >= 0 && nu == MultiPrecision<N>.Truncate(nu)) {
                        long n = (long)nu;

                        a0 = 1;
                        for (int k = 2; k <= n; k++) {
                            a0 *= k;
                        }
                    }
                    else {
                        a0 = MultiPrecision<Double<N>>.Length < 256 ?
                             MultiPrecision<Double<N>>.Gamma(nu.Convert<Double<N>>() + 1) :
                             MultiPrecision<Pow2.N256>.Gamma(nu.Convert<Pow2.N256>() + 1).Convert<Double<N>>();
                    }

                    this.a_table.Add(a0);
                    this.c_table.Add(1 / a0);
                }

                public MultiPrecision<Double<N>> Value(int n) {
                    if (n < 0) {
                        throw new ArgumentOutOfRangeException(nameof(n));
                    }

                    if (n < c_table.Count) {
                        return c_table[n];
                    }

                    for (int k = c_table.Count; k <= n; k++) {
                        MultiPrecision<Double<N>> a = a_table.Last() * (checked(4 * k) * (nu + k));

                        a_table.Add(a);
                        c_table.Add(1 / a);
                    }

                    return c_table[n];
                }
            }

            public class BesselLimitCoef {
                private readonly MultiPrecision<Plus4<N>> squa_nu4;
                private readonly List<MultiPrecision<Plus4<N>>> a_table = new();

                public BesselLimitCoef(MultiPrecision<N> nu) {
                    this.squa_nu4 = 4 * MultiPrecision<Plus4<N>>.Square(nu.Convert<Plus4<N>>());

                    MultiPrecision<Plus4<N>> a1 = (squa_nu4 - 1) / 8;

                    this.a_table.Add(1);
                    this.a_table.Add(a1);
                }

                public MultiPrecision<Plus4<N>> Value(int n) {
                    if (n < 0) {
                        throw new ArgumentOutOfRangeException(nameof(n));
                    }

                    if (n < a_table.Count) {
                        return a_table[n];
                    }

                    for (int k = a_table.Count; k <= n; k++) {
                        MultiPrecision<Plus4<N>> a =
                            a_table.Last() * MultiPrecision<Plus4<N>>.Div(checked(squa_nu4 - (2 * k - 1) * (2 * k - 1)), checked(k * 8));

                        a_table.Add(a);
                    }

                    return a_table[n];
                }
            }

            public class BesselIntegerFiniteTermCoef {
                private readonly List<MultiPrecision<Double<N>>> a_table = new();

                public BesselIntegerFiniteTermCoef(int n) {
                    if (n < 0) {
                        throw new ArgumentOutOfRangeException(nameof(n));
                    }

                    MultiPrecision<Double<N>> a = 1;
                    for (int i = 2; i < n; i++) {
                        a *= i;
                    }

                    this.a_table.Add(a);

                    for (int k = 1; k < n; k++) {
                        a /= checked(4 * k * (n - k));
                        this.a_table.Add(a);
                    }
                }

                public MultiPrecision<Double<N>> Value(int k) {
                    if (k < 0 || k >= a_table.Count) {
                        throw new ArgumentOutOfRangeException(nameof(k));
                    }

                    return a_table[k];
                }
            }

            public class BesselIntegerConvergenceTermCoef {
                private static readonly MultiPrecision<Double<N>> b;

                private readonly int n;
                private readonly List<MultiPrecision<Double<N>>> a_table = new();
                private MultiPrecision<Double<N>> r;

                static BesselIntegerConvergenceTermCoef() {
                    b = 2 * MultiPrecision<Double<N>>.EulerGamma;
                }

                public BesselIntegerConvergenceTermCoef(int n) {
                    if (n < 0) {
                        throw new ArgumentOutOfRangeException(nameof(n));
                    }

                    MultiPrecision<Double<N>> r0 = 1;
                    for (int i = 2; i <= n; i++) {
                        r0 *= i;
                    }

                    MultiPrecision<Double<N>> a0 = (MultiPrecision<Double<N>>.HarmonicNumber(n) - b) / r0;

                    this.n = n;
                    this.r = r0;
                    this.a_table.Add(a0);
                }

                public MultiPrecision<Double<N>> Value(int k) {
                    if (k < 0) {
                        throw new ArgumentOutOfRangeException(nameof(k));
                    }

                    if (k < a_table.Count) {
                        return a_table[k];
                    }

                    for (int i = a_table.Count; i <= k; i++) {
                        r *= checked(4 * i * (n + i));
                        MultiPrecision<Double<N>> a =
                            (MultiPrecision<Double<N>>.HarmonicNumber(i) + MultiPrecision<Double<N>>.HarmonicNumber(n + i) - b) / r;

                        a_table.Add(a);
                    }

                    return a_table[k];
                }
            }
        }
    }
}
