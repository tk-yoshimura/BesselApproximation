﻿using MultiPrecision;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BesselApproximation {

    public static class MultiPrecisionSandbox<N> where N : struct, IConstant {

        public static MultiPrecision<N> BesselJ(MultiPrecision<N> nu, MultiPrecision<N> x) {
            if (MultiPrecision<N>.Abs(nu) > 1024) {
                throw new ArgumentOutOfRangeException(nameof(nu));
            }
            if (nu.IsNaN || x.IsNaN) {
                return MultiPrecision<N>.NaN;
            }

            if (x.Sign == Sign.Minus) {
                if (nu != MultiPrecision<N>.Floor(nu)) {
                    return MultiPrecision<N>.NaN;
                }

                long n = (long)nu;
                return ((n & 1L) == 0) ? BesselJ(nu, MultiPrecision<N>.Abs(x)) : -BesselJ(nu, MultiPrecision<N>.Abs(x));
            }
            if (x.IsZero) {
                return nu.IsZero ? 1 : 0;
            }
            if (nu.Sign == Sign.Minus && nu == MultiPrecision<N>.Floor(nu)) {
                long n = (long)nu;
                return ((n & 1L) == 0) ? BesselJ(MultiPrecision<N>.Abs(nu), x) : -BesselJ(MultiPrecision<N>.Abs(nu), x);
            }

            if (nu - MultiPrecision<N>.Point5 == MultiPrecision<N>.Floor(nu)) {
                long n = (long)MultiPrecision<N>.Floor(nu);

                if (n >= -2 || n < 2) {
                    MultiPrecision<Plus1<N>> x_ex = x.Convert<Plus1<N>>();
                    MultiPrecision<Plus1<N>> envelope = MultiPrecision<Plus1<N>>.Sqrt(2 / (MultiPrecision<Plus1<N>>.PI * x_ex));

                    if (n == -2) {
                        return (envelope * -(MultiPrecision<Plus1<N>>.Cos(x_ex) / x_ex - MultiPrecision<Plus1<N>>.Sin(x_ex))).Convert<N>();
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

            if (x < Consts.BesselJ.ApproxThreshold) {
                return BesselJNearZero(nu, x).Convert<N>();
            }
            else {
                return BesselJLimit(nu, x).Convert<N>();
            }
        }

        private static MultiPrecision<Plus1<N>> BesselJNearZero(MultiPrecision<N> nu, MultiPrecision<N> z) {
            BesselNearZeroCoef<Double<N>> table = Consts.Bessel.NearZeroCoef(nu);

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
            BesselLimitCoef<Plus4<N>> table = Consts.Bessel.LimitCoef(nu);

            MultiPrecision<Plus4<N>> z_ex = z.Convert<Plus4<N>>();
            MultiPrecision<Plus4<N>> v = 1 / z_ex;
            MultiPrecision<Plus4<N>> w = v * v;
                               
            MultiPrecision<Plus4<N>> x = 0, y = 0, p = 1, q = v;

            Sign sign = Sign.Plus;

            for(int k = 0; k <= Consts.BesselJ.LimitApproxTerms; k++) {
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

        private static partial class Consts {
            public static class BesselJ {
                public static MultiPrecision<N> ApproxThreshold { private set; get; }

                public static int LimitApproxTerms { private set; get; }

                static BesselJ() {
                    if (MultiPrecision<N>.Length > 128) {
                        throw new ArgumentOutOfRangeException(nameof(MultiPrecision<N>.Length));
                    }

                    ApproxThreshold = Math.Ceiling(50 + 11.0965 * MultiPrecision<N>.Length);

                    LimitApproxTerms = (int)Math.Ceiling(41 + 11.0190 * MultiPrecision<N>.Length);

#if DEBUG
                    Trace.WriteLine($"BesselJ<{MultiPrecision<N>.Length}> initialized.");
#endif
                }
            }

            public static class Bessel {
                private readonly static Dictionary<MultiPrecision<N>, BesselNearZeroCoef<Double<N>>> nearzero_table = new();
                private readonly static Dictionary<MultiPrecision<N>, BesselLimitCoef<Plus4<N>>> limit_table = new();

                public static BesselNearZeroCoef<Double<N>> NearZeroCoef(MultiPrecision<N> nu) {
                    BesselNearZeroCoef<Double<N>> table;
                    if (nearzero_table.ContainsKey(nu)) {
                        table = nearzero_table[nu];
                    }
                    else {
                        table = new BesselNearZeroCoef<Double<N>>(nu.Convert<Double<N>>());
                        nearzero_table.Add(nu, table);
                    }

                    return table;
                }

                public static BesselLimitCoef<Plus4<N>> LimitCoef(MultiPrecision<N> nu) {
                    BesselLimitCoef<Plus4<N>> table;
                    if (limit_table.ContainsKey(nu)) {
                        table = limit_table[nu];
                    }
                    else {
                        table = new BesselLimitCoef<Plus4<N>>(nu.Convert<Plus4<N>>());
                        limit_table.Add(nu, table);
                    }

                    return table;
                }
            }

            public class BesselNearZeroCoef {
                private readonly MultiPrecision<N> nu;
                private readonly List<MultiPrecision<N>> a_table = new();
                private readonly List<MultiPrecision<N>> c_table = new();

                public BesselNearZeroCoef(MultiPrecision<N> nu) {
                    this.nu = nu;

                    MultiPrecision<N> a0;

                    if (nu >= 0 && nu == MultiPrecision<N>.Floor(nu)) {
                        long n = (long)nu;

                        a0 = 1;
                        for (int k = 2; k <= n; k++) {
                            a0 *= k;
                        }
                    }
                    else {
                        a0 = MultiPrecision<N>.Gamma(nu + 1);
                    }

                    this.a_table.Add(a0);
                    this.c_table.Add(1 / a0);
                }

                public MultiPrecision<N> Value(int n) {
                    if (n < 0) {
                        throw new ArgumentOutOfRangeException(nameof(n));
                    }

                    if (n < c_table.Count) {
                        return c_table[n];
                    }

                    for (int k = c_table.Count; k <= n; k++) {
                        MultiPrecision<N> a = a_table.Last() * (checked(4 * k) * (nu + k));

                        a_table.Add(a);
                        c_table.Add(1 / a);
                    }

                    return c_table[n];
                }
            }

            public class BesselLimitCoef {
                private readonly MultiPrecision<N> squa_nu4;
                private readonly List<MultiPrecision<N>> a_table = new();

                public BesselLimitCoef(MultiPrecision<N> nu) {
                    this.squa_nu4 = 4 * nu * nu;

                    MultiPrecision<N> a1 = (squa_nu4 - 1) / 8;

                    this.a_table.Add(1);
                    this.a_table.Add(a1);
                }

                public MultiPrecision<N> Value(int n) {
                    if (n < 0) {
                        throw new ArgumentOutOfRangeException(nameof(n));
                    }

                    if (n < a_table.Count) {
                        return a_table[n];
                    }

                    for (int k = a_table.Count; k <= n; k++) {
                        MultiPrecision<N> a =
                            a_table.Last() * MultiPrecision<N>.Div(checked(squa_nu4 - (2 * k - 1) * (2 * k - 1)), checked(k * 8));

                        a_table.Add(a);
                    }

                    return a_table[n];
                }
            }
        }
    }
}
