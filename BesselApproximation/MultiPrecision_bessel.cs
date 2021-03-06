using MultiPrecision;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BesselApproximation {

    public static class MultiPrecisionSandbox<N> where N : struct, IConstant {

        public static MultiPrecision<N> BesselJ(MultiPrecision<N> nu, MultiPrecision<N> x) {
            if (MultiPrecision<N>.Abs(nu) > 64) {
                throw new ArgumentOutOfRangeException(
                    nameof(nu),
                    "In the calculation of the Bessel function, nu with an absolute value greater than 64 is not supported."
                );
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

            if (!x.IsFinite) {
                return 0;
            }
            if ((x / 2).IsZero && nu.IsZero) {
                return 1;
            }
            if (nu.Sign == Sign.Minus && nu == MultiPrecision<N>.Truncate(nu)) {
                long n = (long)nu;
                return ((n & 1L) == 0) ? BesselJ(MultiPrecision<N>.Abs(nu), x) : -BesselJ(MultiPrecision<N>.Abs(nu), x);
            }

            if (nu - MultiPrecision<N>.Point5 == MultiPrecision<N>.Floor(nu)) {
                long n = (long)MultiPrecision<N>.Floor(nu);

                if (n >= -2 && n < 2) {
                    MultiPrecision<Plus1<N>> x_ex = x.Convert<Plus1<N>>();
                    MultiPrecision<Plus1<N>> envelope = MultiPrecision<Plus1<N>>.Sqrt(2 / (MultiPrecision<Plus1<N>>.PI * x_ex));

                    if (n == -2) {
                        return -(envelope * (MultiPrecision<Plus1<N>>.Cos(x_ex) / x_ex + MultiPrecision<Plus1<N>>.Sin(x_ex))).Convert<N>();
                    }
                    if (n == -1) {
                        return (envelope * MultiPrecision<Plus1<N>>.Cos(x_ex)).Convert<N>();
                    }
                    if (n == 0) {
                        MultiPrecision<N> y = (envelope * MultiPrecision<Plus1<N>>.Sin(x_ex)).Convert<N>();

                        return y.IsNormal ? y : 0;
                    }
                    if (n == 1) {
                        MultiPrecision<N> y = (envelope * (MultiPrecision<Plus1<N>>.Sin(x_ex) / x_ex - MultiPrecision<Plus1<N>>.Cos(x_ex))).Convert<N>();

                        return y.IsNormal ? y : 0;
                    }
                }
            }

            if (MultiPrecision<N>.Length <= 4) {
                return MultiPrecisionSandbox<Plus1<N>>.BesselJ(nu.Convert<Plus1<N>>(), x.Convert<Plus1<N>>()).Convert<N>();
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
                throw new ArgumentOutOfRangeException(
                    nameof(nu),
                    "In the calculation of the Bessel function, nu with an absolute value greater than 64 is not supported."
                );
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

            if (!x.IsFinite) {
                return 0;
            }
            if (x.IsZero) {
                return MultiPrecision<N>.NegativeInfinity;
            }
            if (nu.Sign == Sign.Minus && nu == MultiPrecision<N>.Truncate(nu)) {
                long n = (long)nu;
                return ((n & 1L) == 0) ? BesselY(MultiPrecision<N>.Abs(nu), x) : -BesselY(MultiPrecision<N>.Abs(nu), x);
            }

            if (nu - MultiPrecision<N>.Point5 == MultiPrecision<N>.Floor(nu)) {
                long n = (long)MultiPrecision<N>.Floor(nu);

                if (n >= -2 && n < 2) {
                    MultiPrecision<Plus1<N>> x_ex = x.Convert<Plus1<N>>();
                    MultiPrecision<Plus1<N>> envelope = MultiPrecision<Plus1<N>>.Sqrt(2 / (MultiPrecision<Plus1<N>>.PI * x_ex));

                    if (n == -2) {
                        MultiPrecision<N> y = -(envelope * (MultiPrecision<Plus1<N>>.Sin(x_ex) / x_ex - MultiPrecision<Plus1<N>>.Cos(x_ex))).Convert<N>();

                        return y.IsNormal ? y : 0;
                    }
                    if (n == -1) {
                        MultiPrecision<N> y = (envelope * MultiPrecision<Plus1<N>>.Sin(x_ex)).Convert<N>();

                        return y.IsNormal ? y : 0;
                    }
                    if (n == 0) {
                        return -(envelope * MultiPrecision<Plus1<N>>.Cos(x_ex)).Convert<N>();
                    }
                    if (n == 1) {
                        return -(envelope * (MultiPrecision<Plus1<N>>.Cos(x_ex) / x_ex + MultiPrecision<Plus1<N>>.Sin(x_ex))).Convert<N>();
                    }
                }
            }

            if (MultiPrecision<N>.Length <= 4) {
                return MultiPrecisionSandbox<Plus1<N>>.BesselY(nu.Convert<Plus1<N>>(), x.Convert<Plus1<N>>()).Convert<N>();
            }

            if (x < Consts.BesselJY.ApproxThreshold) {
                return BesselYNearZero(nu, x);
            }
            else {
                return BesselYLimit(nu, x).Convert<N>();
            }
        }

        public static MultiPrecision<N> BesselI(MultiPrecision<N> nu, MultiPrecision<N> x) {
            if (MultiPrecision<N>.Abs(nu) > 64) {
                throw new ArgumentOutOfRangeException(
                    nameof(nu),
                    "In the calculation of the Bessel function, nu with an absolute value greater than 64 is not supported."
                );
            }
            if (nu.IsNaN || x.IsNaN || x.Sign == Sign.Minus) {
                return MultiPrecision<N>.NaN;
            }

            if (nu.Sign == Sign.Minus && nu == MultiPrecision<N>.Truncate(nu)) {
                return BesselI(MultiPrecision<N>.Abs(nu), x);
            }

            if (nu - MultiPrecision<N>.Point5 == MultiPrecision<N>.Floor(nu)) {
                long n = (long)MultiPrecision<N>.Floor(nu);

                if (n >= -2 && n < 2) {
                    MultiPrecision<Plus1<N>> x_ex = x.Convert<Plus1<N>>();
                    MultiPrecision<Plus1<N>> r = MultiPrecision<Plus1<N>>.Sqrt2 / MultiPrecision<Plus1<N>>.Sqrt(MultiPrecision<Plus1<N>>.PI * x_ex);

                    if (n == -2) {
                        return -(r * (MultiPrecision<Plus1<N>>.Cosh(x_ex) / x_ex - MultiPrecision<Plus1<N>>.Sinh(x_ex))).Convert<N>();
                    }
                    if (n == -1) {
                        return (r * MultiPrecision<Plus1<N>>.Cosh(x_ex)).Convert<N>();
                    }
                    if (n == 0) {
                        MultiPrecision<N> y = (r * MultiPrecision<Plus1<N>>.Sinh(x_ex)).Convert<N>();

                        return y.IsNormal ? y : 0;
                    }
                    if (n == 1) {
                        MultiPrecision<N> y = -(r * (MultiPrecision<Plus1<N>>.Sinh(x_ex) / x_ex - MultiPrecision<Plus1<N>>.Cosh(x_ex))).Convert<N>();

                        return y.IsNormal ? y : 0;
                    }
                }
            }

            if (x < Consts.BesselIK.ApproxThreshold) {
                return BesselINearZero(nu, x).Convert<N>();
            }
            else {
                return BesselILimit(nu, x).Convert<N>();
            }
        }

        public static MultiPrecision<N> BesselK(MultiPrecision<N> nu, MultiPrecision<N> x) {
            if (MultiPrecision<N>.Abs(nu) > 64) {
                throw new ArgumentOutOfRangeException(
                    nameof(nu),
                    "In the calculation of the Bessel function, nu with an absolute value greater than 64 is not supported."
                );
            }
            if (nu.IsNaN || x.IsNaN || x.Sign == Sign.Minus) {
                return MultiPrecision<N>.NaN;
            }

            if (x.IsZero) {
                return MultiPrecision<N>.PositiveInfinity;
            }
            if (nu.Sign == Sign.Minus) {
                return BesselK(MultiPrecision<N>.Abs(nu), x);
            }

            if (nu - MultiPrecision<N>.Point5 == MultiPrecision<N>.Floor(nu)) {
                long n = (long)MultiPrecision<N>.Floor(nu);

                if (n >= 0 && n < 2) {
                    MultiPrecision<Plus1<N>> x_ex = x.Convert<Plus1<N>>();
                    MultiPrecision<Plus1<N>> r = MultiPrecision<Plus1<N>>.Exp(-x_ex) * MultiPrecision<Plus1<N>>.Sqrt(MultiPrecision<Plus1<N>>.PI / (2 * x_ex));

                    if (n == 0) {
                        return r.Convert<N>();
                    }
                    if (n == 1) {
                        return (r * (1 + 1 / x_ex)).Convert<N>();
                    }
                }
            }

            if (x < Consts.BesselIK.ApproxThreshold) {
                return BesselKNearZero(nu, x);
            }
            else {
                return BesselKLimit(nu, x).Convert<N>();
            }
        }

        private static MultiPrecision<Plus1<N>> BesselJNearZero(MultiPrecision<N> nu, MultiPrecision<N> z) {
            Consts.BesselNearZeroCoef table = Consts.Bessel.NearZeroCoef(nu);

            MultiPrecision<Double<N>> z_ex = z.Convert<Double<N>>();
            MultiPrecision<Double<N>> u = 1;
            MultiPrecision<Double<N>> w = z_ex * z_ex;

            MultiPrecision<Double<N>> x = 0;

            Sign sign = Sign.Plus;

            for (int k = 0; k < int.MaxValue; k++, u *= w) {
                MultiPrecision<Double<N>> c = u * table.Value(k);

                if (sign == Sign.Plus) {
                    x += c;
                    sign = Sign.Minus;
                }
                else {
                    x -= c;
                    sign = Sign.Plus;
                }

                if (c.IsZero || x.Exponent - c.Exponent > MultiPrecision<Plus1<N>>.Bits) {
                    break;
                }

                if (k >= MultiPrecision<N>.Bits && Math.Max(x.Exponent, c.Exponent) < -MultiPrecision<N>.Bits * 8) {
                    return 0;
                }
            }

            MultiPrecision<Double<N>> p;
            if (nu == MultiPrecision<N>.Truncate(nu)) {
                int n = (int)nu;

                p = MultiPrecision<Double<N>>.Pow(z_ex / 2, n);
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

            for (int k = 0; k <= Consts.BesselJY.LimitApproxTerms; k++, p *= w, q *= w) {
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

        private static MultiPrecision<N> BesselYNearZero(MultiPrecision<N> nu, MultiPrecision<N> z) {
            int n = (int)MultiPrecision<N>.Round(nu);

            if (nu != n) {
                MultiPrecision<N> dnu = nu - n;

                if (dnu.Exponent >= -16) {
                    return BesselYNonIntegerNu(nu, z);
                }
                if (dnu.Exponent >= -32) {
                    return MultiPrecisionSandbox<Plus1<N>>.BesselYNonIntegerNu(nu.Convert<Plus1<N>>(), z.Convert<Plus1<N>>()).Convert<N>();
                }
                if (dnu.Exponent >= -48) {
                    return MultiPrecisionSandbox<Plus2<N>>.BesselYNonIntegerNu(nu.Convert<Plus2<N>>(), z.Convert<Plus2<N>>()).Convert<N>();
                }
                if (dnu.Exponent >= -80) {
                    return MultiPrecisionSandbox<Plus4<N>>.BesselYNonIntegerNu(nu.Convert<Plus4<N>>(), z.Convert<Plus4<N>>()).Convert<N>();
                }
                if (dnu.Exponent >= -144) {
                    return MultiPrecisionSandbox<Plus8<N>>.BesselYNonIntegerNu(nu.Convert<Plus8<N>>(), z.Convert<Plus8<N>>()).Convert<N>();
                }
                if (dnu.Exponent >= -272) {
                    return MultiPrecisionSandbox<Plus16<N>>.BesselYNonIntegerNu(nu.Convert<Plus16<N>>(), z.Convert<Plus16<N>>()).Convert<N>();
                }

                throw new ArgumentException(
                    "The calculation of the BesselY function value is invalid because it loses digits" +
                    " when nu is extremely close to an integer. (|nu - round(nu)| < 1.32 x 10^-82 and nu != round(nu))",
                    nameof(nu));
            }

            return BesselYIntegerNuNearZero(n, z);
        }

        private static MultiPrecision<N> BesselYNonIntegerNu(MultiPrecision<N> nu, MultiPrecision<N> z) {
            MultiPrecision<Plus1<N>> bessel_j_pos = BesselJNearZero(nu, z);
            MultiPrecision<Plus1<N>> bessel_j_neg = BesselJNearZero(-nu, z);

            (MultiPrecision<Plus1<N>> sin, MultiPrecision<Plus1<N>> cos) = Consts.Bessel.SinCos(nu);

            MultiPrecision<Plus1<N>> y = (bessel_j_pos * cos - bessel_j_neg) / sin;

            return y.Convert<N>();
        }

        private static MultiPrecision<N> BesselYIntegerNuNearZero(int n, MultiPrecision<N> z) {
            Consts.BesselIntegerFiniteTermCoef finite_table = Consts.Bessel.IntegerFiniteTermCoef(n);
            Consts.BesselIntegerConvergenceTermCoef convergence_table = Consts.Bessel.IntegerConvergenceTermCoef(n);

            MultiPrecision<Double<N>> z_ex = z.Convert<Double<N>>();
            MultiPrecision<Double<N>> u = 1;
            MultiPrecision<Double<N>> w = z_ex * z_ex;

            MultiPrecision<Double<N>> r
                = 2 * BesselJNearZero(n, z).Convert<Double<N>>() * MultiPrecision<Double<N>>.Log(z_ex / 2);

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

            return d.Convert<N>();
        }

        private static MultiPrecision<Plus1<N>> BesselYLimit(MultiPrecision<N> nu, MultiPrecision<N> z) {
            Consts.BesselLimitCoef table = Consts.Bessel.LimitCoef(nu);

            MultiPrecision<Plus4<N>> z_ex = z.Convert<Plus4<N>>();
            MultiPrecision<Plus4<N>> v = 1 / z_ex;
            MultiPrecision<Plus4<N>> w = v * v;

            MultiPrecision<Plus4<N>> x = 0, y = 0, p = 1, q = v;

            Sign sign = Sign.Plus;

            for (int k = 0; k <= Consts.BesselJY.LimitApproxTerms; k++, p *= w, q *= w) {
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

        private static MultiPrecision<Plus1<N>> BesselINearZero(MultiPrecision<N> nu, MultiPrecision<N> z) {
            Consts.BesselNearZeroCoef table = Consts.Bessel.NearZeroCoef(nu);

            MultiPrecision<Double<N>> z_ex = z.Convert<Double<N>>();
            MultiPrecision<Double<N>> u = 1;
            MultiPrecision<Double<N>> w = z_ex * z_ex;

            MultiPrecision<Double<N>> x = 0;

            for (int k = 0; k < int.MaxValue; k++, u *= w) {
                MultiPrecision<Double<N>> c = u * table.Value(k);

                x += c;

                if (c.IsZero || x.Exponent - c.Exponent > MultiPrecision<Plus1<N>>.Bits) {
                    break;
                }
            }

            MultiPrecision<Double<N>> p;
            if (nu == MultiPrecision<N>.Truncate(nu)) {
                int n = (int)nu;

                p = MultiPrecision<Double<N>>.Pow(z_ex / 2, n);
            }
            else {
                p = MultiPrecision<Double<N>>.Pow(z_ex / 2, nu.Convert<Double<N>>());
            }

            MultiPrecision<Double<N>> y = x * p;

            return y.Convert<Plus1<N>>();
        }

        private static MultiPrecision<Plus1<N>> BesselILimit(MultiPrecision<N> nu, MultiPrecision<N> z) {
            Consts.BesselLimitCoef table = Consts.Bessel.LimitCoef(nu);

            MultiPrecision<Plus4<N>> z_ex = z.Convert<Plus4<N>>();
            MultiPrecision<Plus4<N>> v = 1 / z_ex;

            MultiPrecision<Plus4<N>> x = 0, p = 1;

            Sign sign = Sign.Plus;

            for (int k = 0; k <= Consts.BesselIK.LimitApproxTerms; k++, p *= v) {
                MultiPrecision<Plus4<N>> c = p * table.Value(k);

                if (sign == Sign.Plus) {
                    x += c;
                    sign = Sign.Minus;
                }
                else {
                    x -= c;
                    sign = Sign.Plus;
                }

                if (c.IsZero || x.Exponent - c.Exponent > MultiPrecision<Plus1<N>>.Bits) {
                    break;
                }
            }

            MultiPrecision<Plus4<N>> r =
                MultiPrecision<Plus4<N>>.Exp(z_ex) / MultiPrecision<Plus4<N>>.Sqrt(2 * MultiPrecision<Plus4<N>>.PI * z_ex);

            MultiPrecision<Plus4<N>> y = r * x;

            return y.Convert<Plus1<N>>();
        }

        private static MultiPrecision<N> BesselKNearZero(MultiPrecision<N> nu, MultiPrecision<N> z) {
            int n = (int)MultiPrecision<N>.Round(nu);

            if (nu != n) {
                MultiPrecision<N> dnu = nu - n;

                if (dnu.Exponent >= -16) {
                    return MultiPrecisionSandbox<Double<Plus2<N>>>.BesselKNonIntegerNu(nu.Convert<Double<Plus2<N>>>(), z.Convert<Double<Plus2<N>>>()).Convert<N>();
                }
                if (dnu.Exponent >= -32) {
                    return MultiPrecisionSandbox<Double<Plus4<N>>>.BesselKNonIntegerNu(nu.Convert<Double<Plus4<N>>>(), z.Convert<Double<Plus4<N>>>()).Convert<N>();
                }
                if (dnu.Exponent >= -48) {
                    return MultiPrecisionSandbox<Double<Plus8<N>>>.BesselKNonIntegerNu(nu.Convert<Double<Plus8<N>>>(), z.Convert<Double<Plus8<N>>>()).Convert<N>();
                }
                if (dnu.Exponent >= -80) {
                    return MultiPrecisionSandbox<Double<Plus16<N>>>.BesselKNonIntegerNu(nu.Convert<Double<Plus16<N>>>(), z.Convert<Double<Plus16<N>>>()).Convert<N>();
                }
                if (dnu.Exponent >= -144) {
                    return MultiPrecisionSandbox<Double<Plus32<N>>>.BesselKNonIntegerNu(nu.Convert<Double<Plus32<N>>>(), z.Convert<Double<Plus32<N>>>()).Convert<N>();
                }
                if (dnu.Exponent >= -272) {
                    return MultiPrecisionSandbox<Double<Plus64<N>>>.BesselKNonIntegerNu(nu.Convert<Double<Plus64<N>>>(), z.Convert<Double<Plus64<N>>>()).Convert<N>();
                }

                throw new ArgumentException(
                    "The calculation of the BesselK function value is invalid because it loses digits" +
                    " when nu is extremely close to an integer. (|nu - round(nu)| < 1.32 x 10^-82 and nu != round(nu))",
                    nameof(nu));
            }

            return MultiPrecisionSandbox<Plus4<N>>.BesselKIntegerNuNearZero(n, z.Convert<Plus4<N>>()).Convert<N>();
        }

        private static MultiPrecision<N> BesselKNonIntegerNu(MultiPrecision<N> nu, MultiPrecision<N> z) {
            MultiPrecision<Plus1<N>> bessel_i_pos = BesselINearZero(nu, z);
            MultiPrecision<Plus1<N>> bessel_i_neg = BesselINearZero(-nu, z);

            (MultiPrecision<Plus1<N>> sin, _) = Consts.Bessel.SinCos(nu);

            MultiPrecision<Plus1<N>> y = MultiPrecision<Plus1<N>>.PI * (bessel_i_neg - bessel_i_pos) / (2 * sin);

            return y.Convert<N>();
        }

        private static MultiPrecision<N> BesselKIntegerNuNearZero(int n, MultiPrecision<N> z) {
            Consts.BesselIntegerFiniteTermCoef finite_table = Consts.Bessel.IntegerFiniteTermCoef(n);
            Consts.BesselIntegerConvergenceTermCoef convergence_table = Consts.Bessel.IntegerConvergenceTermCoef(n);

            MultiPrecision<Double<N>> z_ex = z.Convert<Double<N>>();
            MultiPrecision<Double<N>> u = 1;
            MultiPrecision<Double<N>> w = z_ex * z_ex;

            MultiPrecision<Double<N>> r
                = MultiPrecisionSandbox<Double<N>>.BesselINearZero(n, z_ex).Convert<Double<N>>() * MultiPrecision<Double<N>>.Log(z_ex / 2);

            long r_exponent = r.Exponent;

            MultiPrecision<Double<N>> m = MultiPrecision<Double<N>>.Pow(z_ex / 2, n);

            MultiPrecision<Double<N>> x = 0, y = 0;

            Sign sign = Sign.Plus;

            for (int k = 0; k < int.MaxValue; k++, u *= w) {
                MultiPrecision<Double<N>> c = u * convergence_table.Value(k);

                y += c;

                if (k < n) {
                    if (sign == Sign.Plus) {
                        x += u * finite_table.Value(k);
                        sign = Sign.Minus;
                    }
                    else {
                        x -= u * finite_table.Value(k);
                        sign = Sign.Plus;
                    }
                    continue;
                }

                if (r.Exponent < -MultiPrecision<N>.Bits * 8 && (c.IsZero || y.Exponent - c.Exponent > MultiPrecision<Plus1<N>>.Bits)) {
                    break;
                }

                if (c.IsZero || Math.Min(y.Exponent - c.Exponent, r_exponent - c.Exponent - m.Exponent) > MultiPrecision<Double<N>>.Bits) {
                    break;
                }
            }

            MultiPrecision<Double<N>> d = (x / m + ((n & 1) == 0 ? 1 : -1) * (y * m - 2 * r)) / 2;

            return d.Convert<N>();
        }

        private static MultiPrecision<Plus1<N>> BesselKLimit(MultiPrecision<N> nu, MultiPrecision<N> z) {
            Consts.BesselLimitCoef table = Consts.Bessel.LimitCoef(nu);

            MultiPrecision<Plus4<N>> z_ex = z.Convert<Plus4<N>>();
            MultiPrecision<Plus4<N>> v = 1 / z_ex;

            MultiPrecision<Plus4<N>> x = 0, p = 1;

            for (int k = 0; k <= Consts.BesselIK.LimitApproxTerms; k++, p *= v) {
                MultiPrecision<Plus4<N>> c = p * table.Value(k);

                x += c;

                if (c.IsZero || x.Exponent - c.Exponent > MultiPrecision<Plus1<N>>.Bits) {
                    break;
                }
            }

            MultiPrecision<Plus4<N>> r =
                MultiPrecision<Plus4<N>>.Exp(-z_ex) * MultiPrecision<Plus4<N>>.Sqrt(MultiPrecision<Plus4<N>>.PI / (2 * z_ex));

            MultiPrecision<Plus4<N>> y = r * x;

            return y.Convert<Plus1<N>>();
        }

        private static partial class Consts {
            public static class BesselJY {
                public static MultiPrecision<N> ApproxThreshold { private set; get; }

                public static int LimitApproxTerms { private set; get; }

                static BesselJY() {
                    if (MultiPrecision<N>.Length > 64) {
                        throw new ArgumentOutOfRangeException(nameof(MultiPrecision<N>.Length));
                    }

                    ApproxThreshold = Math.Ceiling(50 + 11.0965 * MultiPrecision<N>.Length);

                    LimitApproxTerms = (int)Math.Ceiling(41 + 11.0190 * MultiPrecision<N>.Length);

#if DEBUG
                    Trace.WriteLine($"BesselJY<{MultiPrecision<N>.Length}> initialized.");
#endif
                }
            }

            public static class BesselIK {
                public static MultiPrecision<N> ApproxThreshold { private set; get; }

                public static int LimitApproxTerms { private set; get; }

                static BesselIK() {
                    if (MultiPrecision<N>.Length > 64) {
                        throw new ArgumentOutOfRangeException(nameof(MultiPrecision<N>.Length));
                    }

                    ApproxThreshold = Math.Ceiling(45 + 11.0965 * MultiPrecision<N>.Length);

                    LimitApproxTerms = (int)Math.Ceiling(80 + 21.2180 * MultiPrecision<N>.Length);

#if DEBUG
                    Trace.WriteLine($"BesselIK<{MultiPrecision<N>.Length}> initialized.");
#endif
                }
            }

            public static class Bessel {
                private readonly static Dictionary<MultiPrecision<N>, BesselNearZeroCoef> nearzero_table = new();
                private readonly static Dictionary<MultiPrecision<N>, BesselLimitCoef> limit_table = new();

                private readonly static Dictionary<int, BesselIntegerFiniteTermCoef> integer_finite_table = new();
                private readonly static Dictionary<int, BesselIntegerConvergenceTermCoef> interger_convergence_table = new();

                private readonly static Dictionary<MultiPrecision<N>, (MultiPrecision<Plus1<N>> sin, MultiPrecision<Plus1<N>> cos)> sincos_table = new();

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

                public static (MultiPrecision<Plus1<N>> sin, MultiPrecision<Plus1<N>> cos) SinCos(MultiPrecision<N> nu) {
                    if (!sincos_table.ContainsKey(nu)) {
                        sincos_table.Add(nu, (MultiPrecision<Plus1<N>>.SinPI(nu.Convert<Plus1<N>>()), MultiPrecision<Plus1<N>>.CosPI(nu.Convert<Plus1<N>>())));
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

                    for (long k = c_table.Count; k <= n; k++) {
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

                    for (long k = a_table.Count; k <= n; k++) {
                        MultiPrecision<Plus4<N>> a =
                            a_table.Last() * MultiPrecision<Plus4<N>>.Div(squa_nu4 - checked((2 * k - 1) * (2 * k - 1)), checked(k * 8));

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

                    for (long k = 1; k < n; k++) {
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

                    for (long i = a_table.Count; i <= k; i++) {
                        r *= checked(4 * i * (n + i));

                        MultiPrecision<Double<N>> a =
                            (MultiPrecision<Double<N>>.HarmonicNumber(checked((int)i))
                            + MultiPrecision<Double<N>>.HarmonicNumber(checked((int)(n + i))) - b) / r;

                        a_table.Add(a);
                    }

                    return a_table[k];
                }
            }
        }
    }
}
