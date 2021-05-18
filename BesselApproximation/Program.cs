using MultiPrecision;
using System;
using System.IO;

namespace BesselApproximation {
    class Program {
        static void Main(string[] args) {
            //BesselLimitConvergenceSummary<Pow2.N4>("../../../../results/bessel_limit_convergence_n4.txt");
            //BesselLimitConvergenceSummary<Pow2.N8>("../../../../results/bessel_limit_convergence_n8.txt");
            //BesselLimitConvergenceSummary<Pow2.N16>("../../../../results/bessel_limit_convergence_n16.txt");
            //BesselLimitConvergenceSummary<Pow2.N32>("../../../../results/bessel_limit_convergence_n32.txt");
            //BesselLimitConvergenceSummary<Pow2.N64>("../../../../results/bessel_limit_convergence_n64.txt");
            //BesselLimitConvergenceSummary<Pow2.N128>("../../../../results/bessel_limit_convergence_n128.txt");
            //BesselLimitConvergenceSummary<Pow2.N256>("../../../../results/bessel_limit_convergence_n256.txt");

            //BesselLimitValuesSummary<Pow2.N4>("../../../../results/bessel_limit_values_n4.txt", 60, 56);
            //BesselLimitValuesSummary<Pow2.N8>("../../../../results/bessel_limit_values_n8.txt", 104, 100);
            //BesselLimitValuesSummary<Pow2.N16>("../../../../results/bessel_limit_values_n16.txt", 192, 172);
            //BesselLimitValuesSummary<Pow2.N32>("../../../../results/bessel_limit_values_n32.txt", 372, 368);
            //BesselLimitValuesSummary<Pow2.N64>("../../../../results/bessel_limit_values_n64.txt", 724, 700);
            //BesselLimitValuesSummary<Pow2.N128>("../../../../results/bessel_limit_values_n128.txt", 1436, 1428);
            //BesselLimitValuesSummary<Pow2.N256>("../../../../results/bessel_limit_values_n256.txt", 2856, 2826);

            //BesselNearZeroConvergenceSummary<Pow2.N4, Expand25<Pow2.N4>>("../../../../results/bessel_nearzero_convergence_n4_m5.txt", 60, 56);
            //BesselNearZeroConvergenceSummary<Pow2.N4, Expand50<Pow2.N4>>("../../../../results/bessel_nearzero_convergence_n4_m6.txt", 60, 56);
            //BesselNearZeroConvergenceSummary<Pow2.N4, Double<Pow2.N4>>("../../../../results/bessel_nearzero_convergence_n4_m8_zoom.txt", 60, 56);

            //BesselNearZeroConvergenceSummary<Pow2.N8, Expand25<Pow2.N8>>("../../../../results/bessel_nearzero_convergence_n8_m10.txt", 104, 100);
            //BesselNearZeroConvergenceSummary<Pow2.N8, Expand50<Pow2.N8>>("../../../../results/bessel_nearzero_convergence_n8_m12.txt", 104, 100);
            //BesselNearZeroConvergenceSummary<Pow2.N8, Double<Pow2.N8>>("../../../../results/bessel_nearzero_convergence_n8_m16_zoom.txt", 104, 100);

            //BesselNearZeroConvergenceSummary<Pow2.N16, Expand25<Pow2.N16>>("../../../../results/bessel_nearzero_convergence_n16_m20.txt", 192, 172);
            //BesselNearZeroConvergenceSummary<Pow2.N16, Expand50<Pow2.N16>>("../../../../results/bessel_nearzero_convergence_n16_m24.txt", 192, 172);
            //BesselNearZeroConvergenceSummary<Pow2.N16, Double<Pow2.N16>>("../../../../results/bessel_nearzero_convergence_n16_m32_zoom.txt", 192, 172);

            //BesselNearZeroConvergenceSummary<Pow2.N32, Expand25<Pow2.N32>>("../../../../results/bessel_nearzero_convergence_n32_m40.txt", 372, 368);
            //BesselNearZeroConvergenceSummary<Pow2.N32, Expand50<Pow2.N32>>("../../../../results/bessel_nearzero_convergence_n32_m48.txt", 372, 368);
            //BesselNearZeroConvergenceSummary<Pow2.N32, Double<Pow2.N32>>("../../../../results/bessel_nearzero_convergence_n32_m64_zoom.txt", 372, 368);

            //BesselNearZeroConvergenceSummary<Pow2.N64, Expand25<Pow2.N64>>("../../../../results/bessel_nearzero_convergence_n64_m80.txt", 724, 700);
            //BesselNearZeroConvergenceSummary<Pow2.N64, Expand50<Pow2.N64>>("../../../../results/bessel_nearzero_convergence_n64_m96.txt", 724, 700);
            //BesselNearZeroConvergenceSummary<Pow2.N64, Double<Pow2.N64>>("../../../../results/bessel_nearzero_convergence_n64_m128_zoom.txt", 724, 700);

            //BesselNearZeroConvergenceSummary<Pow2.N128, Expand25<Pow2.N128>>("../../../../results/bessel_nearzero_convergence_n128_m160.txt", 1436, 1428);
            //BesselNearZeroConvergenceSummary<Pow2.N128, Expand50<Pow2.N128>>("../../../../results/bessel_nearzero_convergence_n128_m196.txt", 1436, 1428);
            //BesselNearZeroConvergenceSummary<Pow2.N128, Double<Pow2.N128>>("../../../../results/bessel_nearzero_convergence_n128_m256_zoom.txt", 1436, 1428);

            //BesselNearZeroConvergenceSummary<Pow2.N256, Expand25<Pow2.N256>>("../../../../results/bessel_nearzero_convergence_n256_m320.txt", 2856, 2826);
            //BesselNearZeroConvergenceSummary<Pow2.N256, Expand50<Pow2.N256>>("../../../../results/bessel_nearzero_convergence_n256_m384.txt", 2856, 2826);
            //BesselNearZeroConvergenceSummary<Pow2.N256, Double<Pow2.N256>>("../../../../results/bessel_nearzero_convergence_n256_m512.txt", 2856, 2826);

            BesselJConvergenceSummary<Pow2.N4>("../../../../results_disused/bessel_j_n4.txt");
            BesselYConvergenceSummary<Pow2.N4>("../../../../results_disused/bessel_y_n4.txt");

            BesselJConvergenceSummary<Pow2.N8>("../../../../results_disused/bessel_j_n8.txt");
            BesselYConvergenceSummary<Pow2.N8>("../../../../results_disused/bessel_y_n8.txt");

            BesselJConvergenceSummary<Pow2.N16>("../../../../results_disused/bessel_j_n16.txt");
            BesselYConvergenceSummary<Pow2.N16>("../../../../results_disused/bessel_y_n16.txt");
            
            BesselJConvergenceSummary<Pow2.N32>("../../../../results_disused/bessel_j_n32.txt");
            BesselYConvergenceSummary<Pow2.N32>("../../../../results_disused/bessel_y_n32.txt");
            
            BesselJConvergenceSummary<Pow2.N64>("../../../../results_disused/bessel_j_n64.txt");
            BesselYConvergenceSummary<Pow2.N64>("../../../../results_disused/bessel_y_n64.txt");

            //BesselNonIntegerSummary<Pow2.N8>("../../../../results/bessel_y_n8_nonint.txt");
            //BesselNonIntegerSummary<Pow2.N16>("../../../../results/bessel_y_n16_nonint.txt");
            //BesselNonIntegerSummary<Pow2.N32>("../../../../results/bessel_y_n32_nonint.txt");

            Console.WriteLine("END");
            Console.Read();
        }

        private static void BesselLimitConvergenceSummary<N>(string filepath) where N : struct, IConstant {
            using (StreamWriter sw = new StreamWriter(filepath)) {
                sw.WriteLine($"bits: {MultiPrecision<N>.Bits}");
                sw.WriteLine("nu,z,terms");

                for (decimal nu = 0; nu <= 4; nu += 1 / 32m) {
                    Console.WriteLine(nu);

                    int z, terms;
                    for (z = 1; z <= 65536; z *= 2) {
                        terms = BesselLimitApprox<N>.BesselTermConvergence(nu, z);

                        Console.WriteLine($"{z},{terms}");

                        if (terms < int.MaxValue) {
                            break;
                        }
                    }
                    while (z >= 4) {
                        z = z * 7 / 8;

                        terms = BesselLimitApprox<N>.BesselTermConvergence(nu, z);

                        Console.WriteLine($"{z},{terms}");

                        if (terms >= int.MaxValue) {
                            break;
                        }
                    }
                    z = z / 4 * 4;

                    while (true) {
                        z += 4;

                        terms = BesselLimitApprox<N>.BesselTermConvergence(nu, z);

                        Console.WriteLine($"{z},{terms}");

                        if (terms < int.MaxValue) {
                            break;
                        }
                    }

                    sw.WriteLine($"{nu},{z},{terms}");
                    Console.WriteLine($"{nu},{z},{terms}");
                }
            }
        }

        private static void BesselLimitValuesSummary<N>(string filepath, int z_init, int terms) where N : struct, IConstant {
            using (StreamWriter sw = new StreamWriter(filepath)) {
                sw.WriteLine($"bits: {MultiPrecision<N>.Bits}");

                for (decimal nu = 0; nu <= 2; nu += 1 / 32m) {
                    sw.WriteLine($"nu: {nu}");

                    for (decimal z = z_init; z <= z_init + 4; z += 0.25m) {
                        MultiPrecision<N> t = BesselLimitApprox<N>.Value(nu, z, terms);

                        sw.WriteLine($"  z: {z}");
                        sw.WriteLine($"  {t}");
                        sw.WriteLine($"  {t.ToHexcode()}");
                    }
                }
            }
        }

        private static void BesselNearZeroConvergenceSummary<N, M>(string filepath, int z_init, int terms) where N : struct, IConstant where M : struct, IConstant {
            using (StreamWriter sw = new StreamWriter(filepath)) {
                sw.WriteLine($"bits: {MultiPrecision<N>.Bits}");

                int min_matchbits = MultiPrecision<N>.Bits;

                for (decimal nu = 0; nu <= 2; nu += 1 / 16m) {
                    sw.WriteLine($"nu: {nu}");

                    for (decimal z = z_init - 4; z <= z_init + 4; z += 1 / 32m) {
                        MultiPrecision<N> t = BesselLimitApprox<N>.Value(nu, z, terms);
                        MultiPrecision<N> s = BesselNearZeroApprox<N, M>.Value(nu, z);

                        sw.WriteLine($"  z: {z}");
                        sw.WriteLine($"  LM : {t}");
                        sw.WriteLine($"  NZ : {s}");
                        sw.WriteLine($"  LM : {t.ToHexcode()}");
                        sw.WriteLine($"  NZ : {s.ToHexcode()}");

                        MultiPrecision<N> err = MultiPrecision<N>.Abs(t - s);

                        sw.WriteLine($"  err : {err}");

                        for (int keepbits = MultiPrecision<N>.Bits; keepbits >= 0; keepbits--) {
                            if (keepbits == 0) {
                                min_matchbits = 0;
                            }
                            else if (MultiPrecision<N>.RoundMantissa(t, MultiPrecision<N>.Bits - keepbits) == MultiPrecision<N>.RoundMantissa(s, MultiPrecision<N>.Bits - keepbits)) {
                                sw.WriteLine($"  matchbits : {keepbits}");

                                if (keepbits < min_matchbits) {
                                    min_matchbits = keepbits;
                                }

                                break;
                            }
                        }

                    }
                }

                sw.WriteLine($"min matchbits : {min_matchbits}");
            }
        }

        private static void BesselJConvergenceSummary<N>(string filepath) where N : struct, IConstant {
            using (StreamWriter sw = new StreamWriter(filepath)) {
                sw.WriteLine($"bits: {MultiPrecision<N>.Bits}");

                int z_threshold = (int)Math.Ceiling(50 + 11.0965 * MultiPrecision<N>.Length);

                sw.WriteLine($"z threshold: {z_threshold}");

                for (decimal nu = -64; nu <= 64; nu += 1 / 8m) {
                    int min_matchbits = MultiPrecision<N>.Bits;

                    sw.WriteLine($"nu: {nu}");

                    for (decimal z = 0; z < Math.Min(8 + 1 / 8m, z_threshold - 8); z += 1 / 8m) {
                        min_matchbits = plot(sw, nu, min_matchbits, z);
                    }
                    for (decimal z = Math.Max(0, z_threshold - 8); z <= z_threshold + 8; z += 1 / 8m) {
                        min_matchbits = plot(sw, nu, min_matchbits, z);
                    }

                    sw.WriteLine($"min matchbits : {min_matchbits}");
                }
            }

            static int plot(StreamWriter sw, decimal nu, int min_matchbits, decimal z) {
                MultiPrecision<N> t = MultiPrecisionSandbox<N>.BesselJ(nu, z);
                MultiPrecision<Plus1<N>> s = MultiPrecisionSandbox<Plus1<N>>.BesselJ(nu, z);

                sw.WriteLine($"  z: {z}");
                sw.WriteLine($"  f : {t}");
                sw.WriteLine($"  {t.ToHexcode()}");
                sw.WriteLine($"  {s.ToHexcode()}");

                MultiPrecision<N> err = MultiPrecision<N>.Abs(t - s.Convert<N>());

                sw.WriteLine($"  err : {err}");

                for (int keepbits = MultiPrecision<N>.Bits; keepbits >= 0; keepbits--) {
                    if (keepbits == 0) {
                        min_matchbits = 0;
                    }
                    else if (MultiPrecision<N>.RoundMantissa(t, MultiPrecision<N>.Bits - keepbits) == MultiPrecision<N>.RoundMantissa(s.Convert<N>(), MultiPrecision<N>.Bits - keepbits)) {
                        sw.WriteLine($"  matchbits : {keepbits}");

                        if (keepbits < min_matchbits) {
                            min_matchbits = keepbits;
                        }

                        break;
                    }
                }

                return min_matchbits;
            }
        }

        private static void BesselYConvergenceSummary<N>(string filepath) where N : struct, IConstant {
            using (StreamWriter sw = new StreamWriter(filepath)) {
                sw.WriteLine($"bits: {MultiPrecision<N>.Bits}");

                int z_threshold = (int)Math.Ceiling(50 + 11.0965 * MultiPrecision<N>.Length);

                sw.WriteLine($"z threshold: {z_threshold}");

                for (decimal nu = -64; nu <= 64; nu += 1 / 8m) {
                    int min_matchbits = MultiPrecision<N>.Bits;

                    sw.WriteLine($"nu: {nu}");

                    for (decimal z = 0; z < Math.Min(8 + 1 / 8m, z_threshold - 8); z += 1 / 8m) {
                        min_matchbits = plot(sw, nu, min_matchbits, z);
                    }
                    for (decimal z = Math.Max(0, z_threshold - 8); z <= z_threshold + 8; z += 1 / 8m) {
                        min_matchbits = plot(sw, nu, min_matchbits, z);
                    }

                    sw.WriteLine($"min matchbits : {min_matchbits}");
                }
            }

            static int plot(StreamWriter sw, decimal nu, int min_matchbits, decimal z) {
                MultiPrecision<N> t = MultiPrecisionSandbox<N>.BesselY(nu, z);
                MultiPrecision<Plus1<N>> s = MultiPrecisionSandbox<Plus1<N>>.BesselY(nu, z);

                sw.WriteLine($"  z: {z}");
                sw.WriteLine($"  f : {t}");
                sw.WriteLine($"  {t.ToHexcode()}");
                sw.WriteLine($"  {s.ToHexcode()}");

                MultiPrecision<N> err = MultiPrecision<N>.Abs(t - s.Convert<N>());

                sw.WriteLine($"  err : {err}");

                for (int keepbits = MultiPrecision<N>.Bits; keepbits >= 0; keepbits--) {
                    if (keepbits == 0) {
                        min_matchbits = 0;
                    }
                    else if (MultiPrecision<N>.RoundMantissa(t, MultiPrecision<N>.Bits - keepbits) == MultiPrecision<N>.RoundMantissa(s.Convert<N>(), MultiPrecision<N>.Bits - keepbits)) {
                        sw.WriteLine($"  matchbits : {keepbits}");

                        if (keepbits < min_matchbits) {
                            min_matchbits = keepbits;
                        }

                        break;
                    }
                }

                return min_matchbits;
            }
        }


        private static void BesselYNonIntegerSummary<N>(string filepath) where N : struct, IConstant {
            using (StreamWriter sw = new StreamWriter(filepath)) {
                sw.WriteLine($"bits: {MultiPrecision<N>.Bits}");

                int z_threshold = (int)Math.Ceiling(50 + 11.0965 * MultiPrecision<N>.Length);

                sw.WriteLine($"z threshold: {z_threshold}");

                MultiPrecision<N>[] test_nu = new MultiPrecision<N>[] {
                    MultiPrecision<N>.Ldexp(1, -1), 
                    MultiPrecision<N>.Ldexp(1, -16), 
                    MultiPrecision<N>.Ldexp(1, -17), 
                    MultiPrecision<N>.Ldexp(1, -32), 
                    MultiPrecision<N>.Ldexp(1, -33), 
                    MultiPrecision<N>.Ldexp(1, -48), 
                    MultiPrecision<N>.Ldexp(1, -49), 
                    MultiPrecision<N>.Ldexp(1, -80), 
                    MultiPrecision<N>.Ldexp(1, -81), 
                    MultiPrecision<N>.Ldexp(1, -144), 
                    MultiPrecision<N>.Ldexp(1, -145), 
                    MultiPrecision<N>.Ldexp(1, -272)
                };

                for (int n = -60; n <= 60; n += 10) {
                    foreach (MultiPrecision<N> nu in test_nu) {

                        int min_matchbits = MultiPrecision<N>.Bits;
                
                        sw.WriteLine($"nu: {n + nu}, delta nu: {nu} {nu.Exponent}");

                        for (decimal z = 0; z < Math.Min(8 + 1 / 8m, z_threshold - 8); z += 1 / 8m) {
                            min_matchbits = plot(sw, n + nu, min_matchbits, z);
                        }
                        for (decimal z = Math.Max(0, z_threshold - 8); z <= z_threshold + 8; z += 1 / 8m) {
                            min_matchbits = plot(sw, n + nu, min_matchbits, z);
                        }

                        sw.WriteLine($"nu: {n - nu}, delta nu: {nu} {nu.Exponent}");

                        for (decimal z = 0; z < Math.Min(8 + 1 / 8m, z_threshold - 8); z += 1 / 8m) {
                            min_matchbits = plot(sw, n - nu, min_matchbits, z);
                        }
                        for (decimal z = Math.Max(0, z_threshold - 8); z <= z_threshold + 8; z += 1 / 8m) {
                            min_matchbits = plot(sw, n - nu, min_matchbits, z);
                        }

                        sw.WriteLine($"min matchbits : {min_matchbits}");
                    }
                }
            }

            static int plot(StreamWriter sw, MultiPrecision<N> nu, int min_matchbits, decimal z) {
                MultiPrecision<N> t = MultiPrecisionSandbox<N>.BesselJ(nu, z);
                MultiPrecision<Plus1<N>> s = MultiPrecisionSandbox<Plus1<N>>.BesselJ(nu.Convert<Plus1<N>>(), z);

                sw.WriteLine($"  z: {z}");
                sw.WriteLine($"  f : {t}");
                sw.WriteLine($"  {t.ToHexcode()}");
                sw.WriteLine($"  {s.ToHexcode()}");

                MultiPrecision<N> err = MultiPrecision<N>.Abs(t - s.Convert<N>());

                sw.WriteLine($"  err : {err}");

                for (int keepbits = MultiPrecision<N>.Bits; keepbits >= 0; keepbits--) {
                    if (keepbits == 0) {
                        min_matchbits = 0;
                    }
                    else if (MultiPrecision<N>.RoundMantissa(t, MultiPrecision<N>.Bits - keepbits) == MultiPrecision<N>.RoundMantissa(s.Convert<N>(), MultiPrecision<N>.Bits - keepbits)) {
                        sw.WriteLine($"  matchbits : {keepbits}");

                        if (keepbits < min_matchbits) {
                            min_matchbits = keepbits;
                        }

                        break;
                    }
                }

                return min_matchbits;
            }
        }
    }
}
