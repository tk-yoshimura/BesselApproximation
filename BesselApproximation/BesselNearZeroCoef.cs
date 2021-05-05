using MultiPrecision;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BesselApproximation {
    public class BesselNearZeroCoef<N> where N : struct, IConstant {
        private readonly MultiPrecision<N> nu;
        private readonly List<MultiPrecision<N>> a_table = new();
        private readonly List<MultiPrecision<N>> c_table = new();

        public BesselNearZeroCoef(MultiPrecision<N> nu) {
            this.nu = nu;

            MultiPrecision<N> a0 = MultiPrecision<N>.Gamma(nu + 1);
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
                MultiPrecision<N> a =
                    a_table.Last() * (checked(4 * k) * (nu + k));

                a_table.Add(a);
                c_table.Add(1 / a);
            }

            return c_table[n];
        }
    }
}
