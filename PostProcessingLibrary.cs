using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DataFilterGenericLibrary
{
    public class PostProcessLibrary
    {
        public static double snr(ref double[] signal, ref double[] noise)
        {
            double sr = rssq(ref signal);
            double nr = rssq(ref noise);
            return 10 * Math.Log10((sr * sr) / (nr * nr));
        }

        private static double rssq(ref double[] signal)
        {
            double s = 0;
            int lx = signal.Length;
            for (int i = 0; i < lx; i++)
            {
                s += signal[i] * signal[i];
            }

            return Math.Sqrt(s);
        }
    }
}
