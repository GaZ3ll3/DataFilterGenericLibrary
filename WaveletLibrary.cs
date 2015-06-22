using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DataFilterGenericLibrary
{
    #region Wavelet Class
    public class Wavelet
    {
        public string waveletMode;
        public Double[] h1;
        public Double[] g1;
        public Double[] h2;
        public Double[] g2;
        public int nCoeff;
        public int offset;
    }
    #endregion

    #region Threshold Enum Type
    public enum Threshold { rigrsure, heursure, sqtwolog, minimaxi };
    #endregion

    #region Padding Enum Type
    public enum Padding {sym, sp0, sp1, zero, per};
    #endregion

    #region Scaling Enum Type
    public enum Scaling {one, sln, mln };
    #endregion

    public class WaveletLibrary
    {
        #region Public Members
        public static double[] 
            wden(ref double[] input , Threshold thres, bool soft, Scaling scale, int n, ref Wavelet w)
        {
            List<double[]> coeff = WaveletLibrary.wavedec(ref input, n, ref w);

            double[] output = new double[2 * ((input.Length + 1) / 2)];

            double thr;

            double[] scal;

            switch (scale)
            {
                case Scaling.one:
                    scal = new double[n];
                    for (int i = 0; i < n; i++)
                    {
                        scal[i] = 1;
                    }
                    break;
                case Scaling.sln:
                    scal = new double[n];
                    double ssc = sln_wnoisest(coeff, 1);
                    Console.WriteLine("ssc : {0}", ssc);
                    for (int i = 0; i < n; i++)
                    {
                        scal[i] = ssc;
                    }
                    break;
                case Scaling.mln:
                    scal = mln_wnoisest(coeff);
                    break;
                default:
                    throw new Exception("Bad Args");
            }


            // rescale



            for (int i = 0; i < n; i++)
            {
                if (thres == Threshold.sqtwolog || thres == Threshold.minimaxi)
                {
                    thr = thrselect(coeff, thres);
                }
                else
                {
                    if (scal[i] < 1e-7 * coeff[i].Max())
                    {
                        thr = 0;
                    }
                    else
                    {
                        double[] tmp = new double[coeff[i].Length];

                        Array.Copy(coeff[i], 0, tmp, 0, tmp.Length);

                        for (int k = 0; k < tmp.Length; k++)
                        {
                            tmp[k] /= scal[i];
                        }
                        thr = thrselect(tmp, thres);
                    }
                }
                for (int j = 0; j < coeff[i].Length; j++)
                {
                    coeff[i][j] = wthr(coeff[i][j], scal[i] * thr, soft);
                }               
            }

            WaveletLibrary.waverec(ref output, coeff, ref w, 0, output.Length);
            return output;
        }

        #endregion


        #region Private Members
        private static double wthr(double x, double thr, bool soft)
        {
            if (soft)
            {
                if (Math.Abs(x) < thr)
                {
                    return 0;
                }
                else
                {
                    return x > 0 ? Math.Abs(x) - thr : thr - Math.Abs(x);
                }
            }
            else
            {
                if (Math.Abs(x) < thr)
                {
                    return 0;
                }
                else
                {
                    return x;
                }
            }
        }

        private static List<double[]>
            wavedec(ref double[] input, int n, ref Wavelet w)
        {
            int s = input.Length;

            List<double[]> c = new List<double[]>();

            //workspaces
            double[] tmp = new double[ s + 2 *(w.g1.Length - 1) + w.g1.Length - 1];
            double[] z = new double[s + 2 * (w.g1.Length - 1) - w.g1.Length + 1];
            //pass by reference's value
            double[] iter = input;

            for (int k = 1; k <= n; ++k)
            {                
                Tuple<double[], double[]> x = dwt( ref z, ref tmp, ref iter, ref w, 0, Padding.sp0);                
                c.Add(x.Item2);               
                iter = x.Item1;                
            }

            c.Add(iter);
            return c;
        }


        private static void
            waverec(ref double[] output, List<double[]> coeff, ref Wavelet w, int level, int n)
        {
            int nmax = coeff.Count() - 1;

            double[] tmp = new double[n];
            double[] g2 = fliplr(w.g1);
            double[] h2 = fliplr(w.h1);

            Array.Copy(coeff[nmax], 0, output, 0, coeff[nmax].Length);

            for (int p = nmax; p > level; p--)
            {
                idwt(ref output, ref tmp, coeff[p - 1], ref g2, ref h2, 0, Padding.sp0);
            }

        }

        private static double thrselect(double[] coeff, Threshold thres)
        {
            int n = coeff.Length;
            double thr;

            switch (thres)
            {
                case Threshold.rigrsure:
                    double[] sx2 = new double[n];

                    int p = 0;


                    for (int j = 0; j < n; j++)
                    {
                        sx2[p++] = coeff[j] * coeff[j];
                    }
                    

                    Array.Sort(sx2);

                    double min_risk = Double.MaxValue;
                    double curr_risk;
                    double cumsum = 0;
                    p = 0;

                    for (int i = 1; i <= n; i++)
                    {
                        cumsum += sx2[i - 1];
                        curr_risk = (n - 2 * (i) + cumsum + (n - i) * sx2[i - 1]) / (double)n;

                        if (curr_risk < min_risk)
                        {
                            min_risk = curr_risk;
                            p = i - 1;
                        }
                    }

                    thr = Math.Sqrt(sx2[p]);

                    break;
                case Threshold.minimaxi:
                    if (n <= 32)
                    {
                        thr = 0;
                    }
                    else
                    {
                        thr = 0.3936 + 0.1829 * (Math.Log(n) / Math.Log(2));
                    }
                    break;
                case Threshold.heursure:

                    double hthr = Math.Sqrt(2 * Math.Log(n));
                    double eta = 0;

                    for (int j = 0; j < coeff.Length; j++)
                    {
                        eta += coeff[j] * coeff[j];
                    }
                    

                    eta -= n;
                    eta /= (double)n;

                    double crit = Math.Pow(Math.Log(n) / Math.Log(2), 1.5) / Math.Sqrt(n);

                    if (eta < crit)
                    {
                        thr = hthr;
                    }
                    else
                    {
                        thr = Math.Min(thrselect(coeff, Threshold.rigrsure), hthr);
                    }

                    break;
                case Threshold.sqtwolog:
                    thr = Math.Sqrt(2 * Math.Log(n));
                    break;
                default:
                    throw new Exception("Bad Args");


            }
            //Console.Write("thr: {0}", thr);
            return thr;
        }

        private static double thrselect(List<double[]> coeff, Threshold thres)
        {
            double thr;

            int n = 0;

            foreach (double[] x in coeff)
            {
                n += x.Length;
            }

            //Console.WriteLine("coeff len : {0} ", n);
            switch (thres)
            {
                case Threshold.rigrsure:
                    double[] sx2 = new double[n];

                    int p = 0;

                    for (int i = 0; i < coeff.Count(); i++)
                    {
                        for (int j = 0; j < coeff[i].Length; j++)
                        {
                            sx2[p++] = coeff[i][j] * coeff[i][j];
                        }
                    }

                    Array.Sort(sx2);

                    double min_risk = Double.MaxValue;
                    double curr_risk;
                    double cumsum = 0;
                    p = 0;

                    for (int i = 1; i <= n; i++)
                    {
                        cumsum += sx2[i - 1];
                        curr_risk = (n - 2 * (i) + cumsum + (n - i) * sx2[i - 1]) / (double)n;

                        if (curr_risk < min_risk)
                        {
                            min_risk = curr_risk;
                            p = i - 1;
                        }
                    }

                    thr = Math.Sqrt(sx2[p]);

                    break;
                case Threshold.minimaxi:
                    if (n <= 32)
                    {
                        thr = 0;
                    }
                    else
                    {
                        thr = 0.3936 + 0.1829 * (Math.Log(n) / Math.Log(2));
                    }
                    break;
                case Threshold.heursure:

                    double hthr = Math.Sqrt(2 * Math.Log(n));
                    double eta = 0;
                    for (int i = 0; i < coeff.Count(); i++)
                    {
                        for (int j = 0; j < coeff[i].Length; j++)
                        {
                            eta += coeff[i][j] * coeff[i][j];
                        }
                    }

                    eta -= n;
                    eta /= (double)n;

                    double crit = Math.Pow(Math.Log(n)/ Math.Log(2), 1.5)/Math.Sqrt(n);

                    if (eta < crit)
                    {
                        thr = hthr;
                    }
                    else 
                    {
                        thr = Math.Min(thrselect(coeff, Threshold.rigrsure), hthr);
                    }

                    break;
                case Threshold.sqtwolog:
                    thr = Math.Sqrt(2 * Math.Log(n));
                    break;
                default:
                    throw new Exception("Bad Args");
                    
                    
            }
            //Console.Write("thr: {0}", thr);
            return thr;
        }

        private static double sln_wnoisest(List<double[]> coeff, int n)
        {
            int lx = coeff.Count() - 1 ;

            if (n > lx  || n < 1)
            {
                throw new Exception("Bad Args");
            }

            int ls = coeff[n - 1].Length;

            double[] span = new double[ls];

            for (int i = 0; i < coeff[n - 1].Length; i++)
            {
                span[i] = Math.Abs(coeff[n - 1][i]);
            }



            Array.Sort(span);

            
            return ls % 2 == 0 ? (span[(ls) / 2 - 1] + span[(ls) / 2 ]) / 2 / 0.6745 : span[(ls + 1) / 2 - 1] / 0.6745;
        }

        private static double[] mln_wnoisest(List<double[]> coeff)
        {
            int lx = coeff.Count() - 1;
            double[] x = new double[lx];

            for (int i = 1; i <= lx; i++)
            {
                int ls = coeff[i - 1].Length;
                double[] span = new double[ls];

                for (int j = 0; j < coeff[i - 1].Length; j++)
                {
                    span[j] = Math.Abs(coeff[i - 1][j]);
                }

                Array.Sort(span);

                x[i - 1] = ls % 2 == 0 ? (span[(ls) / 2 - 1] + span[(ls) / 2 ]) / 2 / 0.6745 : span[(ls + 1) / 2 - 1] / 0.6745;
            }

            return x;
        }

        private static Tuple<double[], double[]>
            dwt(ref double[] z, ref double[] tmp, ref double[] input, ref Wavelet w, int shift, Padding pad)
        {
            int lf = w.g1.Length;
            int lx = input.Length;
            
            shift = shift % 2; // usually it is 0

            int first = 2 - shift;

            bool flagPer = (pad == Padding.per);

            int lenEXT, last;

            if (!flagPer)
            {
                lenEXT = lf - 1; last = lx + lf - 1;
            }
            else 
            {
                lenEXT = lf / 2; last = lx % 2 == 0 ? lx : lx - 1; 
            }

            double[] a_coeff = new double[1 + (last - first) / 2];
            double[] d_coeff = new double[1 + (last - first) / 2];
            double[] y = new double[lx + 2 * (w.g1.Length - 1)];

            wextend(ref y, ref input, pad, lenEXT);

            int p = 0;

            conv2(ref z, ref tmp, ref y, ref w.h1, y.Length, 1, w.h1.Length, 1, "valid");
            
            for (int i = first - 1; i <= last - 1; i += 2)
            {
                d_coeff[p++] = z[i];
            }

            conv2(ref z , ref tmp, ref y, ref w.g1, y.Length, 1, w.g1.Length, 1, "valid");
            
            p = 0;
            for (int i = first - 1; i <= last - 1; i += 2)
            {
                a_coeff[p++] = z[i];
            }

            return Tuple.Create(a_coeff, d_coeff);
        }

        //to be efficient, app is as long as the length of input at least.
        private static void idwt(ref double[] app, ref double[] tmp, double[] det, ref double[] g2, ref double[] h2, int shift, Padding pad)
        {

            bool flagPer = (pad == Padding.per);
            if (flagPer)
            {
                throw new Exception("Not Implmented");
            }

            int s = flagPer ? 2 * det.Length : 2 * det.Length - g2.Length + 2;

            
            if (tmp.Length < s)
            {
                throw new Exception("Bad Args");
            }
            //clean up
            for (int i = 0; i < s; i++)
            {
                tmp[i] = 0;
            }

            upsconv1(ref tmp, ref app, ref g2, det.Length, shift, pad);
            upsconv1(ref tmp, ref det, ref h2, det.Length, shift, pad);

            Array.Copy(tmp,0 , app, 0, s);
        }

        private static double[] fliplr(double[] x)
        {
            int lx = x.Length;
            double[] y = new double[lx];
            for (int i = 0; i < lx; i++)
            {
                y[i] = x[lx - 1 - i];
            }

            return y;
        }

        private static void upsconv1(ref double[] res, ref double[] x, ref double[] f, int _lx, int shift, Padding pad)
        {
             //edit app.
            bool flagPer = (pad == Padding.per);

            if (flagPer) 
            {
                throw new Exception("Not Implemented");
            }
            int lx = 2 * _lx;
            int lf = f.Length;

            int s = flagPer ? lx : lx - lf + 2;

            int lz = x.Length % 2==0 ? lx : lx - 1;

            double[] z = new double[lz];
            double[] tmp = new double[lz + lf - 1];

            for (int i = 0 ; i < lz ; i += 2)
            {
                z[i] = x[i/2];
            }
            _conv2(ref tmp, ref z, ref f, lz, 1, lf, 1, 1);

            int sx = tmp.Length;
            double d = (double)(sx - s)/ 2.0;

            int first = (int)Math.Floor(d);
            int last = sx - 1 - (int)Math.Ceiling(d);


            if (shift == 1)
            {
                first = (int)Math.Ceiling(d);
                last = sx -1 - (int)Math.Floor(d);
            }

            for (int i = 0; i <= last - first; i++)
            {
                res[i] += tmp[first + i];
            }
        }


        private static void wextend(ref double[] y, ref double[] x, Padding pad, int lenEXT)
        {
            //todo
            int lx = x.Length;
            switch (pad)
            {
                case Padding.sp0:
                    getSP0_ext(ref y, x[0], x[lx - 1], lenEXT);
                    Array.Copy(x, 0, y, lenEXT, lx);
                    break;
            }
        }

        // 1d method, always both sides
        private static void getSP0_ext(ref double[] y, double x_l, double x_r, int lf)
        {
            int n = y.Length;
            for (int i = 0; i < lf; i++)
            {
                y[i] = x_l;
                y[n - 1 - i] = x_r;
            }
        }

        private static void
            _conv2(ref double[] c, ref double[] a, ref double[] b, int ma , int na, int mb, int nb, int plusminus)
        {
            int p, q;
            double w;
            int mc, nc;
            int k, l, i, j;
            int r;

            int flops = 0;

            mc = ma + mb - 1;
            nc = na + nb - 1;

            for (j = 0; j < nc; j++)
            {
                for (i = 0; i < mc; i++)
                {
                    c[i + j * mc] = 0;
                }
            }

                r = 0;
            for (j = 0; j < nb; ++j)
            {
                for (i = 0; i < mb; ++i)
                {
                    w = b[r++];

                    if (w != 0.0)
                    {
                        p = i + j * mc;
                        

                        for (l = 0, q = 0; l < na; l++)
                        {
                            for (k = 0; k < ma; k++)
                            {
                                c[p++] += a[q++] * w * plusminus;
                            }

                            p += mb - 1;
                        }
                        flops += 2 * ma * na;
                    }
                }
            }
        }

        private static void
            subMatrix(ref double[] b, ref double[] a, int ma, int na, int rstart, int rend, int cstart, int cend)
        {
            int p, q;
            int i, j;
            int mb, nb, step;

            mb = rend - rstart + 1;
            nb = cend - cstart + 1;

            step = ma - mb;

            q = 0;
            p = rstart + cstart * ma;

            for (j = 0; j < nb; ++j)
            {
                for (i = 0; i < mb; ++i)
                {
                    b[q++] = a[p++];
                }
                p += step;
            }
        }

        private static void
            conv2(ref double[] c, ref double[] tmp, ref double[] a, ref double[] b, int ma, int na, int mb, int nb, string shape)
        {
            int switched;

            int mc;
            int nc;

            int tmpm = ma + mb - 1;
            int tmpn = na + nb - 1;
            int code ; 

            switch (shape[0])
            {
                case 's':
                    code = 0;
                    break;
                case 'f':
                    code = 1;
                    break;
                case 'v':
                    code = 2;
                    break;
                default:
                    throw new Exception("Bad Args");
            }

            if (ma * na > mb * nb)
            {
                _conv2(ref tmp,ref  a, ref b, ma, na, mb, nb, 1);

                switched = 0;
            }
            else
            {
                _conv2(ref tmp, ref b, ref a, mb, nb, ma, na, 1);
                switched = 1;
            }

            switch (code)
            {
                case 0:
                    Array.Copy(tmp, 0, c, 0, Math.Min(tmp.Length, c.Length));
                    break;

                case 1:
                    if (switched == 1)
                    {
                        mc = mb;
                        nc = nb;
                        mb = ma;
                        nb = na;
                    }
                    else
                    {
                        mc = ma;
                        nc = na;
                    }

                    subMatrix(ref c, ref tmp, tmpm, tmpn, mb / 2, mb / 2 + mc - 1, nb / 2, nb / 2 + nc - 1);

                    break;
                case 2:
                    mc = ma - mb + 1;
                    nc = na - nb + 1;

                    subMatrix(ref c, ref tmp, tmpm, tmpn, mb - 1, mb + mc - 2, nb  - 1, nb  + nc - 2);
                    break;
            }
        }


        #endregion

    }

    #region Daubechies Wavelet Settings
    public class Daubechies
    {
        /// <summary>
        /// daubechies_init
        /// 
        /// usage:
        /// initialize the wavelet struct with Daubechies 
        /// wavelets.
        /// </summary>
        /// <param name="member"></param>
        /// <param name="centered"></param>
        /// <returns></returns>
        public static Wavelet daubechies_init(int member, bool centered)
        {
            Wavelet wavelet = new Wavelet();

            switch (member)
            {
                case 4:
                    wavelet.waveletMode = "db4";
                    wavelet.h1 = h_4;
                    wavelet.g1 = g_4;
                    wavelet.h2 = h_4;
                    wavelet.g2 = g_4;
                    break;

                case 6:
                    wavelet.waveletMode = "db6";
                    wavelet.h1 = h_6;
                    wavelet.g1 = g_6;
                    wavelet.h2 = h_6;
                    wavelet.g2 = g_6;
                    break;

                case 8:
                    wavelet.waveletMode = "db8";
                    wavelet.h1 = h_8;
                    wavelet.g1 = g_8;
                    wavelet.h2 = h_8;
                    wavelet.g2 = g_8;
                    break;

                case 10:
                    wavelet.waveletMode = "db10";
                    wavelet.h1 = h_10;
                    wavelet.g1 = g_10;
                    wavelet.h2 = h_10;
                    wavelet.g2 = g_10;
                    break;

                case 12:
                    wavelet.waveletMode = "db12";
                    wavelet.h1 = h_12;
                    wavelet.g1 = g_12;
                    wavelet.h2 = h_12;
                    wavelet.g2 = g_12;
                    break;

                case 14:
                    wavelet.waveletMode = "db14";
                    wavelet.h1 = h_14;
                    wavelet.g1 = g_14;
                    wavelet.h2 = h_14;
                    wavelet.g2 = g_14;
                    break;

                case 16:
                    wavelet.waveletMode = "db16";
                    wavelet.h1 = h_16;
                    wavelet.g1 = g_16;
                    wavelet.h2 = h_16;
                    wavelet.g2 = g_16;
                    break;

                case 18:
                    wavelet.waveletMode = "db18";
                    wavelet.h1 = h_18;
                    wavelet.g1 = g_18;
                    wavelet.h2 = h_18;
                    wavelet.g2 = g_18;
                    break;

                case 20:
                    wavelet.waveletMode = "db20";
                    wavelet.h1 = h_20;
                    wavelet.g1 = g_20;
                    wavelet.h2 = h_20;
                    wavelet.g2 = g_20;
                    break;

                default:
                    wavelet.waveletMode = "db4";
                    wavelet.h1 = h_4;
                    wavelet.g1 = g_4;
                    wavelet.h2 = h_4;
                    wavelet.g2 = g_4;
                    break;
            }


            wavelet.nCoeff = member;
            if (!centered)
            {
                wavelet.offset = 0;
            }
            else
            {
                wavelet.offset = (member >> 1);
            }

            return wavelet;

        }

        public static double[] h_4 = new double[]
        { 
            0.48296291314453414337487159986,
            0.83651630373780790557529378092,
            0.22414386804201338102597276224,
            -0.12940952255126038117444941881
        };

        public static double[] g_4 = new double[]
        { 
            -0.12940952255126038117444941881,
            -0.22414386804201338102597276224,
            0.83651630373780790557529378092,
            -0.48296291314453414337487159986
        };

        public static double[] h_6 = new double[] 
        { 
            -0.33267055295008261599851158914,
            0.80689150931109257649449360409,
            -0.45987750211849157009515194215,
            -0.13501102001025458869638990670,
            0.08544127388202666169281916918,
            0.03522629188570953660274066472
        };

        public static double[] g_6 = new double[]
        { 
            0.03522629188570953660274066472,
            -0.08544127388202666169281916918,
            -0.13501102001025458869638990670,
            0.45987750211849157009515194215,
            0.80689150931109257649449360409,
            0.33267055295008261599851158914
        };

        public static double[] h_8 = new double[]
        {
            0.23037781330889650086329118304,
            0.71484657055291564708992195527,
            0.63088076792985890788171633830,
            -0.02798376941685985421141374718,
            -0.18703481171909308407957067279,
            0.03084138183556076362721936253,
            0.03288301166688519973540751355,
            -0.01059740178506903210488320852
        };

        public static double[] g_8 = new double[]
        {
            -0.01059740178506903210488320852,
            -0.03288301166688519973540751355,
            0.03084138183556076362721936253,
            0.18703481171909308407957067279,
            -0.02798376941685985421141374718,
            -0.63088076792985890788171633830,
            0.71484657055291564708992195527,
            -0.23037781330889650086329118304
        };

        public static double[] h_10 = new double[] 
        {
            0.16010239797419291448072374802,
            0.60382926979718967054011930653,
            0.72430852843777292772807124410,
            0.13842814590132073150539714634,
            -0.24229488706638203186257137947,
            -0.03224486958463837464847975506,
            0.07757149384004571352313048939,
            -0.00624149021279827427419051911,
            -0.01258075199908199946850973993,
            0.00333572528547377127799818342
        };

        public static double[] g_10 = new double[] 
        { 
            0.00333572528547377127799818342,
            0.01258075199908199946850973993,
            -0.00624149021279827427419051911,
            -0.07757149384004571352313048939,
            -0.03224486958463837464847975506,
            0.24229488706638203186257137947,
            0.13842814590132073150539714634,
            -0.72430852843777292772807124410,
            0.60382926979718967054011930653,
            -0.16010239797419291448072374802
        };

        public static double[] h_12 = new double[]
        { 
            0.11154074335010946362132391724,
            0.49462389039845308567720417688,
            0.75113390802109535067893449844,
            0.31525035170919762908598965481,
            -0.22626469396543982007631450066,
            -0.12976686756726193556228960588,
            0.09750160558732304910234355254,
            0.02752286553030572862554083950,
            -0.03158203931748602956507908070,
            0.00055384220116149613925191840,
            0.00477725751094551063963597525,
            -0.00107730108530847956485262161
        };

        public static double[] g_12 = new double[] 
        { 
            -0.00107730108530847956485262161,
            -0.00477725751094551063963597525,
            0.00055384220116149613925191840,
            0.03158203931748602956507908070,
            0.02752286553030572862554083950,
            -0.09750160558732304910234355254,
            -0.12976686756726193556228960588,
            0.22626469396543982007631450066,
            0.31525035170919762908598965481,
            -0.75113390802109535067893449844,
            0.49462389039845308567720417688,
            -0.11154074335010946362132391724
        };

        public static double[] h_14 = new double[]
        { 
            0.07785205408500917901996352196,
            0.39653931948191730653900039094,
            0.72913209084623511991694307034,
            0.46978228740519312247159116097,
            -0.14390600392856497540506836221,
            -0.22403618499387498263814042023,
            0.07130921926683026475087657050,
            0.08061260915108307191292248036,
            -0.03802993693501441357959206160,
            -0.01657454163066688065410767489,
            0.01255099855609984061298988603,
            0.00042957797292136652113212912,
            -0.00180164070404749091526826291,
            0.00035371379997452024844629584
        };

        public static double[] g_14 = new double[] 
        { 
            0.00035371379997452024844629584,
            0.00180164070404749091526826291,
            0.00042957797292136652113212912,
            -0.01255099855609984061298988603,
            -0.01657454163066688065410767489,
            0.03802993693501441357959206160,
            0.08061260915108307191292248036,
            -0.07130921926683026475087657050,
            -0.22403618499387498263814042023,
            0.14390600392856497540506836221,
            0.46978228740519312247159116097,
            -0.72913209084623511991694307034,
            0.39653931948191730653900039094,
            -0.07785205408500917901996352196
        };

        public static double[] h_16 = new double[]
        {
            0.05441584224310400995500940520,
            0.31287159091429997065916237551,
            0.67563073629728980680780076705,
            0.58535468365420671277126552005,
            -0.01582910525634930566738054788,
            -0.28401554296154692651620313237,
            0.00047248457391328277036059001,
            0.12874742662047845885702928751,
            -0.01736930100180754616961614887,
            -0.04408825393079475150676372324,
            0.01398102791739828164872293057,
            0.00874609404740577671638274325,
            -0.00487035299345157431042218156,
            -0.00039174037337694704629808036,
            0.00067544940645056936636954757,
            -0.00011747678412476953373062823
        };

        public static double[] g_16 = new double[] 
        {
            -0.00011747678412476953373062823,
            -0.00067544940645056936636954757,
            -0.00039174037337694704629808036,
            0.00487035299345157431042218156,
            0.00874609404740577671638274325,
            -0.01398102791739828164872293057,
            -0.04408825393079475150676372324,
            0.01736930100180754616961614887,
            0.12874742662047845885702928751,
            -0.00047248457391328277036059001,
            -0.28401554296154692651620313237,
            0.01582910525634930566738054788,
            0.58535468365420671277126552005,
            -0.67563073629728980680780076705,
            0.31287159091429997065916237551,
            -0.05441584224310400995500940520
        };

        public static double[] h_18 = new double[] 
        {
            0.03807794736387834658869765888,
            0.24383467461259035373204158165,
            0.60482312369011111190307686743,
            0.65728807805130053807821263905,
            0.13319738582500757619095494590,
            -0.29327378327917490880640319524,
            -0.09684078322297646051350813354,
            0.14854074933810638013507271751,
            0.03072568147933337921231740072,
            -0.06763282906132997367564227483,
            0.00025094711483145195758718975,
            0.02236166212367909720537378270,
            -0.00472320475775139727792570785,
            -0.00428150368246342983449679500,
            0.00184764688305622647661912949,
            0.00023038576352319596720521639,
            -0.00025196318894271013697498868,
            0.00003934732031627159948068988
        };

        public static double[] g_18 = new double[]
        {
            0.00003934732031627159948068988,
            0.00025196318894271013697498868,
            0.00023038576352319596720521639,
            -0.00184764688305622647661912949,
            -0.00428150368246342983449679500,
            0.00472320475775139727792570785,
            0.02236166212367909720537378270,
            -0.00025094711483145195758718975,
            -0.06763282906132997367564227483,
            -0.03072568147933337921231740072,
            0.14854074933810638013507271751,
            0.09684078322297646051350813354,
            -0.29327378327917490880640319524,
            -0.13319738582500757619095494590,
            0.65728807805130053807821263905,
            -0.60482312369011111190307686743,
            0.24383467461259035373204158165,
            -0.03807794736387834658869765888
        };

        public static double[] h_20 = new double[]
        {
            -0.02667005790055555358661744877,//
            0.18817680007769148902089297368,
            -0.52720118893172558648174482796,//
            0.68845903945360356574187178255,
            -0.28117234366057746074872699845,//
            -0.24984642432731537941610189792,
            0.19594627437737704350429925432,//
            0.12736934033579326008267723320,
            -0.09305736460357235116035228984,//
            -0.07139414716639708714533609308,
            0.02945753682187581285828323760,//
            0.03321267405934100173976365318,
            -0.00360655356695616965542329142,//
            -0.01073317548333057504431811411,
            -0.00139535174705290116578931845,//
            0.00199240529518505611715874224,
            0.00068585669495971162656137098,//
            -0.00011646685512928545095148097,
            -0.00009358867032006959133405013,//
            -0.00001326420289452124481243668
        };

        public static double[] g_20 = new double[]
        { 
            -0.00001326420289452124481243668,
            0.00009358867032006959133405013,//
            -0.00011646685512928545095148097,
            -0.00068585669495971162656137098,//
            0.00199240529518505611715874224,
            0.00139535174705290116578931845,//
            -0.01073317548333057504431811411,
            0.00360655356695616965542329142,//
            0.03321267405934100173976365318,
            -0.02945753682187581285828323760,//
            -0.07139414716639708714533609308,
            0.09305736460357235116035228984,//
            0.12736934033579326008267723320,
            -0.19594627437737704350429925432,//
            -0.24984642432731537941610189792,
            0.28117234366057746074872699845,//
            0.68845903945360356574187178255,
            0.52720118893172558648174482796,//
            0.18817680007769148902089297368,
            0.02667005790055555358661744877//
        };




    }
    #endregion

    #region Haar Wavelet Settings
    public class Haar
    {
        public static Wavelet haar_init(bool centered)
        {
            Wavelet wavelet = new Wavelet();

            wavelet.waveletMode = "haar";
            wavelet.h1 = ch2;
            wavelet.g1 = cg2;
            wavelet.h2 = ch2;
            wavelet.g2 = cg2;

            wavelet.nCoeff = 2;

            if (!centered)
            {
                wavelet.offset = 0;
            }
            else
            {
                wavelet.offset = 1;
            }

            return wavelet;
        }


        public static double[] ch2 = new double[] { 1.0 / Math.Sqrt(2), 1.0 / Math.Sqrt(2) };
        public static double[] cg2 = new double[] { 1.0 / Math.Sqrt(2), -1.0 / Math.Sqrt(2) };
    }
    #endregion
}
