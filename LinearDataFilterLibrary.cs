using System;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.Diagnostics;

namespace DataFilterGenericLibrary
{
    public class LinearDataFilterLibrary
    {

        #region Public Methods

        #region Moving Average
        public static double[] movingAverage(int q, double[] input)
        {
            int inputLength = input.Length;

            double[] output = new double[inputLength];

            for (int _dataIndex = 0; _dataIndex < inputLength; _dataIndex ++)
            {

                Double _average = input[_dataIndex];

                for (int _shift = 1; _shift <= q / 2; _shift++)
                {
                    int _lIndex = _dataIndex - _shift;
                    int _rIndex = _dataIndex + _shift;

                    if (_rIndex > inputLength - 1)
                    {
                        _rIndex = 2 * (inputLength - 1) - _rIndex;
                    }

                    if (_lIndex < 0)
                    {
                        _lIndex = -_lIndex;
                    }

                    if ((q & 1) == 0 && _shift == q / 2)
                    {
                        _average += 0.5 * (input[_lIndex] + input[_rIndex]);
                    }
                    else
                    {
                        _average += input[_lIndex] + input[_rIndex];
                    }

                }

                _average /= q;
                output[_dataIndex] = _average;
            }
            return output;
        }
        #endregion

        #region Binomial Average
        public static double[] binomialAverage(int q, double[] input)
        {

            int inputLength = input.Length;
            double[] output = new double[inputLength];

            double[] coefficient = new double[q + 1];

            coefficient[0] = 1.0 / Math.Pow(2.0, 2 * q);

            for (int _coefIndex = 1; _coefIndex <= q; _coefIndex++)
            {
                coefficient[_coefIndex] = coefficient[_coefIndex - 1] * (2 * q - _coefIndex + 1) / (double)_coefIndex;
            }

            for (int _dataIndex = 0; _dataIndex < inputLength; _dataIndex++)
            {
                Double _average = input[_dataIndex] * coefficient[q];

                for (int _shift = 1; _shift <= q; _shift++)
                {
                    int _lIndex = _dataIndex - _shift;
                    int _rIndex = _dataIndex + _shift;

                    if (_rIndex > inputLength - 1)
                    {
                        _rIndex = 2 * (inputLength - 1) - _rIndex;
                    }

                    if (_lIndex < 0)
                    {
                        _lIndex = -_lIndex;
                    }

                    _average += coefficient[q - _shift] * (input[_lIndex] + input[_rIndex]);

                }

                output[_dataIndex] = _average;
            }

            return output;
        }
        #endregion

        #region Exponential Smooth
        public static double[] exponentialSmooth(double t, double[] input)
        {
            int inputLength = input.Length;
            double[] output = new double[inputLength];

            output[0] = input[0];

            for (int _dataIndex = 1; _dataIndex < inputLength; _dataIndex++)
            {
                output[_dataIndex] = (t * input[_dataIndex] + (1 - t) * output[_dataIndex - 1]);
            }

            return output;
        }
        #endregion

        #region Fourier Smooth
        public static double[] fourierSmooth(int bandWidth, double[] input)
        {
            int inputLength = input.Length;
            double[] output = new double[inputLength];

            Tuple<double[], double[]> freqOutput = fft(input, inputLength);
            

            int freqLength = freqOutput.Item1.Length;

            for (int _freqIndex = bandWidth + 1; _freqIndex < freqLength - bandWidth - 1; _freqIndex++)
            {
                freqOutput.Item1[_freqIndex] = 0.0;
                freqOutput.Item2[_freqIndex] = 0.0;
            }

            ifft(freqOutput.Item1, freqOutput.Item2, freqLength);

            for (int _dataIndex = 0; _dataIndex < inputLength; _dataIndex++ )
            {
                output[_dataIndex] = freqOutput.Item1[_dataIndex];
            }

            return output;
            
        }
        #endregion

        #region Savitzky-Golay Smooth
        public static double[] savgolSmooth(int nl, int nr, int m, double[] input)
        {
            double[] coeffs = savgolCoefficient(nl, nr, m);

            int inputLength = input.Length;

            double[] output = new double[inputLength];

            for (int _dataIndex = 0; _dataIndex < inputLength; _dataIndex++)
            {
                for (int _shift = -nl; _shift <= nr; _shift++)
                {
                    int _Index = _dataIndex + _shift;

                    if (_Index > inputLength - 1)
                    {
                        _Index = 2 * (inputLength - 1) - _Index;
                    }

                    if (_Index < 0)
                    {
                        _Index = -_Index;
                    }

                    output[_dataIndex] += coeffs[nl + _shift] * input[_Index];
                }
            }
            return output;

        }
        #endregion

        #region Local Regression 
        public static double[] localRegression(int k, int m,  double[] input)
        {
            if ((k % 2) == 0)
            {
                // can workaround to plus 1.
                throw new Exception(@"Bad Args.");
            }


            int half = (k - 1)/2;
            int inputLength = input.Length;

            double[] output = new double[inputLength];
            double[] weights = new double[k];
            double[] coeffs = new double[k];

            int nr, nl;

            for (int _dataIndex = 0; _dataIndex <= half; _dataIndex++)
            {
                nl = _dataIndex;
                nr = k - 1 - nl;

                setLocalRegressionWeights(nl, nr, weights);
                localRegressionCoeffs(nl, nr, m, input, weights, coeffs);

                for (int _shift = -nl; _shift <= nr; _shift++)
                {
                    output[_dataIndex] += coeffs[nl + _shift] * input[_dataIndex + _shift];
                }
            }

            setLocalRegressionWeights(half, half, weights);
            localRegressionCoeffs(half, half, m, input, weights, coeffs);

            for (int _dataIndex = half + 1; _dataIndex < inputLength - half; _dataIndex++)
            {
                for (int _shift = -half; _shift <= half; _shift++)
                {
                    output[_dataIndex] += coeffs[half + _shift] * input[_dataIndex + _shift];
                }
            }

            for (int _dataIndex = inputLength - half; _dataIndex < inputLength; _dataIndex ++ )
            {
                nr = inputLength - 1 - _dataIndex;
                nl = k - 1 - nr;

                setLocalRegressionWeights(nl, nr, weights);
                localRegressionCoeffs(nl, nr, m, input, weights, coeffs);

                for (int _shift = -nl; _shift <= nr; _shift++)
                {
                    output[_dataIndex] += coeffs[nl + _shift] * input[_dataIndex + _shift];
                }

            }
            return output;
        }
        #endregion

        #region Robust Local Regression
        public static double[] robustLocalRegression(int k, int m, double[] input)
        {
            int inputLength = input.Length;

            double[] residual = new double[inputLength];
            double[] output = localRegression(k, m, input);

            for (int _iter = 0; _iter < 5; _iter ++)
            {
                
                try
                {
                    robustLocalRegressionResidualAnalysis(k, m, input, output, residual);
                }
                catch
                {
                    return output;
                }                
            }

            return output;

        }
        #endregion

        #region Spline Smooth
        public static void splineSmooth(double[] input)
        {

        }
        #endregion

        #region RDP method
        public static List<int> rdp(double epsilon, double[] input)
        {

            int inputLength = input.Length;

            List<int> indexList = new List<int>();

            indexList.Add(0);
            indexList.Add(inputLength - 1);

            recursiveRDP(epsilon, 0, inputLength, input, ref indexList);
            indexList.Sort();

            return indexList;

        }
        #endregion

        #endregion

        #region Private Methods

        #region Robust Local Regression Residual Analysis
        private static void robustLocalRegressionResidualAnalysis(int k, int m, double[] input,
            double[] output, double[] residual)
        {
            int inputLength = input.Length;

            if (inputLength != output.Length)
            {
                throw new Exception(@"Bad Args.");
            }

            if ((k % 2) == 0)
            {
                // can workaround to plus 1.
                throw new Exception(@"Bad Args.");
            }


            int half = (k - 1) / 2;

            for (int _dataIndex = 0; _dataIndex < inputLength; _dataIndex++)
            {
                residual[_dataIndex] = input[_dataIndex] - output[_dataIndex];
                output[_dataIndex] = 0;
            }

            // parallel section
            Parallel.For(0, inputLength, _dataIndex =>
            {
                int nr, nl;

                double[] _weights = new double[k];
                double[] _coeffs = new double[k];

                nl = _dataIndex <= half ? _dataIndex : half;
                if (_dataIndex >= inputLength - half)
                {
                    nl = k - inputLength + _dataIndex;
                }
                nr = k - 1 - nl;

                setRobustLocalRegressionWeights(nl, nr, _dataIndex, _weights, residual);
                localRegressionCoeffs(nl, nr, m, input, _weights, _coeffs);

                for (int _shift = -nl; _shift <= nr; _shift++)
                {
                    output[_dataIndex] += _coeffs[nl + _shift] * input[_dataIndex + _shift];
                }
            });
           

        }

        #endregion

        #region Check Power of Two
        private static bool isPowerOfTwo(int N)
        {
            return (N != 0) && ((N & (N - 1)) == 0);
        }
        #endregion

        #region Compute next Power of Two
        private static int nextPowerOfTwo(int N)
        {
            N--;
            N |= (N >> 1);
            N |= (N >> 2);
            N |= (N >> 4);
            N |= (N >> 8);
            N |= (N >> 16);
            return (N + 1);
        }
        #endregion

        #region FFT
        private static Tuple<double[], double[]> fft(double[] input, int N)
        {
            int inputLength = input.Length;

            if (inputLength != N)
            {
                throw new System.Exception(@"Bad Args.");
            }

            int Ntmp = isPowerOfTwo(N) ? N : nextPowerOfTwo(N);


            double[] realOutput = new double[2 * Ntmp];
            double[] imagOutput = new double[2 * Ntmp];


            for (int _dataIndex = inputLength; _dataIndex < Ntmp; _dataIndex++)
            {
                realOutput[_dataIndex] = input[N - 1];
                realOutput[2 * Ntmp - 1 - _dataIndex] = input[N - 1];
            }


            for (int _dataIndex = 0; _dataIndex < inputLength; _dataIndex++)
            {
                realOutput[_dataIndex] = input[_dataIndex];
                realOutput[2 * Ntmp - 1 - _dataIndex] = input[_dataIndex];
            }



            N = 2 * Ntmp;

            int half = N >> 1;
            
            // REVERSE BIT
            int _swapIndex = 0;
            for (int _dataIndex = 0; _dataIndex < N - 1; _dataIndex++)
            {

                if (_dataIndex < _swapIndex)
                {
                    Double tmp;
                    tmp = realOutput[_dataIndex];
                    realOutput[_dataIndex] = realOutput[_swapIndex];
                    realOutput[_swapIndex] = tmp;
                }

                int stepSize = half;
                while (stepSize <= _swapIndex)
                {
                    _swapIndex -= stepSize;
                    stepSize >>= 1;
                }
                _swapIndex += stepSize;
            }

            // FFT
            Double realPart = -1.0;
            Double imagPart = 0.0;

            int _nextLevelIncrement = 1;

            for (int _levelIndex = 1; _levelIndex < N; _levelIndex <<= 1)
            {
                int _currLevelIncrement = _nextLevelIncrement;
                _nextLevelIncrement <<= 1;

                Double realCoeff = 1.0;
                Double imagCoeff = 0.0;

                for (int _currLevelIncrementIndex = 0; _currLevelIncrementIndex < _currLevelIncrement; _currLevelIncrementIndex++)
                {
                    for (int _jumpIndex = _currLevelIncrementIndex; _jumpIndex < N; _jumpIndex += _nextLevelIncrement)
                    {
                        int _currPos = _jumpIndex + _currLevelIncrement;

                        Double realTmp = realCoeff * realOutput[_currPos] - imagCoeff * imagOutput[_currPos];
                        Double imagTmp = realCoeff * imagOutput[_currPos] + imagCoeff * realOutput[_currPos];

                        realOutput[_currPos] = realOutput[_jumpIndex] - realTmp;
                        imagOutput[_currPos] = imagOutput[_jumpIndex] - imagTmp;

                        realOutput[_jumpIndex] += realTmp;
                        imagOutput[_jumpIndex] += imagTmp;
                    }

                    Double tmp = realCoeff * realPart - imagCoeff * imagPart;
                    imagCoeff = realCoeff * imagPart + imagCoeff * realPart;
                    realCoeff = tmp;
                }
                imagPart = Math.Sqrt((1.0 - realPart) / 2.0);
                imagPart = -imagPart;
                realPart = Math.Sqrt((1.0 + realPart) / 2.0);
            }

            for (int i = 0; i < N; i++)
            {
                realOutput[i] /= N;
                imagOutput[i] /= N;
            }

            return Tuple.Create(realOutput, imagOutput);
        }
        #endregion

        #region iFFT
        private static void ifft(double[] realInput, double[] imagInput, int N)
        {
            int half = N >> 1;

            // REVERSE BIT
            int _swapIndex = 0;
            for (int _dataIndex = 0; _dataIndex < N - 1; _dataIndex++)
            {

                if (_dataIndex < _swapIndex)
                {
                    Double tmp;
                    tmp = realInput[_dataIndex];
                    realInput[_dataIndex] = realInput[_swapIndex];
                    realInput[_swapIndex] = tmp;
                    tmp = imagInput[_dataIndex];
                    imagInput[_dataIndex] = imagInput[_swapIndex];
                    imagInput[_swapIndex] = tmp;
                }

                int stepSize = half;
                while (stepSize <= _swapIndex)
                {
                    _swapIndex -= stepSize;
                    stepSize >>= 1;
                }
                _swapIndex += stepSize;
            }

            // FFT
            Double realPart = -1.0;
            Double imagPart = 0.0;

            int _nextLevelIncrement = 1;

            for (int _levelIndex = 1; _levelIndex < N; _levelIndex <<= 1)
            {
                int _currLevelIncrement = _nextLevelIncrement;
                _nextLevelIncrement <<= 1;

                Double realCoeff = 1.0;
                Double imagCoeff = 0.0;

                for (int _currLevelIncrementIndex = 0; _currLevelIncrementIndex < _currLevelIncrement; _currLevelIncrementIndex++)
                {
                    for (int _jumpIndex = _currLevelIncrementIndex; _jumpIndex < N; _jumpIndex += _nextLevelIncrement)
                    {
                        int _currPos = _jumpIndex + _currLevelIncrement;

                        Double realTmp = realCoeff * realInput[_currPos] - imagCoeff * imagInput[_currPos];
                        Double imagTmp = realCoeff * imagInput[_currPos] + imagCoeff * realInput[_currPos];

                        realInput[_currPos] = realInput[_jumpIndex] - realTmp;
                        imagInput[_currPos] = imagInput[_jumpIndex] - imagTmp;

                        realInput[_jumpIndex] += realTmp;
                        imagInput[_jumpIndex] += imagTmp;
                    }

                    Double tmp = realCoeff * realPart - imagCoeff * imagPart;
                    imagCoeff = realCoeff * imagPart + imagCoeff * realPart;
                    realCoeff = tmp;
                }
                imagPart = Math.Sqrt((1.0 - realPart) / 2.0);
                realPart = Math.Sqrt((1.0 + realPart) / 2.0);
            }
        }
        #endregion

        #region LU Decomposition
        private static Tuple<double[][], int[]> ludcmp(double[][] A)
        {
            int n = A.Length - 1;

            int[] pi = new int[n + 1];

            double p = 0, aki = 0, akpi = 0;
            int kp = 0, pik = 0, pikp = 0;

            for (int j = 0; j <= n; j++)
            {
                pi[j] = j;
            }

            for (int k = 0; k <= n; k++)
            {
                p = 0;
                for (int i = k; i <= n; i++)
                {
                    if (Math.Abs(A[i][k]) > p)
                    {
                        p = Math.Abs(A[i][k]);
                        kp = i;
                    }
                }
                if (p == 0)
                {
                    throw new System.Exception(@"Singularity Detected.");
                }

                pik = pi[k];
                pikp = pi[kp];
                pi[k] = pikp;
                pi[kp] = pik;

                for (int i = 0; i <= n; i++)
                {
                    aki = A[k][i];
                    akpi = A[kp][i];
                    A[k][i] = akpi;
                    A[kp][i] = aki;
                }

                for (int i = k + 1; i <= n; i++)
                {
                    A[i][k] = A[i][k] / A[k][k];
                    for (int j = k + 1; j <= n; j++)
                    {
                        A[i][j] = A[i][j] - (A[i][k] * A[k][j]);
                    }
                }
            }
            return Tuple.Create(A, pi);
        }
        #endregion

        #region LU Backward Substitution
        private static double[] lubksb(double[][] A, int[] idx, double[] b)
        {
            int n = A.Length;

            double[] x = new double[n];

            for (int i = 0; i < n; i++)
            {
                x[i] = b[idx[i]];
            }

            for(int i = 1; i < n; i++)
            {
                double sum = x[i];
                for (int j = 0; j < i; j++)
                {
                    sum -= A[i][j] * x[j];
                }
                x[i] = sum;
            }


            x[n - 1] /= A[n - 1][n - 1];

            for (int i =  n - 2; i >=0; i--)
            {
                double sum = x[i];
                for (int j = i + 1; j < n; j++)
                {
                    sum -= A[i][j] * x[j];
                }
                x[i] = sum / A[i][i];
            }

            return x;
        }
        #endregion

        #region Savitzky-Golay Coefficients
        private static double[] savgolCoefficient(int nl, int nr, int m)
        {
            if (nl < 0 || nr < 0 || nl + nr < m)
            {
                throw new System.Exception(@"Bad Args.");
            }

            double[][] A = new double[m + 1][];
            for (int i = 0; i < m + 1; i++)
            {
                A[i] = new double[m + 1];
            }

            double[] b = new double[m + 1];

            double sum;

            for (int i = 0; i <= m; i++)
            {
                for (int j = 0; j <= m; j++)
                {
                    sum = (i == 0 && j == 0) ? 1.0 : 0.0;
                    for (int k = 1; k <= nr; k++)
                    {
                        sum += Math.Pow(k, i + j);
                    }

                    for (int k = 1; k <= nl; k++)
                    {
                        sum += Math.Pow(-k, i + j);
                    }

                    A[i][j] = sum;
                }
            }

            Tuple<double[][], int[]> output = ludcmp(A);

            b[0] = 1;
            b = lubksb(A, output.Item2, b);

            double[] coeffs = new double[nl + nr + 1];

            for (int i = -nl; i <= nr; i++)
            {
                sum = b[0];
                for (int j = 1; j <= m; j++)
                {
                    sum += b[j] * Math.Pow(i, j);
                }
                coeffs[i + nl] = sum;
            }
            return coeffs;
        }
        #endregion

        #region Shortest Distance From A point to Segment
        private static double
            shortestDistanceToSegment(int currIndex, double currValue,
            int headIndex, double headValue, int tailIndex, double tailValue)
        {
            return Math.Abs((headValue - tailValue) * currIndex + (tailIndex - headIndex) *
                currValue + (headIndex * tailValue - headValue * tailIndex)) / Math.Sqrt(Math.Pow((headValue - tailValue), 2)
                + Math.Pow(tailIndex - headIndex, 2));
        }
        #endregion

        #region Recursive RDP sub-routine
        private static void recursiveRDP(double epsilon, int first, int last, double[] input, ref List<int> indexList)
        {
            double dmax = 0.0;
            int index = first;

            for (int _dataIndex = first; _dataIndex < last; _dataIndex++)
            {
                double d = shortestDistanceToSegment(_dataIndex, input[_dataIndex], first, input[first], last - 1, input[last - 1]);

                if (d > dmax)
                {
                    index = _dataIndex;
                    dmax = d;
                }
            }

            if (dmax > epsilon)
            {
                indexList.Add(index);
                recursiveRDP(epsilon, first, index, input, ref indexList);
                recursiveRDP(epsilon, index, last, input, ref indexList);
            }
        }
        #endregion

        #region Local Regression Coefficients
        private static void localRegressionCoeffs(int nl, int nr, int m, double[] input, double[] weights, double[] coeffs)
        {
            if (nl < 0 || nr < 0 || nl + nr < m)
            {
                Console.WriteLine("{0}, {1}, {2}", nl, nr, m);
                throw new Exception(@"Bad Args.");
            }



            double[][] A = new double[m + 1][];
            for (int i = 0; i < m + 1; i++)
            {
                A[i] = new double[m + 1];
            }

            double[] b = new double[m + 1];

            double sum;

            for (int i = 0; i <= m; i++)
            {
                for (int j = 0; j <= m; j++)
                {
                    sum = (i == 0 && j == 0) ? 1.0 : 0.0;
                    for (int k = 1; k <= nr; k++)
                    {
                        sum += weights[k + nl] * Math.Pow(k, i + j);
                    }

                    for (int k = 1; k <= nl; k++)
                    {
                        sum += weights[nl - k] * Math.Pow(-k, i + j);
                    }

                    A[i][j] = sum;
                }
            }

            Tuple<double[][], int[]> output = ludcmp(A);

            b[0] = 1;
            b = lubksb(A, output.Item2, b);

            for (int i = -nl; i <= nr; i++)
            {
                sum = b[0];
                for (int j = 1; j <= m; j++)
                {
                    sum += b[j] * Math.Pow(i, j);
                }
                coeffs[i + nl] = sum * weights[i + nl];
            }
        }
        #endregion

        #region Set Local Regression Weights
        private static void setLocalRegressionWeights(int nl, int nr, double[] weights)
        {
            double dmax = (double)Math.Max(nr, nl) + 1;

            for (int _wIndex = 0; _wIndex < nr + nl + 1; _wIndex++)
            {
                weights[_wIndex] = Math.Pow(1 - Math.Pow((Math.Abs(nl - _wIndex)) / dmax, 3), 3);
            }
        }
        #endregion

        #region set Robust Local Regression Weights
        private static void setRobustLocalRegressionWeights(int nl, int nr, int ind, double[] weights, double[] residual)
        {
            for (int _dataIndex = -nl; _dataIndex <= nr; _dataIndex++)
            {
                weights[nl + _dataIndex] = Math.Abs(residual[ind + _dataIndex]);
            }

            Array.Sort(weights);
            double median = weights[weights.Length / 2];


            for (int _dataIndex = -nl; _dataIndex <= nr; _dataIndex++)
            {
                weights[nl + _dataIndex] =
                    Math.Abs(residual[ind + _dataIndex]) < 6 * median ?
                    Math.Pow((1 - Math.Pow(residual[nl + _dataIndex] / (6 * median), 2)), 2) : 0;
            }
        }
        #endregion

        #endregion
    }

}
