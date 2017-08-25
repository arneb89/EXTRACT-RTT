using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FITS_READER
{
    class Extractor
    {
        private static Image image = null;
        
        public static double[][] Flux = null;
        public static double[][] Background = null;
        public static double[][] FluxDebiased = null;

        public static double[][] fluxes_min = null;
        public static double[][] fluxes_min0 = null;
        public static double[][] fluxes_min1 = null;

        public static void Extract(Image im, double[][] pos_ord, double[][] pos_min, int aperture_wing)
        {
            image = im;
            int pix_count = im.NAXIS2;
            int ord_count = pos_ord.Length;

            int[][] pos_min_int = new int[pos_min.Length][];
            for (int i = 0; i < pos_min_int.Length; i++)
                pos_min_int[i] = new int[pix_count];

            for (int i = 0; i < pos_min_int.Length; i++)
                for (int j = 0; j < pos_min_int[0].Length; j++)
                    pos_min_int[i][j] = (int)Math.Round(pos_min[i][j], 0);

            MinFluxes0(pos_min_int);
            for (int i = 0; i < fluxes_min0.Length; i++)
            {
                fluxes_min0[i] = Filter(fluxes_min0[i], 3);
            }

            Flux = new double[ord_count][];
            for (int i = 0; i < Flux.Length; i++)
                Flux[i] = new double[pix_count];

            Background = new double[ord_count][];
            for (int i = 0; i < Background.Length; i++)
                Background[i] = new double[pix_count];

            FluxDebiased = new double[ord_count][];
            for (int i = 0; i < FluxDebiased.Length; i++)
                FluxDebiased[i] = new double[pix_count];

            int mode = 1;
            Console.Write("Extracted orders: ");
            for (int ord = 0; ord < ord_count; ord++)
            {
                Console.Write("[{0}]", ord);
                for (int i = 0; i < pix_count; i++)
                {
                    int n1, n2;
                    double bkg = 0;

                    n1 = pos_min_int[ord][i];
                    n2 = pos_min_int[ord + 1][i];

                    if (aperture_wing != -1)
                    {
                        int ord_pos_int = (int)Math.Round(pos_ord[ord][i], 0);
                        n1 = ord_pos_int - aperture_wing;
                        n2 = ord_pos_int + aperture_wing;
                    }

                    bkg = 0.5 * (fluxes_min0[ord][i] + fluxes_min0[ord + 1][i]) * (n2 - n1 + 1);

                    Background[ord][i] = bkg;

                    if (mode == 1)
                    {
                        double sum = 0;
                        for (int j = n1; j <= n2; j++)
                        {
                            sum += im[j, i];
                        }
                        Flux[ord][i] = sum;
                    }
                    FluxDebiased[ord][i] = Flux[ord][i] - bkg;
                    if (FluxDebiased[ord][i] < 0) FluxDebiased[ord][i] = 0;
                }
            }
            Console.Write("\r\n");
        }

        public static void ExtractOptSlit(Image im, double[][] pos_ord, double[][] pos_min, int aperture_wing)
        {
            image = im;
            int pix_count = im.NAXIS2;
            int ord_count = pos_ord.Length;

            int[][] pos_min_int = new int[pos_min.Length][];
            for (int i = 0; i < pos_min_int.Length; i++)
                pos_min_int[i] = new int[pix_count];

            for (int i = 0; i < pos_min_int.Length; i++)
                for (int j = 0; j < pos_min_int[0].Length; j++)
                    pos_min_int[i][j] = (int)Math.Round(pos_min[i][j], 0);

            int wing = 5;
            int degree = 3;
            MinFluxes0(pos_min_int);
            MinFluxes1(pos_min_int, wing, degree);
            BackGroundFit(pos_min_int);

            Flux = new double[ord_count][];
            for (int i = 0; i < Flux.Length; i++)
                Flux[i] = new double[pix_count];
            Background = new double[ord_count][];
            for (int i = 0; i < Background.Length; i++)
                Background[i] = new double[pix_count];
            FluxDebiased = new double[ord_count][];
            for (int i = 0; i < FluxDebiased.Length; i++)
                FluxDebiased[i] = new double[pix_count];
            int mode = 1;
            Console.Write("Extracted orders: ");
            for (int ord = 0; ord < ord_count; ord++)
            {
                Console.Write("[{0}]", ord);
                for (int i = 0; i < pix_count; i++)
                {
                    int n1, n2;
                    double bkg = 0;

                    n1 = pos_min_int[ord][i];
                    n2 = pos_min_int[ord + 1][i];

                    if (aperture_wing != -1)
                    {
                        int ord_pos_int = (int)Math.Round(pos_ord[ord][i], 0);
                        n1 = ord_pos_int - aperture_wing;
                        n2 = ord_pos_int + aperture_wing;
                    }

                    bkg = 0.5 * (fluxes_min[ord][i] + fluxes_min[ord + 1][i]) / (n2 - n1 + 1);

                    if (mode == 1)
                    {
                        double sum = 0;
                        for (int j = n1; j <= n2; j++)
                        {
                            sum += im[j, i];
                        }
                        Flux[ord][i] = sum;
                    }
                    if (mode == 2)
                    {

                    }
                    FluxDebiased[ord][i] = Flux[ord][i] - bkg;
                    if (FluxDebiased[ord][i] < 0) FluxDebiased[ord][i] = 0;
                }
            }
            Console.Write("\r\n");
        }

        //private static double[] polynom(double x)
        //{
        //    double[] ans = new double[poly_degree + 1];
        //    ans[0] = 1.0;
        //    for (int i = 1; i < poly_degree + 1; i++) ans[i] = x * ans[i - 1];
        //    return ans;
        //}

        private static double ApproxAlongColumn(int column, int row, int wing, int degree)
        {
            Fitter fitter = new Fitter();
            double f1;
            double[] x = new double[1 + 2 * wing];
            double[] y = new double[1 + 2 * wing];
            double[] sigma = new double[1 + 2 * wing];

            int k1 = row - wing;
            int k2 = row + wing;

            for (int i = 0; i < x.Length; i++)
            {
                x[i] = i + k1;
                y[i] = image[i + k1, column];
                if (y[i] > 1) sigma[i] = Math.Sqrt(y[i]);
                else sigma[i] = 1;
            }

            double[] coeffs = fitter.WightedPolynom(x, y, sigma, degree);
            f1 = Polynom(coeffs, (double)row);

            return f1;
        }

        private static double Polynom(double[] coeffs, double x)
        {
            double sum = 0;
            for (int i = 0; i < coeffs.Length; i++)
            {
                sum += Math.Pow(x, i) * coeffs[i];
            }
            return sum;
        }

        public static void MinFluxes1(int[][] pos_min, int nwing, int degree)
        {
            int traces_num = pos_min.Length;
            int pixels_num = pos_min[0].Length;

            fluxes_min1 = new double[traces_num][];
            for (int i = 0; i < traces_num; i++)
            {
                fluxes_min1[i] = new double[pixels_num];
            }

            for (int i = 0; i < traces_num; i++)
            {
                for (int j = 0; j < pixels_num; j++)
                {
                    fluxes_min1[i][j] = ApproxAlongColumn(j, pos_min[i][j], 6, 3);
                }
            }
        }

        public static void MinFluxes0(int[][] pos_min)
        {
            int traces_num = pos_min.Length;
            int pixels_num = pos_min[0].Length;

            fluxes_min0 = new double[traces_num][];
            for (int i = 0; i < traces_num; i++)
            {
                fluxes_min0[i] = new double[pixels_num];
            }
            for (int i = 0; i < traces_num; i++)
            {
                for (int j = 0; j < pixels_num; j++)
                {
                    fluxes_min0[i][j] = image[pos_min[i][j], j];
                }
            }
        }

        private static void BackGroundFit(int[][] pos_min)
        {
            int traces_num = pos_min.Length;
            int pixels_num = pos_min[0].Length;

            double[] xx = new double[pixels_num];
            for (int i = 0; i < pixels_num; i++) xx[i] = (double)i / pixels_num;

            double[] yy = new double[pixels_num];

            fluxes_min = new double[traces_num][];
            for (int i = 0; i < traces_num; i++)
            {
                fluxes_min[i] = new double[pixels_num];
            }

            for (int i = 0; i < traces_num; i++)
            {
                for (int j = 0; j < pixels_num; j++)
                    yy[j] = image[pos_min[i][j], j];


                double yyMax = yy.Max();
                for (int j = 0; j < pixels_num; j++)
                    yy[j] = yy[j] / yyMax;

                RobustFitter.FitPolynom(xx, yy, 4, 1);

                for (int j = 0; j < pixels_num; j++)
                {
                    fluxes_min[i][j] = RobustFitter.FuncValue(xx[j]) * yyMax;
                }
            }
        }

        public static double[] Filter(double[] ff, int wing)
        {
            double[] filtered = new double[ff.Length];
            int xx1, xx2, nn;
            double ave;
            for (int i = 0; i < ff.Length; i++)
            {
                ave = 0;
                xx1 = i - wing;
                xx2 = i + wing;
                nn = 2 * wing + 1;
                if (xx1 < 0)
                {
                    nn = Math.Abs(xx1) + wing + 1;
                    xx1 = 0;
                }
                if (xx2 > ff.Length - 1)
                {
                    nn = Math.Abs(ff.Length - 1 - xx2) + wing + 1;
                    xx2 = ff.Length - 1;
                }
                for (int j = xx1; j <= xx2; j++)
                {
                    ave += ff[j];
                }
                ave = ave / nn;
                filtered[i] = ave;
            }
            return filtered;
        }
    }
}
