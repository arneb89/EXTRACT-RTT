using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace FITS_READER
{
    class Saver
    {
        public static void SaveOrderedDistribution(double[][] x, string file_name)
        {
            StreamWriter sw = new StreamWriter(file_name);
            for (int i = 0; i < x[0].Length; i++)
            {
                sw.Write((i + 1).ToString() + "\t");
                for (int j = 0; j < x.Length; j++)
                {
                    sw.Write(x[j][i].ToString() + "\t");
                }
                sw.Write("\r\n");
            }
            sw.Close();
        }

        public static void SaveColumn(double[] x, string file_name)
        {
            StreamWriter sw = new StreamWriter(file_name);
            for (int i = 0; i < x.Length; i++)
                sw.WriteLine("{0}", x[i]);
            sw.Close();
        }
    }
}
