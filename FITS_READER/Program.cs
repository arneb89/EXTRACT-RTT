using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace FITS_READER
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.Beep();
            Console.ForegroundColor = ConsoleColor.Yellow;
            Console.WriteLine("*******************************************************");
            Console.WriteLine("*                    RTT-150 CES                      *");
            Console.WriteLine("*            SPECTRA REDUCTION PIPELINE               *");
            Console.WriteLine("*******************************************************");
            Console.ForegroundColor = ConsoleColor.White;

            //InitFile.ReadInitFile("INIT.dat");
            Init.Initialize("INIT.dat");
            if (Init.ErrorString != "")
            {
                Console.WriteLine("Error in INIT file!");
                Console.WriteLine(Init.ErrorString);
                goto STOP;
            }

            string dir_main = (string)Init.Value("DIR_MAIN");

            if (!Directory.Exists(dir_main))
            {
                try
                {
                    Directory.CreateDirectory(dir_main);
                }
                catch
                {
                    Console.WriteLine("Cannot create main directory.\r\n");
                    goto STOP;
                }
            }

            string[] bias_files = null, flat_files = null, thar_files = null, obj_files = null;
            Console.WriteLine("Searching for Bias, Flat, ThAr and Object files...");
            try
            {
                bias_files = ImageSelector.GetFilesByHDU(
                    (string)Init.Value("DIR_BIAS"), (string)Init.Value("FMASK_BIAS"), "imagetyp", "*");
                flat_files = ImageSelector.GetFilesByHDU(
                    (string)Init.Value("DIR_FLAT"), (string)Init.Value("FMASK_FLAT"), "imagetyp", "*");
                thar_files = ImageSelector.GetFilesByHDU(
                    (string)Init.Value("DIR_CLBR"), (string)Init.Value("FMASK_CLBR"), "imagetyp", "*");
                obj_files = ImageSelector.GetFilesByHDU(
                    (string)Init.Value("DIR_OBJ"), (string)Init.Value("FMASK_OBJ"), "imagetyp", "*");
            }
            catch
            {
                Console.WriteLine("Error in files searching...");
                goto STOP;
            }

            Console.WriteLine("Bias files number: {0}", bias_files.Length);
            Console.WriteLine("Flat files number: {0}", flat_files.Length);
            Console.WriteLine("ThAr files number: {0}", thar_files.Length);
            Console.WriteLine("Object files number: {0}", obj_files.Length);

            if (bias_files.Length == 0 || flat_files.Length == 0 ||
                thar_files.Length == 0 || obj_files.Length == 0)
            {
                Console.WriteLine("Some necessary files were not found...");
                goto STOP;
            }

            Image[] biases = new Image[bias_files.Length];
            Image[] flats = new Image[flat_files.Length];
            Image[] thars = new Image[thar_files.Length];
            Image[] objects = new Image[obj_files.Length];

            Console.WriteLine("Loading Bias images...");
            for (int i = 0; i < bias_files.Length; i++)
            {
                Console.WriteLine("->" + bias_files[i]);
                biases[i] = new Image();
                biases[i].LoadImage(bias_files[i]);
            }
            Console.WriteLine("Bias images averaging...");
            Image bias_aver = ImageCombinator.Median(biases);
            for (int i = 0; i < bias_files.Length; i++) biases[i] = null;

            Console.WriteLine("Loading Flat images...");
            for (int i = 0; i < flat_files.Length; i++)
            {
                Console.WriteLine("->" + flat_files[i]);
                flats[i]=new Image();
                flats[i].LoadImage(flat_files[i]);
            }
            Console.WriteLine("Flat images averaging...");
            Image flat_aver = ImageCombinator.Median(flats);
            for (int i = 0; i < flat_files.Length; i++) flats[i] = null;

            Console.WriteLine("Loading Th-Ar images...");
            for (int i = 0; i < thar_files.Length; i++)
            {
                Console.WriteLine("->" + thar_files[i]);
                thars[i] = new Image();
                thars[i].LoadImage(thar_files[i]); 
            }
            Console.WriteLine("ThAr images averaging...");
            Image Thar_aver = ImageCombinator.Median(thars);
            for (int i = 0; i < thar_files.Length; i++) thars[i] = null;
           
            Console.WriteLine("Bias substraction from averanged Flat image...");
            Image flat_aver_db = flat_aver - bias_aver;

            Console.WriteLine("Bias substraction from averanged ThAr image...");
            Image Thar_aver_db = Thar_aver - bias_aver;

            Console.WriteLine("Replasing pixels values:\r\n negative -> 1 for flat;\r\n negative -> 0 for Th-Ar;");
            for(int i=0; i<flat_aver_db.NAXIS1; i++)
                for (int j = 0; j < flat_aver_db.NAXIS2; j++)
                {
                    if (flat_aver_db[i, j] <= 0) flat_aver_db[i, j] = 1;
                    if (Thar_aver_db[i, j] < 0) Thar_aver_db[i, j] = 0;
                }

            Console.WriteLine("Search for order locations...");
            int polinim_degree_order = 3;
            Locator.Locate(ref flat_aver_db, polinim_degree_order);

            double[][] pos_ord = Locator.Ord_Pos;
            double[][] pos_min = Locator.Min_Pos;

            double[][] fluxes;

            Console.WriteLine("Extraction of the flat spectra...");
            Extractor.Extract(flat_aver_db, pos_ord, pos_min, -1);
            fluxes = Extractor.FluxDebiased;
            Saver.SaveOrderedDistribution(fluxes, dir_main + "\\Flat_Orders.dat");

            // Flat spectra normalization;
            Console.ForegroundColor = ConsoleColor.Yellow;
            Console.WriteLine("Flat spectra normalization...");
            Console.ForegroundColor = ConsoleColor.White;

            Normator.Norm(fluxes, 12, 0, 13);
            Saver.SaveOrderedDistribution(Normator.OrdersNorm, dir_main + "\\Flat_Orders_Norm.dat");
            Saver.SaveOrderedDistribution(Normator.FitCurves, dir_main + "\\Flat_Orders_Fit.dat");

            // Th-Ar spectra extraction;
            Console.ForegroundColor = ConsoleColor.Yellow;
            Console.WriteLine("Extraction of the Th-Ar spectrum...");
            Console.ForegroundColor = ConsoleColor.White;

            Extractor.Extract(Thar_aver_db, pos_ord, pos_min,-1);

            fluxes = Extractor.FluxDebiased;

            Saver.SaveOrderedDistribution(fluxes, dir_main + "\\thar_extacted.txt");
            Saver.SaveOrderedDistribution(Extractor.Background, dir_main + "\\thar_bkg.dat");
            Saver.SaveOrderedDistribution(Extractor.Flux, dir_main + "\\thar_flx_and_bkg.dat");
            Saver.SaveOrderedDistribution(Extractor.fluxes_min0, dir_main + "\\Back.0.dat");

            // Wavelength calibration;
            double[][] lambdas;
            Console.ForegroundColor = ConsoleColor.Yellow;
            Console.WriteLine("Wavelength calibration...");
            Console.ForegroundColor = ConsoleColor.White;
            
            EcheData.LoadDispCurves("CalibrationCurves.dat");

            //int oo = 9, ox = 5;
            //double cutLimit = 3, fluxLimit = 500;
            //int iterNum = 20;

            int oo = (int)Init.Value("WAVE_OO");
            int ox = (int)Init.Value("WAVE_OX");
            int iterNum = (int)Init.Value("WAVE_NINER");
            double cutLimit = (double)Init.Value("WAVE_REJ");
            double fluxLimit = 500;

            Console.WriteLine("OO = {0}; OX = {1}; Reject = {2}; IterNum = {3}; FluxLimit = {4}", 
                oo, ox, cutLimit, iterNum, fluxLimit);
            WLCalibration.Calibrate(fluxes, oo, ox, cutLimit, iterNum, fluxLimit);
            lambdas = WLCalibration.Lambdas;
            Saver.SaveOrderedDistribution(lambdas, dir_main + "\\lambdas.dat");

            // processing object images;
            Console.ForegroundColor = ConsoleColor.Yellow;
            Console.WriteLine("Processing object spectra...");
            Console.ForegroundColor = ConsoleColor.White;

            // observatory coordinates and altitude reading;
            IIDType obs_lat = (IIDType)Init.Value("OBS_LAT");
            IIDType obs_long = (IIDType)Init.Value("OBS_LONG");
            double obs_alt = (double)Init.Value("OBS_ALT");
            double obs_lat_deg = obs_lat.DH + obs_lat.MM / 60.0 + obs_lat.SS / 3600.0;
            double obs_long_deg = obs_long.DH + obs_long.MM / 60.0 + obs_long.SS / 3600.0;
            int time_zone = (int)Init.Value("TIME_ZONE");
            string time_type = (string)Init.Value("TIME_TYPE");

            string ra_str = (string)Init.Value("OBJ_RASTR");
            string de_str = (string)Init.Value("OBJ_DESTR");
            bool read_coord_fits = (bool)Init.Value("OBJ_INFITS");

            IIDType ra = null, de = null;
            double ra_hh = 0, de_deg = 0;
            string time_str = (string)Init.Value("TIME_STR");
            string date_str = (string)Init.Value("DATE_STR");
            string expo_str = (string)Init.Value("EXPO_STR");

            if (!read_coord_fits)
            {
                ra = (IIDType)Init.Value("OBJ_RA");
                de = (IIDType)Init.Value("OBJ_DEC");
                ra_hh = ra.DH + ra.MM / 60.0 + ra.SS / 3600.0;
                de_deg = ra.DH + ra.MM / 60.0 + ra.SS / 3600.0;
            }

            StreamWriter output = new StreamWriter(dir_main + "\\output.dat");

            FITSHeaderReader fhr = new FITSHeaderReader();
            Image object_image = null;
            for (int i = 0; i < obj_files.Length; i++)
            {
                // Lading object image;
                Console.WriteLine("Loading Object image {0}", obj_files[i]);
                object_image = new Image();
                object_image.LoadImage(obj_files[i]);

                // Reading FITS-header data;
                if (read_coord_fits)
                {
                    ra = fhr.ReadAsIID(obj_files[i], ra_str);
                    de = fhr.ReadAsIID(obj_files[i], de_str);
                    ra_hh = ra.DH + ra.MM / 60.0 + ra.SS / 3600.0;
                    de_deg = ra.DH + ra.MM / 60.0 + ra.SS / 3600.0;
                }

                IIDType time = fhr.ReadAsIID(obj_files[i], time_str);
                IIDType date = fhr.ReadAsIID(obj_files[i], date_str);
                double exposure = fhr.ReadAsDouble(obj_files[i], expo_str);

                DateTime dt = new DateTime(date.DH, date.MM, (int)date.SS,
                    time.DH, time.MM, (int)time.SS);
                if (time_type == "LT") dt.AddHours(-(double)time_zone);
                dt.AddSeconds(0.5 * exposure);

                // Calculate JD and HJD;
                double jd = DateConvertors.JDCnv(dt.Year, dt.Month, dt.Day,
                    (double)dt.Hour + dt.Minute / 60.0 + dt.Second / 3600.0);
                double epoch = dt.Year + dt.Month / 12.0 + dt.Day / 365.0 +
                    dt.Hour / (24 * 365.0) + dt.Minute / (24 * 365.0 * 60.0);
                Precess.Correct(ref ra_hh, ref de_deg, 2000.0, epoch);
                double hjd = DateConvertors.helio_jd(jd, ra_hh, de_deg);

                // Calculate components of observer motion to the object direction;
                double vdiurnal = 0, vbar = 0, vhel = 0, corr = 0;
                HelCorr.DoCorrection(obs_long_deg, obs_lat_deg, obs_alt, 360-ra_hh, de_deg, jd, 
                    ref vdiurnal, ref vbar, ref vhel);
                corr = vdiurnal + vbar + vhel;
                
                // Bias substraction;
                Console.WriteLine("Bias substraction...");
                object_image = object_image - bias_aver;
                for (int j = 0; j < object_image.NAXIS1; j++)
                    for (int k = 0; k < object_image.NAXIS2; k++)
                        if (object_image[j, k] < 0) object_image[j, k] = 0;

                // Image optimization;
                Smooth.DoSmooth(ref object_image);

                // Create directory for extracted spectra;
                string dirName = obj_files[i].Substring(0, obj_files[i].IndexOf(".fit"));
                int first = obj_files[i].LastIndexOf("\\");
                int last = obj_files[i].IndexOf(".fit");
                string fileName = dirName + "\\" + obj_files[i].Substring(first, last - first);
                Directory.CreateDirectory(dirName);

                // save object data;
                output.WriteLine("{0}\t{1}\t{2}\t{3}"+
                    "\t{4:0.0000000}\t{5:0.0000000}\t{6:0.000}\t{7:0.000}\t{8:0.000}\t{9:0.000}",
                    obj_files[i].Substring(first + 1), ra.ToString(), de.ToString(),
                    dt.ToString(), jd, hjd, vdiurnal, vbar, vhel, corr);
                output.Flush();

                Console.WriteLine("Extraction spectra from {0} ", obj_files[i]);
                Extractor.Extract(object_image, pos_ord, pos_min, -1);
                fluxes = Extractor.FluxDebiased;
                double[][] backgr = Extractor.Background;
                Saver.SaveOrderedDistribution(backgr, fileName + "_bkg.dat");
                Saver.SaveOrderedDistribution(fluxes, fileName + "_x.dat");
                for (int k = 0; k < fluxes.Length; k++)
                {
                    for (int j = 0; j < fluxes[0].Length; j++)
                    {
                        fluxes[k][j] = fluxes[k][j] / Normator.OrdersNorm[k][j];
                    }
                }
                Saver.SaveOrderedDistribution(fluxes, fileName + "_x_n.dat");
                

                for (int k = 0; k < fluxes.Length; k++)
                {
                    StreamWriter sw = new StreamWriter(fileName + 
                        string.Format("_{0:000}", k) + ".dat");
                    for (int j = 0; j < fluxes[k].Length; j++)
                    {
                        sw.WriteLine(
                            string.Format("{0:0000.0000}\t{1:0.00000E000}", 
                            lambdas[k][j], fluxes[k][j]).Replace(",", "."));
                    }
                    sw.Close();
                }
            }

            output.Close();

        STOP:
            Console.ForegroundColor = ConsoleColor.Yellow;
            Console.WriteLine("End of the pipeline. Press any key to exit...");
            Console.Beep();
            Console.ReadKey();
        }
    }
}
