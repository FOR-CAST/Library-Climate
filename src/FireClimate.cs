using System;
using System.Collections.Generic;

namespace Landis.Library.Climate
{
    // Canadian Forest Fire Weather Index System.
    // Equation numbers refer to Van Wagner, C.E. and T.L. Pickett (1985) / Van Wagner, C.E. (1987).
    // Constants follow the `cffdrs` R package (v1.9.2) so that codes calculated here are directly
    // comparable with codes calculated in R. Where cffdrs rounds a constant differently from the
    // published tables the cffdrs value is used and the difference is noted.
    public partial class Climate
    {
        #region fields

        // day-length adjustment for the Drought Code, Van Wagner (1987) Table 3
        private static readonly double[] Lf = { -1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6 };

        // effective day length for the Duff Moisture Code, Van Wagner (1987) Table 1
        private static readonly double[] Le = { 6.5, 7.5, 9.0, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8.0, 7.0, 6.0 };

        // scaling constant relating the Fine Fuel Moisture Code to moisture content, Eqs. 1 & 25.
        // 250 * 59.5 / 101 = 147.27723, which the published equations round to 147.2.
        private const double FfmcCoefficient = 250.0 * 59.5 / 101.0;

        // ClimateRecord.Precip is [cm]; the FWI equations are defined for daily rainfall in [mm]
        private const double precipTransformation = 10.0;

        #endregion

        #region private methods

        private static void CalculateDailyFireWeather(List<ClimateRecord> yearRecords)
        {
            double fineFuelMoistureCodeYesterday = ConfigParameters.FineFuelMoistureCode_Yesterday;
            double duffMoistureCodeYesterday = ConfigParameters.DuffMoistureCode_Yesterday;
            double droughtCodeYesterday = ConfigParameters.DroughtCode_Yesterday;

            // assume ConfigParameters.SpringStart and ConfigParameters.WinterStart are one-based values
            // convert to zero-based values, i.e. shift the range [1, 365] to [0, 364]
            // Therefore, a ConfigParameters.SpringStart of 30 will start on day 29 in the input log
            var start = Math.Max(0, ConfigParameters.SpringStart - 1);
            var end = Math.Min(364, ConfigParameters.WinterStart - 1);

            // loop over days from SpringStart to WinterStart
            for (var d = start; d <= end; ++d)
            {
                var record = yearRecords[d];

                // the FWI equations take daily rainfall in [mm]; ClimateRecord.Precip is [cm]
                var precip = record.Precip * precipTransformation;

                // RH is derived from dew point or specific humidity and can saturate at
                // exactly 100 percent; back it off as cffdrs does so the equilibrium
                // moisture content equations stay in range
                var rh = record.RH >= 100.0 ? 99.9999 : record.RH;

                var mo = Calculate_mo(fineFuelMoistureCodeYesterday, precip);
                var Ed = Calculate_Ed(rh, record.Temp);
                var Ew = Calculate_Ew(rh, record.Temp);
                var ko = Calculate_ko(rh, record.WindSpeed);
                var kd = Calculate_kd(ko, record.Temp);
                var kl = Calculate_kl(rh, record.WindSpeed);
                var kw = Calculate_kw(kl, record.Temp);
                var m = Calculate_m(mo, Ed, kd, Ew, kw);
                var fineFuelMoistureCode = Calculate_FineFuelMoistureCode(m);
                var re = Calculate_re(precip);
                var Mo = Calculate_Mo(duffMoistureCodeYesterday);
                var b = Calculate_b(duffMoistureCodeYesterday);
                var Mr = Calculate_Mr(re, b, Mo);
                var Pr = Calculate_Pr(Mr);
                var month = Climate.MonthOfYear(d);
                var K = Calculate_K(record.Temp, rh, Le[month]);
                var duffMoistureCode = Calculate_DuffMoistureCode(precip, Pr, K, duffMoistureCodeYesterday);
                var rd = Calculate_rd(precip);
                var Qo = Calculate_Qo(droughtCodeYesterday);
                var Qr = Calculate_Qr(Qo, rd);
                var Dr = Calculate_Dr(Qr);
                var V = Calculate_V(record.Temp, Lf[month]);
                var droughtCode = Calculate_DroughtCode(precip, Dr, V, droughtCodeYesterday);
                var windFunction_ISI = Calculate_WindFunction_ISI(record.WindSpeed);

                // Eq. 25 - the ISI is driven by the moisture content implied by today's (bounded)
                // FFMC, not by the unbounded moisture content that produced it
                var fineFuelMoistureFunction_ISI = Calculate_FineFuelMoistureFunction_ISI(Calculate_FineFuelMoisture(fineFuelMoistureCode));

                var initialSpreadIndex = Calculate_InitialSpreadIndex(windFunction_ISI, fineFuelMoistureFunction_ISI);
                var buildUpIndex = Calculate_BuildUpIndex(duffMoistureCode, droughtCode);
                var fD = Calculate_fD(buildUpIndex);
                var B = Calculate_B(initialSpreadIndex, fD);

                record.DuffMoistureCode = duffMoistureCodeYesterday = duffMoistureCode;
                record.DroughtCode = droughtCodeYesterday = droughtCode;
                record.BuildUpIndex = buildUpIndex;
                record.FineFuelMoistureCode = fineFuelMoistureCodeYesterday = fineFuelMoistureCode;
                record.FireWeatherIndex = Calculate_FireWeatherIndex(B);
            }
        }

        // Eq. 1 - moisture content [%] implied by a Fine Fuel Moisture Code
        private static double Calculate_FineFuelMoisture(double fineFuelMoistureCode) => FfmcCoefficient * (101.0 - fineFuelMoistureCode) / (59.5 + fineFuelMoistureCode);

        // Eqs. 1, 2, 3a & 3b - yesterday's moisture content after the rainfall phase
        private static double Calculate_mo(double fineFuelMoistureCode, double precip)
        {
            var mo = Calculate_FineFuelMoisture(fineFuelMoistureCode);

            if (precip > 0.5)
            {
                // Eq. 2 - rainfall reaching the fuel bed, reduced by canopy interception
                var rf = precip - 0.5;

                // Eq. 3a
                var wetting = 42.5 * rf * Math.Exp(-100.0 / (251.0 - mo)) * (1.0 - Math.Exp(-6.93 / rf));

                // Eq. 3b - additional wetting when the fuel is already very wet
                if (mo > 150.0)
                    wetting += 0.0015 * (mo - 150.0) * (mo - 150.0) * Math.Sqrt(rf);

                mo += wetting;
            }

            // the moisture content of pine litter tops out at about 250 percent
            return Math.Min(250.0, mo);
        }

        private static double Calculate_Ed(double rh, double temp) => 0.942 * Math.Pow(rh, 0.679) + 11.0 * Math.Exp((rh - 100.0) / 10.0) + 0.18 * (21.1 - temp) * (1.0 - Math.Exp(-0.115 * rh));

        private static double Calculate_Ew(double rh, double temp) => 0.618 * Math.Pow(rh, 0.753) + 10.0 * Math.Exp((rh - 100.0) / 10.0) + 0.18 * (21.1 - temp) * (1.0 - Math.Exp(-0.115 * rh));  //selfs

        private static double Calculate_ko(double rh, double windSpeed) => 0.424 * (1.0 - Math.Pow(rh / 100.0, 1.7)) + 0.0694 * Math.Pow(windSpeed, 0.5) * (1.0 - Math.Pow((rh / 100.0), 8));

        private static double Calculate_kd(double ko, double temp) => ko * 0.581 * Math.Exp(0.0365 * temp);

        private static double Calculate_kl(double rh, double windSpeed) => 0.424 * (1.0 - Math.Pow((100.0 - rh) / 100.0, 1.7)) + 0.0694 * Math.Pow(windSpeed, 0.5) * (1.0 - Math.Pow(((100.0 - rh) / 100.0), 8));

        private static double Calculate_kw(double kl, double temp) => kl * 0.581 * Math.Exp(0.0365 * temp);

        private static double Calculate_m(double mo, double Ed, double kd, double Ew, double kw)
        {
            if (mo > Ed)
                return Ed + (mo - Ed) * Math.Pow(10.0, -kd);

            return (mo < Ed && mo < Ew) ? Ew - (Ew - mo) * Math.Pow(10.0, -kw) : mo;
        }

        // Eq. 10 - the FFMC scale runs from 0 to 101
        private static double Calculate_FineFuelMoistureCode(double m) => Math.Min(101.0, Math.Max(0.0, 59.5 * (250.0 - m) / (FfmcCoefficient + m)));

        // Eq. 11
        private static double Calculate_re(double precip) => precip > 1.5 ? 0.92 * precip - 1.27 : 0.0;

        // Eq. 12. cffdrs writes this as 20 + 280 / exp(0.023 * DMC); the published form is
        // 20 + exp(5.6348 - DMC / 43.43), which differs by about 0.1 percent in Mo.
        private static double Calculate_Mo(double duffMoistureCode) => 20.0 + 280.0 / Math.Exp(0.023 * duffMoistureCode);

        // Eqs. 13a, 13b & 13c
        private static double Calculate_b(double duffMoistureCode)
        {
            if (duffMoistureCode <= 33.0)
            { 
                return 100.0 / (0.5 + 0.3 * duffMoistureCode);
            }

            if (duffMoistureCode > 65.0)
            {
                return 6.2 * Math.Log(duffMoistureCode) - 17.2;
            }

            return 14.0 - 1.3 * Math.Log(duffMoistureCode);
        }

        // Eq. 14
        private static double Calculate_Mr(double re, double b, double Mo) => Mo + 1000.0 * re / (48.77 + b * re);

        // Eq. 15
        private static double Calculate_Pr(double Mr) => Math.Max(0.0, 43.43 * (5.6348 - Math.Log(Mr - 20.0)));

        // Eq. 16 - the log drying rate, scaled by 1e-6 because Eq. 17 below multiplies by 100
        private static double Calculate_K(double temp, double rh, double Le) => temp < -1.1 ? 0.0 : 1.894 * (temp + 1.1) * (100.0 - rh) * Le * Math.Pow(10.0, -6.0);

        // Eq. 17 - the DMC scale has no upper bound but cannot go below zero
        private static double Calculate_DuffMoistureCode(double precip, double Pr, double K, double duffMoistureCodeYesterday) => Math.Max(0.0, (precip > 1.5 ? Pr : duffMoistureCodeYesterday) + 100.0 * K);

        // Eq. 19
        private static double Calculate_rd(double precip) => precip > 2.8 ? 0.83 * precip - 1.27 : 0.0;

        // Eq. 20
        private static double Calculate_Qo(double droughtCode) => 800.0 * Math.Exp(-droughtCode / 400.0);

        // Eq. 21
        private static double Calculate_Qr(double Qo, double rd) => Qo + 3.937 * rd;

        // Eq. 22
        private static double Calculate_Dr(double Qr) => Math.Max(0.0, 400.0 * Math.Log(800.0 / Qr));

        // Eq. 18 - twice the day-length adjusted drying factor, which cannot be negative.
        // Without the lower bound the winter values of Lf drive the Drought Code below zero.
        private static double Calculate_V(double temp, double Lf) => Math.Max(0.0, temp < -2.8 ? Lf : 0.36 * (temp + 2.8) + Lf);

        // Eq. 23 - the DC scale has no upper bound but cannot go below zero
        private static double Calculate_DroughtCode(double precip, double Dr, double V, double droughtCodeYesterday) => Math.Max(0.0, (precip > 2.8 ? Dr : droughtCodeYesterday) + 0.5 * V);

        // Eq. 24
        private static double Calculate_WindFunction_ISI(double windSpeed) => Math.Exp(0.05039 * windSpeed);

        // Eq. 26
        private static double Calculate_FineFuelMoistureFunction_ISI(double m) => 91.9 * Math.Exp(-0.1386 * m) * (1.0 + Math.Pow(m, 5.31) / (4.93 * Math.Pow(10.0, 7.0)));

        // Eq. 27
        private static double Calculate_InitialSpreadIndex(double windFunction_ISI, double fineFuelMoistureFunction_ISI) => 0.208 * windFunction_ISI * fineFuelMoistureFunction_ISI;

        // Eqs. 27a & 27b
        private static double Calculate_BuildUpIndex(double duffMoistureCode, double droughtCode)
        {
            if (duffMoistureCode <= 0.0 && droughtCode <= 0.0)
                return 0.0;

            // Eq. 27a
            var buildUpIndex = 0.8 * duffMoistureCode * droughtCode / (duffMoistureCode + 0.4 * droughtCode);

            if (buildUpIndex >= duffMoistureCode)
                return buildUpIndex;

            // Eq. 27b. The correction term is (0.0114 * DMC) ^ 1.7, not 0.0114 * DMC ^ 1.7.
            var p = duffMoistureCode <= 0.0 ? 0.0 : (duffMoistureCode - buildUpIndex) / duffMoistureCode;

            return Math.Max(0.0, duffMoistureCode - p * (0.92 + Math.Pow(0.0114 * duffMoistureCode, 1.7)));
        }

        // Eq. 28
        private static double Calculate_fD(double buildUpIndex) => buildUpIndex <= 80.0 ? 0.626 * Math.Pow(buildUpIndex, 0.809) + 2.0 : 1000.0 / (25.0 + 108.64 * Math.Exp(-0.023 * buildUpIndex));

        // Eq. 29
        private static double Calculate_B(double initialSpreadIndex, double fD) => 0.1 * initialSpreadIndex * fD;

        // Eqs. 30a & 30b
        private static double Calculate_FireWeatherIndex(double B) => B > 1.0 ? Math.Exp(2.72 * Math.Pow(0.434 * Math.Log(B), 0.647)) : B; 

        #endregion
    }
}
