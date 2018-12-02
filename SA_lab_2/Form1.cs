using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Windows.Forms;
using System.IO;
using ZedGraph;

namespace SA_lab_2
{

    public partial class Form1 : Form
    {
        const int N = 3;
        int n, dy;
        double[,] vlnovFunc;
        const double e = 0.00000001;
        bool Approximate = false, Flag = true, Vlasn = false;
        Stream XStream, XForecast;
        Stream YStream;
        Stream OStream;
        Graphics gPanel;
        double[][,] X = new double[N][,];
        double[][,] Z;
        double[][] Y, MinX, MaxX, ZY, logZY;
        double[] MinY, MaxY, diff;
        int[] d;
        double[][][][] PolyCoef;
        double[][][] ksi = new double[N][][];
        public Form1()
        {
            InitializeComponent();
            openFileDialog2.Multiselect = true;
            gPanel = panel1.CreateGraphics();

        }
        public class Polynom
        {
            double[] A;
            public Polynom(double[] a)
            {
                A = (double[])a.Clone();
            }
            public double Method(double x)
            {
                double res = 0;
                for (int i = 0; i < A.Length; i++)
                {
                    double t = A[i];
                    for (int j = 0; j < i; j++) t *= x;
                    res += t;
                }
                return res;
            }
        }
        double Chebishev(int P, double x)
        {
            if (P == 0) return 1;
            if (P == 1) return x;
            return 2 * x * Chebishev(P - 1, x) - Chebishev(P - 2, x);
        }
        double SChebishev(int P, double x)
        {
            if (P == 0) return 0.5;
            return Chebishev(P, 2 * x - 1);
        }
        double SChebishev2(int P, double x)
        {
            if (P == 0) return 0.5;
            return Chebishev2(P, x);
        }
        double Chebishev2(int P, double x)
        {
            if (P == 0) return 1;
            if (P == 1) return 2 * x;
            return 2 * x * Chebishev2(P - 1, x) - Chebishev2(P - 2, x);
        }
        double SChebishev3(int P, double x)
        {
            if (P == 0) return 0.5;
            return Chebishev3(P, x);
        }
        double Chebishev3(int P, double x)
        {
            if (P == 0) return 1;
            if (P == 1) return x;
            return x * Chebishev3(P - 1, x) - Chebishev3(P - 2, x);
        }
        double[][] PolynomChebishev(int p)
        {
            double[][] T = new double[p + 1][];
            for (int i = 0; i < p + 1; i++)
            {
                T[i] = new double[p + 1];
            }

            T[0][0] = 1;
            T[1][1] = 0.5;
            for (int i = 2; i < p + 1; i++)
            {
                for (int j = 0; j < p + 1; j++)
                {
                    T[i][j] = -T[i - 2][j];
                    if (j > 0)
                        T[i][j] += T[i - 1][j - 1];
                }
            }
            return T;
        }
        double[][] PolynomChebishev2(int p)
        {
            double[][] T = new double[p + 1][];
            for (int i = 0; i < p + 1; i++)
            {
                T[i] = new double[p + 1];
            }

            T[0][0] = 1;
            T[1][1] = 2;
            for (int i = 2; i < p + 1; i++)
            {
                for (int j = 0; j < p + 1; j++)
                {
                    T[i][j] = -T[i - 2][j];
                    if (j > 0)
                        T[i][j] += 2 * T[i - 1][j - 1];
                }
            }
            return T;
        }
        double[][] PolynomChebishev3(int p)
        {
            double[][] T = new double[p + 1][];
            for (int i = 0; i < p + 1; i++)
            {
                T[i] = new double[p + 1];
            }

            T[0][0] = 1;
            T[1][1] = 1;
            for (int i = 2; i < p + 1; i++)
            {
                for (int j = 0; j < p + 1; j++)
                {
                    T[i][j] = -T[i - 2][j];
                    if (j > 0)
                        T[i][j] += T[i - 1][j - 1];
                }
            }

            return T;
        }
        public delegate double funk(int P, double x);
        public delegate double[][] Pfunk(int P);
        private double[] minus(double[] a, double[] b)
        {
            double[] t = new double[a.Length];
            for (int i = 0; i < a.Length; i++)
            {
                t[i] = a[i] - b[i];
            }
            return t;
        }
        private double norm(double[] a)
        {
            double t = 0;
            for (int i = 0; i < a.Length; i++)
            {
                t += a[i] * a[i];
            }
            t = Math.Sqrt(t);
            return t;
        }
        private void xload_Click(object sender, EventArgs e)
        {
            if (openFileDialog1.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    if ((XStream = openFileDialog1.OpenFile()) != null)
                    {
                        Xinput.Text = openFileDialog1.SafeFileName;
                    }
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Unable to read a file. " + ex.Message);
                }
            }
        }
        double[] SpragenGradKvadrFunction(double[] x, double[][] A, double[] b)
        {
            double[][] At = new double[A[0].Length][];
            for (int i = 0; i < At.Length; i++)
                At[i] = new double[A.Length];
            At = Transpose(A);
            b = MultiMatrVect(At, b);
            A = MultMatr(At, A);
            double[] r = new double[b.Length];
            double[] x1 = new double[x.Length];
            double[] r1 = new double[b.Length];
            double[] z = new double[b.Length];
            double[] temp = new double[A.Length];
            r1 = minus(b, MultiMatrVect(A, x));
            for (int i = 0; i < r.Length; i++)
                z[i] = r1[i];
            double alfa, bet;
            int count = 0;
            do
            {
                count++;
                temp = MultiMatrVect(A, z);
                alfa = Skalar(r1, r1) / Skalar(temp, z);
                for (int i = 0; i < x.Length; i++)
                {
                    x[i] = x1[i];
                    r[i] = r1[i];
                    x1[i] = x[i] + alfa * z[i];
                    r1[i] = r[i] - alfa * temp[i];
                }
                bet = Skalar(r1, r1) / Skalar(r, r);
                for (int i = 0; i < x.Length; i++)
                    z[i] = r1[i] + bet * z[i];
                double t = MaxVec(minus(MultiMatrVect(A, x1), b));
            }
            while ((MaxVec(minus(MultiMatrVect(A, x1), b)) > e) && (norm(minus(x, x1)) > Math.Pow(e, 2)));
            return x;
        }
        public double MaxVec(double[] x)
        {
            double max = Math.Abs(x[0]);
            for (int i = 1; i < x.Length; i++)
                if (max < Math.Abs(x[i]))
                    max = Math.Abs(x[i]);
            return max;
        }
        public double[] MultiMatrVect(double[][] A, double[] x)
        {
            double[] mult = new double[A.Length];
            for (int i = 0; i < A.Length; i++)
            {
                for (int j = 0; j < x.Length; j++)
                    mult[i] += A[i][j] * x[j];
            }
            return mult;
        }

        public double Skalar(double[] x, double[] y)
        {
            double skal = 0;
            for (int i = 0; i < x.Length; i++)
                skal += x[i] * y[i];
            return skal;
        }

        public double[][] Transpose(double[][] A)
        {
            double[][] Atr = new double[A[0].Length][];
            for (int i = 0; i < A[0].Length; i++)
                Atr[i] = new double[A.Length];

            for (int i = 0; i < A.Length; i++)
                for (int j = 0; j < Atr.Length; j++)
                    Atr[j][i] = A[i][j];

            return Atr;
        }
        double[][] MultMatr(double[][] At, double[][] A)
        {
            double[][] mult = new double[At.Length][];
            for (int i = 0; i < At.Length; i++)
            {
                mult[i] = new double[A[0].Length];
                for (int j = 0; j < A[0].Length; j++)
                    for (int k = 0; k < A.Length; k++)
                        mult[i][j] += At[i][k] * A[k][j];
            }
            return mult;
        }
        private void outfilebutton_Click(object sender, EventArgs e)
        {
            saveFileDialog1.Filter = "Text files (*.txt)|*.txt";
            if (saveFileDialog1.ShowDialog() == DialogResult.OK)
            {
                if ((OStream = saveFileDialog1.OpenFile()) != null)
                {
                    outfile.Text = saveFileDialog1.FileName;
                    OStream.Close();
                }
            }
        }
        private void addbutton_Click(object sender, EventArgs e)
        {

            if (openFileDialog2.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    if ((YStream = openFileDialog2.OpenFile()) != null)
                    {
                        Yinput.Text = openFileDialog2.SafeFileName;
                    }
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Unable to read a file. " + ex.Message);
                }
            }
        }
        private void ReadData()
        {
            if (PolinoType.SelectedIndex < 0)
                PolinoType.SelectedIndex = 0;

            textBox1.Text = null;
            // Error Checking
            if (Xinput.Text == null || outfile.Text == null || Yinput.Text == null)
            {
                MessageBox.Show("Enter input files");
                return;
            }
            OStream = saveFileDialog1.OpenFile();
            // Initializing


            n = Convert.ToInt32(Range.Value);
            d = new int[N];
            d[0] = Convert.ToInt32(dim1.Value);
            X[0] = new double[n, d[0]];
            d[1] = Convert.ToInt32(dim2.Value);
            X[1] = new double[n, d[1]];
            d[2] = Convert.ToInt32(dim3.Value);
            X[2] = new double[n, d[2]];
            dy = Convert.ToInt32(dimy.Value);
            this.numericUpDown1.Maximum = new decimal(new int[] {
                dy,
                0,
                0,
                0});
            Y = new double[dy][];
            for (int i = 0; i < dy; i++)
                Y[i] = new double[n];
            PolyCoef = new double[dy][][][];
            for (int m = 0; m < dy; m++)
                PolyCoef[m] = new double[N][][];

            // Reading
            StreamReader Xr = new StreamReader(XStream);
            StreamReader Yr = new StreamReader(YStream);
            String t;
            double[] arr;
            try
            {
                t = Xr.ReadToEnd();
                arr = t.Split(' ', '\t', '\r', '\n', '	').Where(l => !string.IsNullOrEmpty(l)).Select(q => double.Parse(q)).ToArray();
            }
            catch
            {
                MessageBox.Show("Incorrect format of input files");
                return;
            }
            for (int i = 0, k = 0; i < n; i++)
            {
                k++;
                for (int l = 0; l < N; l++)
                    for (int j = 0; j < d[l]; j++, k++)
                        X[l][i, j] = arr[k];
            }
            try
            {
                t = Yr.ReadToEnd();
                arr = t.Split(' ', '\t', '\r', '\n', ' ').Where(l => !string.IsNullOrEmpty(l)).Select(q => double.Parse(q)).ToArray();
            }
            catch
            {
                MessageBox.Show("Incorrect format of files Y ");
                return;
            }
            int p = 0;
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < dy; i++, p++)
                {
                    Y[i][j] = arr[p];
                }
            }
        }
        private void MinMaxFound()
        {
            Result.Text += "Max and min values of vectors:\r\n";
            Result.Refresh();
            MinX = new double[N][];
            for (int i = 0; i < N; i++)
                MinX[i] = new double[d[i]];
            MaxX = new double[N][];
            for (int i = 0; i < N; i++)
                MaxX[i] = new double[d[i]];
            MinY = new double[dy];
            MaxY = new double[dy];

            for (int i = 0; i < N; i++)
                for (int j = 0; j < d[i]; j++)
                    MinX[i][j] = MaxX[i][j] = X[i][0, j];
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < d[i]; j++)
                {
                    for (int q0 = 0; q0 < n; q0++)
                    {
                        if (MinX[i][j] > X[i][q0, j]) MinX[i][j] = X[i][q0, j];
                        if (MaxX[i][j] < X[i][q0, j]) MaxX[i][j] = X[i][q0, j];
                    }
                }
            }

            for (int i = 0; i < dy; i++)
                MinY[i] = MaxY[i] = Y[i][0];
            for (int i = 0; i < dy; i++)
            {
                foreach (double x in Y[i]) if (MinY[i] > x) MinY[i] = x;
                foreach (double x in Y[i]) if (MaxY[i] < x) MaxY[i] = x;
            }
            for (int i = 0; i < dy; i++)
                Result.Text += "min Y" + (i + 1) + " = " + MinY[i] + "\tmax Y" + (i + 1) + " = " + MaxY[i] + "\r\n";
        }
        private void Normalizing()
        {
            Result.Refresh();
            Z = new double[N][,];
            ZY = new double[dy][];
            logZY = new double[dy][];
            for (int i = 0; i < N; i++)
            {
                Z[i] = new double[n, d[i]];
                for (int k = 0; k < n; k++)
                    for (int w = 0; w < d[i]; w++)
                        Z[i][k, w] = (X[i][k, w] - MinX[i][w]) / (MaxX[i][w] - MinX[i][w]);

            }
            for (int i = 0; i < dy; i++)
            {
                ZY[i] = new double[n];
                logZY[i] = new double[ZY[i].Length];
                for (int k = 0; k < n; k++)
                {
                    ZY[i][k] = (Y[i][k] - MinY[i]) / (MaxY[i] - MinY[i]);
                    logZY[i][k] = ZY[i][k] + 1;
                    if (logZY[i][k] == 0)
                        logZY[i][k] += e;
                    logZY[i][k] = Math.Log(logZY[i][k]);
                }
            }
        }
        private double[] Bq()
        {
            Result.Refresh();
            Result.Text += "\r\n";
            double[] B = new double[n];
            double Maxy, Miny;

            if (Bq_arifmet.Checked)
            {
                for (int q0 = 0; q0 < n; q0++)
                {
                    Miny = ZY[0][q0];
                    Maxy = ZY[0][q0];
                    for (int i = 1; i < dy; i++)
                    {
                        if (ZY[i][q0] > Maxy) Maxy = ZY[i][q0];
                        if (ZY[i][q0] < Miny) Miny = ZY[i][q0];
                    }
                    //B[q0] = (Maxy + Miny) / 2.0;
                    B[q0] = (Maxy + Miny) / 2.0 + 1;
                    if (B[q0] == 0)
                        B[q0] += e;
                    B[q0] = Math.Log(B[q0]);
                }
            }
            else
            {
                for (int q0 = 0; q0 < n; q0++)
                {
                    Miny = ZY[0][q0];
                    Maxy = ZY[0][q0];
                    for (int i = 1; i < dy; i++)
                    {
                        if (ZY[i][q0] > Maxy) Maxy = ZY[i][q0];
                        if (ZY[i][q0] < Miny) Miny = ZY[i][q0];
                    }
                    //B[q0] = Maxy - Miny;
                    B[q0] = Math.Log(Maxy - Miny + 1);
                }
            }
            Result.Update();
            return B;
        }
        private double[][,] LambdaSearch(double[] B, int[] P, double[][,] lamb)
        {
            funk f;
            Result.Refresh();
            f = new funk(SChebishev2);
            if (PolinoType.SelectedIndex == 2)
                f = new funk(SChebishev3);
            if (PolinoType.SelectedIndex == 0)
            {
                f = new funk(SChebishev);
                if (!checkBox1.Checked)
                {
                    double[] megalamb = new double[P[0] * d[0] + P[1] * d[1] + P[2] * d[2]];
                    double[][] T = new double[n][];
                    int k = 0;
                    for (int q0 = 0; q0 < n; q0++)
                    {
                        T[q0] = new double[P[0] * d[0] + P[1] * d[1] + P[2] * d[2]];
                        k = 0;
                        for (int i = 0; i < N; i++)
                        {
                            ksi[i] = new double[d[i]][];
                            for (int j = 0; j < d[i]; j++)
                            {
                                ksi[i][j] = new double[n];
                                for (int p = 0; p < P[i]; p++)
                                {
                                    if (Vlasn)
                                        T[q0][k] = Math.Cosh(f(p + 1, (2 + Z[i][q0, j]) / 4));
                                    else
                                        T[q0][k] = 1 + f(p + 1, (2 + Z[i][q0, j]) / 4);
                                    if (T[q0][k] <= 0)
                                        T[q0][k] = e;
                                    T[q0][k] = Math.Log(T[q0][k]);
                                    k++;
                                }
                            }
                        }
                    }

                    megalamb = SpragenGradKvadrFunction(megalamb, T, B);

                    k = 0;
                    for (int i = 0; i < N; i++)
                        for (int j = 0; j < d[i]; j++)
                            for (int p = 1; p <= P[i]; p++, k++)
                                lamb[i][j, p] = megalamb[k];

                    k = 0;
                    for (int q = 0; q < n; q++)
                    {
                        k = 0;
                        for (int i = 0; i < N; i++)
                        {
                            for (int j = 0; j < d[i]; j++)
                            {
                                for (int p = 0; p <= P[i]; p++)
                                {
                                    ksi[i][j][q] += lamb[i][j, p] * T[q][k];
                                    if (p != 0)
                                        k++;
                                }
                                ksi[i][j][q] = Math.Exp(ksi[i][j][q]) - 1;
                            }
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < N; i++)
                    {
                        ksi[i] = new double[d[i]][];
                        double[] megalamb = new double[P[i] * d[i]];
                        double[][] T = new double[n][];
                        int k = 0;
                        for (int q0 = 0; q0 < n; q0++)
                        {
                            T[q0] = new double[P[i] * d[i]];
                            k = 0;
                            for (int j = 0; j < d[i]; j++)
                            {
                                ksi[i][j] = new double[n];
                                for (int p = 0; p < P[i]; p++)
                                {
                                    if (Vlasn)
                                        T[q0][k] = Math.Cosh(f(p + 1, (2 + Z[i][q0, j]) / 4));
                                    else
                                        T[q0][k] = 1 + f(p + 1, (2 + Z[i][q0, j]) / 4);
                                    if (T[q0][k] <= 0)
                                        T[q0][k] = e;
                                    T[q0][k] = Math.Log(T[q0][k]);
                                    k++;
                                }
                            }
                        }
                        megalamb = SpragenGradKvadrFunction(megalamb, T, B);

                        k = 0;
                        for (int j = 0; j < d[i]; j++)
                        {
                            lamb[i][j, 0] = 0;
                            for (int p = 1; p <= P[i]; p++, k++)
                            {
                                lamb[i][j, p] = megalamb[k];
                            }
                        }
                        k = 0;
                        for (int q0 = 0; q0 < n; q0++)
                        {
                            k = 0;
                            for (int j = 0; j < d[i]; j++)
                            {
                                for (int p = 0; p <= P[i]; p++)
                                {
                                    ksi[i][j][q0] += lamb[i][j, p] * T[q0][k];
                                    if (p != 0)
                                        k++;
                                }
                                ksi[i][j][q0] = Math.Exp(ksi[i][j][q0]) - 1;
                            }
                        }
                    }
                }
            }
            else
            {
                if (!checkBox1.Checked)
                {
                    double[] megalamb = new double[(P[0] + 1) * d[0] + (P[1] + 1) * d[1] + (P[2] + 1) * d[2]];
                    double[][] T = new double[n][];
                    int k = 0;
                    for (int q0 = 0; q0 < n; q0++)
                    {
                        T[q0] = new double[(P[0] + 1) * d[0] + (P[1] + 1) * d[1] + (P[2] + 1) * d[2]];
                        k = 0;
                        for (int i = 0; i < N; i++)
                        {
                            ksi[i] = new double[d[i]][];
                            for (int j = 0; j < d[i]; j++)
                            {
                                ksi[i][j] = new double[n];
                                for (int p = 0; p <= P[i]; p++)
                                {
                                    if (p == 0)
                                        if (Vlasn)
                                            T[q0][k] = Math.Log(Math.Cosh(1));
                                        else
                                            T[q0][k] = Math.Log(1.5);
                                    else
                                    {
                                        if (Vlasn)
                                            T[q0][k] = Math.Cosh(f(p, Z[i][q0, j]));
                                        else
                                            T[q0][k] = 1 + f(p, Z[i][q0, j]);
                                        if (T[q0][k] <= 0)
                                        {
                                            T[q0][k] = e;
                                        }
                                        T[q0][k] = Math.Log(T[q0][k]);
                                    }
                                    k++;
                                }
                            }
                        }
                    }

                    megalamb = SpragenGradKvadrFunction(megalamb, T, B);

                    k = 0;
                    for (int i = 0; i < N; i++)
                        for (int j = 0; j < d[i]; j++)
                            for (int p = 0; p <= P[i]; p++, k++)
                                lamb[i][j, p] = megalamb[k];

                    k = 0;
                    for (int q = 0; q < n; q++)
                    {
                        k = 0;
                        for (int i = 0; i < N; i++)
                        {
                            for (int j = 0; j < d[i]; j++)
                            {
                                for (int p = 0; p <= P[j]; p++, k++)
                                {
                                    {
                                        ksi[i][j][q] += lamb[i][j, p] * T[q][k];
                                    }
                                    ksi[i][j][q] = Math.Exp(ksi[i][j][q]) - 1;
                                }
                            }
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < N; i++)
                    {
                        ksi[i] = new double[d[i]][];
                        double[] megalamb = new double[(P[i] + 1) * d[i]];
                        double[][] T = new double[n][];
                        int k = 0;
                        for (int q0 = 0; q0 < n; q0++)
                        {
                            T[q0] = new double[(P[i] + 1) * d[i]];
                            k = 0;
                            for (int j = 0; j < d[i]; j++)
                            {
                                ksi[i][j] = new double[n];
                                for (int p = 0; p <= P[i]; p++)
                                {
                                    if (p == 0)
                                        if (Vlasn)
                                            T[q0][k] = Math.Log(Math.Cosh(1));
                                        else
                                            T[q0][k] = Math.Log(1.5);
                                    else
                                    {
                                        if (Vlasn)
                                            T[q0][k] = Math.Cosh(f(p, Z[i][q0, j]));
                                        else
                                            T[q0][k] = 1 + f(p, Z[i][q0, j]);
                                        if (T[q0][k] <= 0)
                                        {
                                            T[q0][k] = e;
                                        }
                                        T[q0][k] = Math.Log(T[q0][k]);
                                    }
                                    k++;
                                }
                            }
                        }
                        megalamb = SpragenGradKvadrFunction(megalamb, T, B);

                        k = 0;
                        for (int j = 0; j < d[i]; j++)
                        {
                            for (int p = 0; p <= P[i]; p++, k++)
                            {
                                lamb[i][j, p] = megalamb[k];
                            }
                        }
                        k = 0;
                        for (int q0 = 0; q0 < n; q0++)
                        {
                            k = 0;
                            for (int j = 0; j < d[i]; j++)
                            {
                                for (int p = 0; p <= P[i]; p++, k++)
                                {
                                    ksi[i][j][q0] += lamb[i][j, p] * T[q0][k];
                                }
                                ksi[i][j][q0] = Math.Exp(ksi[i][j][q0]) - 1;
                            }
                        }
                    }
                }
            }
            if (Flag)
            {
                Result.Text += "\r\n Matrix λ:\r\n";
                for (int i = 0; i < N; i++)
                {
                    Result.Text += "\t ||λ" + (i + 1) + "||:\r\n";
                    for (int j = 0; j < d[i]; j++)
                    {
                        for (int p = 0; p <= P[i]; p++)
                        {
                            Result.Text += "\t" + (lamb[i][j, p]).ToString("F4");
                        }
                        Result.Text += "\r\n";
                    }
                }
            }
            return lamb;
        }
        private double[][][][][] PsiShow(double[][][,] lamb, int[] P, double[][][][][] coef)
        {
            if (Flag)
            {
                for (int m = 0; m < dy; m++)
                {
                    Result.Text += "\r\n ψ" + (m + 1) + ":\r\n";
                    for (int i = 0; i < N; i++)
                    {
                        for (int j = 0; j < d[i]; j++)
                        {
                            Result.Text += "ψ" + (i + 1) + "," + (j + 1) + " = exp (";
                            for (int p = 0; p <= P[i]; p++)
                            {
                                if (p > 0)
                                {
                                    if (lamb[m][i][j, p] > 0)
                                    {
                                        Result.Text += "+";
                                    }
                                }
                                if (Vlasn)
                                    Result.Text += lamb[m][i][j, p].ToString("F4") + " *ln(1 + Csh (T" + p + "(x" + (i + 1) + "," + (j + 1) + ") ) ) ";
                                else
                                    Result.Text += lamb[m][i][j, p].ToString("F4") + " *ln(1 + T" + p + "(x" + (i + 1) + "," + (j + 1) + ") ) ";
                            }
                            Result.Text += " - 1 ) \r\n";
                        }
                    }
                }
            }
            if (!Vlasn)
            {
                Result.Text += "\r\nIn the form of polynomial\r\n";
                Pfunk Pf = new Pfunk(PolynomChebishev);
                if (PolinoType.SelectedIndex == 1)
                    Pf = new Pfunk(PolynomChebishev2);
                if (PolinoType.SelectedIndex == 2)
                    Pf = new Pfunk(PolynomChebishev3);
                for (int m = 0; m < dy; m++)
                {
                    for (int i = 0; i < N; i++)
                    {
                        for (int j = 0; j < d[i]; j++)
                        {
                            double[][] tmp = new double[P[i] + 1][];
                            coef[m][i][j] = Pf(P[i]);
                        }
                    }
                }

                if (Flag)
                {
                    for (int m = 0; m < dy; m++)
                    {
                        Result.Text += (m + 1) + ":\r\n\r\n";
                        for (int i = 0; i < N; i++)
                        {
                            for (int j = 0; j < d[i]; j++)
                            {
                                coef[m][i][j][0][0] = 1.5;
                                Result.Text += "ψ" + (i + 1) + "," + (j + 1) + "(x" + (i + 1) + "," + (j + 1) + ") =  ";
                                for (int k = 0; k <= P[i]; k++)
                                {
                                    if (k > 0)
                                        Result.Text += " * ";
                                    for (int p = 0; p < k + 1; p++)
                                    {
                                        if (p < 2)
                                        {
                                            if (p == 0)
                                            {
                                                if (k > 0)
                                                    coef[m][i][j][k][0] += 1;
                                                Result.Text += " ( " + coef[m][i][j][k][0].ToString("F4") + "  ";
                                            }
                                            else
                                            {
                                                Result.Text += " + " + coef[m][i][j][k][1].ToString("F4") + " * x(" + (i + 1) + "," + (j + 1) + ")  ";
                                            }
                                        }
                                        else
                                            Result.Text += " + " + coef[m][i][j][k][p].ToString("F4") + " * x(" + (i + 1) + "," + (j + 1) + ") ^" + p + "  ";
                                    }
                                    Result.Text += ") ^" + lamb[m][i][j, k].ToString("F4") + "  ";
                                }
                                Result.Text += "- 1 \r\n";
                            }
                        }
                    }
                }
            }
            return coef;
        }
        private double[][][] ASearch(double[][][] a)
        {
            for (int m = 0; m < dy; m++)
            {
                for (int i = 0; i < N; i++)
                {

                    double[] megalamb = new double[d[i]];
                    double[][] T = new double[n][];
                    int k = 0;
                    for (int q0 = 0; q0 < n; q0++)
                    {
                        T[q0] = new double[d[i]];
                        k = 0;
                        for (int j = 0; j < d[i]; j++)
                        {
                            if (Vlasn)
                                T[q0][k] = Math.Cosh(ksi[i][j][q0]);
                            else
                                T[q0][k] = ksi[i][j][q0] + 1;
                            if (T[q0][k] <= 0)
                            {
                                T[q0][k] = e;
                            }
                            T[q0][k] = Math.Log(T[q0][k]);
                            k++;
                        }
                    }
                    megalamb = SpragenGradKvadrFunction(megalamb, T, logZY[m]);
                    for (int j = 0; j < d[i]; j++)
                        a[m][i][j] = megalamb[j];
                }
            }
            if (Flag)
            {
                Result.Text += "\r\n Matrix ||a||:\r\n";

                for (int m = 0; m < dy; m++)
                {
                    Result.Text += "Matrix ||a" + (m + 1) + "||:\r\n";
                    for (int i = 0; i < N; i++)
                    {
                        for (int j = 0; j < d[i]; j++)
                        {
                            Result.Text += "\t" + a[m][i][j].ToString("F4");
                        }
                        Result.Text += "\r\n";
                    }
                }
                Result.Text += "\r\n";
            }
            return a;
        }
        private void ФimShow(double[][][] a)
        {
            Result.Text += "\r\n Ф:\r\n";
            Result.Text += "In the form ψ \r\n";
            for (int i = 0; i < N; i++)
            {

                for (int m = 0; m < dy; m++)
                {
                    Result.Text += "Ф" + (i + 1) + "," + (m + 1) + " = ";
                    for (int j = 0; j < d[i]; j++)
                    {
                        if (Vlasn)
                            Result.Text += "(Csh (ψ" + (i + 1) + "," + (j + 1) + ") + 1) ^" + a[m][i][j].ToString("F4") + "  ";
                        else
                            Result.Text += "(ψ" + (i + 1) + "," + (j + 1) + " + 1) ^" + a[m][i][j].ToString("F4") + "  ";
                    }
                    Result.Text += "- 1  \r\n";
                }
            }
        }
        private double[][] CSearch(double[][][] a, double[][] c)
        {
            diff = new double[dy];
            vlnovFunc = new double[dy, n];
            for (int m = 0; m < dy; m++)
            {

                double[][] T = new double[n][];
                int k = 0;
                for (int q0 = 0; q0 < n; q0++)
                {
                    T[q0] = new double[N];
                    k = 0;
                    for (int i = 0; i < N; i++, k++)
                    {
                        double tmp = 1;
                        for (int j = 0; j < d[i]; j++)
                            if (Vlasn)
                                tmp *= Math.Pow(Math.Cosh(ksi[i][j][q0]), a[m][i][j]);
                            else
                                tmp *= Math.Pow(1 + ksi[i][j][q0], a[m][i][j]);
                        if (Vlasn)
                            tmp = Math.Cosh(tmp - 1);
                        if (tmp <= 0)
                        {
                            tmp = e;
                        }
                        T[q0][k] = Math.Log(tmp);
                    }
                }
                c[m] = SpragenGradKvadrFunction(c[m], T, logZY[m]);

                k = 0;
                for (int q = 0; q < n; q++)
                {
                    k = 0;
                    for (int i = 0; i < N; i++, k++)
                    {
                        vlnovFunc[m, q] += c[m][i] * T[q][k];
                    }
                    vlnovFunc[m, q] = (Math.Exp(vlnovFunc[m, q]) - 1);
                    vlnovFunc[m, q] = ZY[m][q] - (ZY[m][q] - vlnovFunc[m, q]) / 4;
                    if (Math.Abs(vlnovFunc[m, q] - ZY[m][q]) > diff[m])
                        diff[m] = Math.Abs(vlnovFunc[m, q] - ZY[m][q]);
                }
            }

            if (Flag)
            {
                textBox1.Text = "";
                for (int m = 0; m < dy; m++)
                {
                    textBox1.Text += "Max in Y" + (m + 1) + ": " + (Math.Abs(diff[m]) /** (MaxY[m] - MinY[m])*/).ToString("F4") + "\r\n";
                }
                Result.Text += "\r\n Matrix ||c||:\r\n";
                for (int m = 0; m < dy; m++)
                {
                    for (int i = 0; i < N; i++)
                        Result.Text += "\t" + c[m][i].ToString("F4");
                    Result.Text += "\r\n";
                }
                Result.Text += "\r\n";
            }
            return c;
        }
        private void ФShow(double[][] c, int[] P)
        {
            Result.Text += "\r\n Ф:\r\n";
            Result.Text += "In the form Фі: \r\n";
            for (int m = 0; m < dy; m++)
            {
                Result.Text += "Ф" + (m + 1) + " =  ";
                for (int i = 0; i < N; i++)
                {
                    if (Vlasn)
                        Result.Text += "(1 + Csh (Ф" + (i + 1) + "," + (m + 1) + ") ) ^" + c[m][i].ToString("F4") + "  ";
                    else
                        Result.Text += "(1 + Ф" + (i + 1) + "," + (m + 1) + " ) ^" + c[m][i].ToString("F4") + "  ";
                }
                Result.Text += " - 1 \r\n";
            }
        }
        private void ReNorm(double[][] c, int[] P)
        {
            Result.Text += "\r\n  ";
            Result.Text += "\r\n Ф without normalization for different Yi:\r\n";
            double temp;
            for (int m = 0; m < dy; m++)
            {
                Result.Text += "Ф" + (m + 1) + " =  " + (MaxY[m] - MinY[m]).ToString("F4");
                for (int i = 0; i < N; i++)
                {
                    Result.Text += "(1 + Ф" + (i + 1) + "," + (m + 1) + " ) ^" + c[m][i].ToString("F4") + "  ";
                }
                temp = -MaxY[m] + 2 * MinY[m];
                if (temp > 0)
                    Result.Text += " + ";
                Result.Text += temp.ToString("F4") + "\r\n";
            }
        }
        private void Start_Click(object sender, EventArgs e)
        {
            Xinput.Enabled = false;
            Yinput.Enabled = false;
            addbutton.Enabled = false;
            xload.Enabled = false;
            Range.Enabled = false;
            dim1.Enabled = false;
            dim2.Enabled = false;
            dim3.Enabled = false;
            dimy.Enabled = false;
            StPolynom.Enabled = true;

            if (VlasniiMethod.Checked)
                Vlasn = true;
            else
                Vlasn = false;
            Flag = true;
            int[] P = new int[N];
            P[0] = Convert.ToInt32(Rankx1.Value);
            P[1] = Convert.ToInt32(Rankx2.Value);
            P[2] = Convert.ToInt32(Rankx3.Value);
            for (int i = 0; i < N; i++)
                Result.Text += "P" + (i + 1) + " = " + P[i] + "  ";
            Result.Text += "\r\n";
            Result.Text += "Polynomial " + PolinoType.Text.ToString() + "\r\n\r\n";

            ReadData();


            StreamWriter Out = new StreamWriter(OStream);
            MinMaxFound();


            Normalizing();

            double[][][,] lamb = new double[dy][][,];
            for (int m = 0; m < dy; m++)
            {
                lamb[m] = new double[N][,];
                for (int i = 0; i < N; i++)
                {
                    lamb[m][i] = new double[d[i], P[i] + 1];
                }
            }

            double[][] B = new double[dy][];
            for (int m = 0; m < dy; m++)
            {
                B[m] = new double[n];
                for (int i = 0; i < n; i++)
                {
                    if (radioButton_normy.Checked)
                        B[m][i] = logZY[m][i];
                    //B[m][i] = ZY[m][i];
                    else
                        B[m] = Bq();
                }
            }

            for (int m = 0; m < dy; m++)
            {
                Result.Text += (m + 1) + ":\r\n\r\n";
                lamb[m] = LambdaSearch(B[m], P, lamb[m]);
            }

            double[][][][][] coef = new double[dy][][][][];
            for (int m = 0; m < dy; m++)
            {
                coef[m] = new double[N][][][];
                for (int i = 0; i < N; i++)
                {
                    coef[m][i] = new double[d[i]][][];
                    for (int j = 0; j < d[i]; j++)
                    {
                        coef[m][i][j] = new double[P[i] + 1][];
                        for (int p = 0; p < P[i] + 1; p++)
                        {
                            coef[m][i][j][p] = new double[P[i] + 1];
                        }
                    }
                }
            }
            coef = PsiShow(lamb, P, coef);

            double[][][] a = new double[dy][][];
            for (int i = 0; i < dy; i++) a[i] = new double[N][];
            for (int i = 0; i < dy; i++)
                for (int j = 0; j < N; j++)
                    a[i][j] = new double[d[j]];
            a = ASearch(a);


            ФimShow(a);

            double[][] c = new double[dy][];
            for (int m = 0; m < dy; m++)
            {
                c[m] = new double[N];
            }
            c = CSearch(a, c);

            ФShow(c, P);

            ReNorm(c, P);

            Approximate = true;
            Draw(Convert.ToInt32(numericUpDown1.Value - 1));

            XStream.Close();
            XStream = openFileDialog1.OpenFile();
            YStream.Close();
            YStream = openFileDialog2.OpenFile();


            Out.Write(Result.Text);
            Out.Close();
            OStream.Close();
        }
        void Draw(int m)
        {

            gPanel.Clear(Color.White);
            Font segoeUI = new Font("Segoe UI", 6);

            Pen p = new Pen(Color.Black, 2);
            gPanel.DrawLine(p, new Point(3, 3), new Point(3, panel1.Height - 3));
            gPanel.DrawLine(p, new Point(3, panel1.Height - 3), new Point(panel1.Width - 3, panel1.Height - 3));

            p = new Pen(Color.Green, 1);
            gPanel.DrawString("Approximation", segoeUI, new SolidBrush(Color.Green), new PointF(panel1.Width - 100 + 33, 20), new StringFormat());
            gPanel.DrawLine(p, new Point(panel1.Width - 100, 30),
            new Point(panel1.Width - 100 + 30, 30));

            for (int q0 = 0; q0 < n - 1; q0++)
            {
                gPanel.DrawLine(p, new Point(q0 * (panel1.Width / n) + 3, (int)(panel1.Height * (1 - vlnovFunc[m, q0])) - 3),
                new Point((q0 + 1) * (panel1.Width / n) + 3, (int)(panel1.Height * (1 - vlnovFunc[m, q0 + 1])) - 3));
            }

            p = new Pen(Color.Black, 1);
            gPanel.DrawString("Function", segoeUI, new SolidBrush(Color.Black), new PointF(panel1.Width - 100 + 33, 30), new StringFormat());
            gPanel.DrawLine(p, new Point(panel1.Width - 100, 40),
            new Point(panel1.Width - 100 + 30, 40));

            for (int q0 = 0; q0 < n - 1; q0++)
            {
                gPanel.DrawLine(p, new Point(q0 * (panel1.Width / n) + 3, (int)(panel1.Height * (1 - ZY[m][q0])) - 3),
                new Point((q0 + 1) * (panel1.Width / n) + 3, (int)(panel1.Height * (1 - ZY[m][q0 + 1])) - 3));
            }

        }
        private void dimy_ValueChanged(object sender, EventArgs e)
        {
            numericUpDown1.Maximum = dimy.Value;
        }
        private void numericUpDown1_ValueChanged(object sender, EventArgs e)
        {
            if (Approximate)
                Draw(Convert.ToInt32(numericUpDown1.Value - 1));
        }
        private void Form1_Load(object sender, EventArgs e)
        {

        }
        private void StPolynom_Click(object sender, EventArgs e)
        {
            int amount = 7;
            int[] P = new int[N];
            int[][] Pmax = new int[dy][];
            double[] max = new double[dy];
            if (VlasniiMethod.Checked)
                Vlasn = true;
            else
                Vlasn = false;
            for (int m = 0; m < dy; m++)
            {
                Pmax[m] = new int[N];
                max[m] = 1;
            }
            Flag = false;
            for (int p1 = 1; p1 < amount; p1++)
            {
                for (int p2 = 1; p2 < amount; p2++)
                {
                    for (int p3 = 1; p3 < amount; p3++)
                    {
                        P[0] = p1;
                        P[1] = p2;
                        P[2] = p3;
                        double[][][,] lamb = new double[dy][][,];
                        for (int m = 0; m < dy; m++)
                        {
                            lamb[m] = new double[N][,];
                            for (int i = 0; i < N; i++)
                            {
                                lamb[m][i] = new double[d[i], P[i] + 1];
                            }
                        }

                        double[][] B = new double[dy][];
                        for (int m = 0; m < dy; m++)
                        {
                            B[m] = new double[n];
                            for (int i = 0; i < n; i++)
                            {
                                if (radioButton_normy.Checked)
                                    B[m][i] = logZY[m][i];
                                else
                                    B[m] = Bq();
                            }
                        }
                        for (int m = 0; m < dy; m++)
                        {
                            lamb[m] = LambdaSearch(B[m], P, lamb[m]);
                        }
                        double[][][] a = new double[dy][][];
                        for (int i = 0; i < dy; i++) a[i] = new double[N][];
                        for (int i = 0; i < dy; i++)
                            for (int j = 0; j < N; j++)
                                a[i][j] = new double[d[j]];
                        a = ASearch(a);

                        double[][] c = new double[dy][];
                        for (int m = 0; m < dy; m++)
                        {
                            c[m] = new double[N];
                        }
                        c = CSearch(a, c);

                        for (int m = 0; m < dy; m++)
                            if (diff[m] < max[m])
                            {
                                max[m] = diff[m];
                                for (int i = 0; i < N; i++)
                                    Pmax[m][i] = P[i];
                            }
                    }
                }
            }
            Result.Text = "The best power: \r\n";
            for (int m = 0; m < dy; m++)
            {
                for (int i = 0; i < N; i++)
                    Result.Text += "P" + (i + 1) + " = " + Pmax[m][i] + " \r\n";
                Result.Text += "Max y" + (m + 1) + " = " + max[m] + "\r\n";
            }
            using (System.IO.StreamWriter file =
                      new System.IO.StreamWriter(@outfile.Text))
            {
                foreach (char line in Result.Text)
                {
                    file.Write(line);
                }
            }
        }

        private void button1_Click(object sender, EventArgs e)
        {

        }

        private void but_downForecast_Click(object sender, EventArgs e)
        {
            if (FileForecast.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    if ((XStream = openFileDialog1.OpenFile()) != null)
                    {
                        //Box_fileForcast.Text = openFileDialog1.SafeFileName;
                    }
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Unable to read a file. " + ex.Message);
                }
            }
        }
    }
}

