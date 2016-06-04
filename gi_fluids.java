/* 
 
 Program to calculate the transit along GI tract.
 Implementing a Poisson point process to determine the 
 time between emptying events of random fluid packet
 
 A. Talattof
 (C) 03.02.2016
 
 */

import java.util.*;
import java.io.*;
import java.util.Date;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import org.apache.commons.math3.stat.inference.*;
import org.apache.commons.math3.stat.descriptive.*;
import org.apache.commons.math3.stat.StatUtils;

class gi_fluids {
    
//    static double a = 20.0;
//    static double b = 0.10583135331707194;
//    static double c = 34;
//    static double tau = 51;
//    static double a = 44.48223503966336;
//    static double b = 0.04351428717154171/60;
//    static double c = 3.4364748697778102;
//    static double tau = 38.66479268738968*60;
    static double a = 40.0;
    static double b = 0.045/60;
    static double c = 1.5;
    static double tau = 59*60.0;
    static double alpha = 0.2;
    static double beta = 1000;
    static double gamma = 0.25;
    static double VS = 10.0;
    static double s = 0.01;
    static double p1 = 0.11;
    static double p2 = 8.0*Math.pow(10,-7);
    static double p3 = 8.4;
    static double phi = 1.0/3600.0;
    
    
//    static double fwd = 29.46296368613547;
//    static double mix = 0.7535128891400147;
    static double fwd = 35;
    static double mix = 0.85;
    
    public static double lmbd_stomach(double vol, double t) {
        
        double Time = t%(120*60) - tau;
        
        //double lmbd = c*Math.exp(-alpha*(vol-VS)/(float) VS)*(1/(1+a*Math.exp(-b*Time)))  + 0.00001;
        //double lmbd = ((c/2)/(1 + 100*Math.exp(-.5*(vol - VS)/VS)))*(1/(1+a*Math.exp(-b*Time)))  + 0.000001;
        double vol_scaling = ((c/2)/(1 + 35*Math.exp(-.5*(vol - VS)/VS)));
        double sum = 0.0,lmbd = 0.0;
        for (int k = 1; k < 26; k++) {
            sum += Math.pow(-1,k)*Math.sin(-phi*Math.PI*k*Time)/k + p1;
        }
        lmbd = vol_scaling*(p2*Math.pow(sum,p3) + s);
        return lmbd;
    }
    
    public static double lmbd_SI(double vol, double t) {
        
        double Time = t%(120*60) - tau;
        
        // double lmbd = c*Math.exp(-alpha*(150.0-VS)/(float) VS)*(1-1/(1+a*Math.exp(-b*Time)))  + 0.00001;
        // double lmbd = ((c/2)/(1 + 100*Math.exp(-.5*(150.0 - VS)/VS)))*(1/(1+a*Math.exp(-b*Time)))  + 0.000001;
        double vol_scaling = ((c/2)/(1 + 35*Math.exp(-.5*(50 - VS)/VS)));
        double sum = 0.0,lmbd = 0.0;
        for (int k = 1; k <= 25; k++) {
            sum += Math.pow(-1,k)*Math.sin(-phi*Math.PI*k*Time)/k + p1;
        }
        lmbd = vol_scaling*(p2*Math.pow(sum,p3) + s);
        return lmbd;
    }
    
    public static double[] linspace(double min, double max, int points) {
        double[] d = new double[points];
        for (int i = 0; i < points; i++){
            d[i] = min + i * (max - min) / (points - 1);
        }
        return d;
    }

    public static double[] gastric_emptying(int t0, int TT, double Dt, double X0,
                                    int N, double tD) {
        
        double alpha = 0.2;
        double beta = 1000;
        double gamma = 0.25;
        
        double[] t = linspace(t0,TT,N+1);
        
        double packet_vol = 8.0;
        
        double[] X = new double[N+1];
        X[0] = X0;
        
        double t_s = 0;
        double t_int = 0;
        
        Random randomno = new Random();
        
        double Xfwd = 0;
        
        for (int i = 1; i < N+1; i++ ) {
            
            // Stomach
            double pkt_fwd = 0;
            
            double lambda_stomach = lmbd_stomach(X[i-i],t[i-1]+tD);
            double t_s_next = Math.log(1-randomno.nextDouble())/(-lambda_stomach);
            
            if (t_s + t_s_next < t[i]) {
                t_s = t[i];
            }
            if (Math.abs(t[i]-(t_s+t_s_next))<=Dt) {
                t_s = t_s + t_s_next;
                double var = 2/(1+Math.exp(-1.5*(X[i-1]-10)));
                double Z = randomno.nextGaussian()*var;
                
                pkt_fwd = ((Z + packet_vol)/(1+beta*Math.exp(-gamma*(X[i-1] + (Z+packet_vol)))));
            }
            Xfwd = pkt_fwd;
            if (X[i-1]-pkt_fwd<0) {
                Xfwd = 0;
            }
            X[i] = X[i-1] - Xfwd;
        }
        return X;
    }
    
    public static Pair poisson_proc(int t0, int TT, double Dt, double[] X0, double[] Y0,
                                      int N, double tD, int NCOMPS) {
        
        double alpha = 0.2;
        double beta = 1000;
        double gamma = 0.25;
        
        double[] t = linspace(t0,TT,N+1);
        
        double[] KShifts = new double[NCOMPS];
        double[] KFactors = new double[NCOMPS];
        for (int i = 0; i < NCOMPS; i++) {
            KShifts[i] = i*0.3/15;
        }
        KFactors[0] = 1;
        for (int i = 1; i < NCOMPS; i++) {
            KFactors[i] = (float) 1+Math.log((fwd-(float) i*1.2)/((float) fwd));
        }
        
        // Gastric secretion per interval (ml/min)
        double gastric_secretion = 1.0/60;
        // Duodenal secretion per interval (ml/min)
        double duod_secretion = 0.333/60;
        
        // Fluid reabsorption absorption per interval (ml/mn)
        // Stomach - 0 ml/min
        // Duodenum - 8 ml/min
        // Jejunum - 12 ml/min
        // Ileum - 2 ml/mn
        // Colon = 0 ml/min
        double stom = 0;
        double duod = 8/60/3;
        double jej = 12/60/4;
        double ileum = 2/60/3;
        double col = 0;
        double[] fluid_reabsoprtion = {stom,
            duod,duod,duod,
            jej,jej,jej,jej,
            ileum,ileum,ileum,
            col,col,col,
            0};
        for (int i=0;i<fluid_reabsoprtion.length;i++) {
            fluid_reabsoprtion[i] = fluid_reabsoprtion[i];
        }
        
        double packet_vol = 8.0;
        
        double[][] X = new double[N+1][NCOMPS];
        for (int i=0;i<NCOMPS;i++){
            X[0][i] = X0[i];
        }
        double[][] Y = new double[N+1][NCOMPS];
        for (int i=0;i<NCOMPS;i++){
            Y[0][i] = Y0[i];
        }
        
        double t_s = 0;
        double t_int = 0;
        double t_mix = 0;
        
        Random randomno = new Random();
        
        double[] Xfwd = new double[NCOMPS];
        double[] Xmix = new double[NCOMPS];
        double[] pct_fwd = new double[NCOMPS];
        double[] pct_back = new double[NCOMPS];
        
        Xfwd[NCOMPS-1] = 0;
        Xmix[NCOMPS-1] = 0;
        Xmix[1] = 0;
        Xmix[4] = 0;
        Xmix[8] = 0;
        Xmix[11] = 0;
        
        for (int i = 1; i < N+1; i++ ) {
            
            // Stomach
            int IND = 0;
            double pkt_fwd = 0;
            double pkt_mix = 0;
            double f_abs = 0;

            double lambda_stomach = lmbd_stomach(X[i-i][IND],t[i-1]+tD);
            double t_s_next = Math.log(1-randomno.nextDouble())/(-lambda_stomach);
            
            if (t_s + t_s_next < t[i]) {
                t_s = t[i];
            }
            if (Math.abs(t[i]-(t_s+t_s_next))<=Dt) {
                t_s = t_s + t_s_next;
                double var = 2/(1+Math.exp(-1.5*(X[i-1][IND]-10)));
                double Z = randomno.nextGaussian()*var;
                
                pkt_fwd = ((Z + packet_vol)/(1+beta*Math.exp(-gamma*(X[i-1][IND] + (Z+packet_vol)))));
            }
            Xfwd[IND] = pkt_fwd;
            if (X[i-1][IND]-pkt_fwd+gastric_secretion<0) {
                Xfwd[IND] = 0;
            }
            X[i][IND] = X[i-1][IND] - Xfwd[IND] + gastric_secretion;
            
            // Amount of phenol red transiting is proportional to fraction of volume that transits forward/backward
            pct_fwd[IND] = pkt_fwd/X[i-1][IND];
            if (pct_fwd[IND] < 0.0001 || pct_fwd[IND] > 1/0.0001 || X[i-1][IND] == 0) {
                pct_fwd[IND] = 0;
            } else if (Y[i-1][IND]*((float) 1-pct_fwd[IND]) < 0 || ((float) 1-pct_fwd[IND]) < 0) {
                pct_fwd[IND] = 0;
            }
            Y[i][IND] = Y[i-1][IND]*((float) 1-pct_fwd[IND]);
        
            // Forward Transit in Intestinal Compartments
            for (int j = 1; j<14; j++) {
                IND = j;
                pkt_fwd = 0;
                
                double lambda_int = lmbd_SI(X[i-1][IND],t[i-1]+tD-KShifts[IND])/(float) KFactors[IND];
                double t_int_next = Math.log(1-randomno.nextDouble())/(-lambda_int);
                
                if (t_int+t_int_next<t[i]) {
                    t_int = t[i];
                }
                
                if (Math.abs(t[i]-(t_int+t_int_next))<=Dt) {
                    t_int += t_int_next;
                    double var = 2/(1+Math.exp(-1.5*(X[i-1][IND]-10)));
                    double Z = randomno.nextGaussian()*var;
                    pkt_fwd = ((Z + packet_vol)/(1+beta*Math.exp(-gamma*(X[i-1][IND] + (Z+packet_vol)))));
                }
                Xfwd[IND] = pkt_fwd;
                f_abs = (fluid_reabsoprtion[IND]<X[i-1][IND]) ? fluid_reabsoprtion[IND] : 0;
                if (X[i-1][IND]-pkt_fwd+f_abs<0) {
                    Xfwd[IND] = 0;
                }
            
                X[i][IND] = X[i-1][IND] - Xfwd[IND] + Xfwd[IND-1] - f_abs;
                if (j<4) { X[i][IND] = X[i][IND] + duod_secretion; }
                // Amount of phenol red transiting is proportional to fraction of volume that transits forward/backward
                if (X[i-1][IND] == 0) { pct_fwd[IND] = 0; }
                else { pct_fwd[IND] = pkt_fwd/X[i-1][IND]; }
                if (pct_fwd[IND] < 0.0001 || pct_fwd[IND] > 1/0.0001) {
                    pct_fwd[IND] = 0;
                } else if (Y[i-1][IND]*((float) 1-pct_fwd[IND]) + Y[i-1][IND-1]*(pct_fwd[IND-1]) < 0 || ((float) 1-pct_fwd[IND]) < 0) {
                    pct_fwd[IND] = 0;
                }
                Y[i][IND] = Y[i-1][IND]*((float) 1-pct_fwd[IND]) + Y[i-1][IND-1]*(pct_fwd[IND-1]);
                
            }
            X[i][NCOMPS-1] = X[i-1][NCOMPS-1] + Xfwd[NCOMPS-2];
            Y[i][NCOMPS-1] = Y[i-1][NCOMPS-1] + Y[i-1][NCOMPS-2]*pct_fwd[NCOMPS-2];

            // Backward Transit in Intestinal Comparments
            for (int j = 1; j<14; j++) {
                IND = j;
                pkt_mix = 0;
                
                double lambda_int = lmbd_SI(X[i][IND],t[i-1]+tD-KShifts[IND])/((float) KFactors[IND]*mix);
                double t_mix_next = Math.log(1-randomno.nextDouble())/(-lambda_int);
                
                if (t_mix+t_mix_next<t[i]) {
                    t_mix = t[i];
                }
                
                if (Math.abs(t[i]-(t_mix+t_mix_next))<=Dt) {
                    t_mix += t_mix_next;
                    double var = 2/(1+Math.exp(-1.5*(X[i][IND]-10)));
                    double Z = randomno.nextGaussian()*var;
                    pkt_mix = ((Z + (float) packet_vol)/(1+beta*Math.exp(-gamma*(X[i][IND] + (Z+packet_vol)))));
                }
                Xmix[IND] = pkt_mix;
                if (X[i][IND]-pkt_mix<0 || IND==1 || IND==4 || IND==8 || IND==11) {
                    Xmix[IND] = 0;
                }
                
                X[i][IND] = X[i][IND] - Xmix[IND] + Xmix[IND+1];
                // Amount of phenol red transiting is proportional to fraction of volume that transits forward/backward
                if (X[i][IND] == 0) { pct_back[IND] = 0; }
                else { pct_back[IND] = pkt_mix/X[i][IND]; }
                if (pct_back[IND] < 0.0001 || pct_back[IND] > 1/0.0001) {
                    pct_back[IND] = 0;
                } else if (Y[i][IND]*((float) 1-pct_back[IND]) + Y[i-1][IND+1]*(pct_back[IND+1]) < 0 || ((float) 1-pct_back[IND]) < 0) {
                    pct_back[IND] = 0;
                }
                
                Y[i][IND] = Y[i][IND]*((float) 1-pct_back[IND]) + Y[i-1][IND+1]*(pct_back[IND+1]);
            }
        }
        return new Pair(X,Y);
    }
    
    public static double[][] append(double[][] a, double[][] b) {
        double[][] result = new double[a.length + b.length][];
        System.arraycopy(a, 0, result, 0, a.length);
        System.arraycopy(b, 0, result, a.length, b.length);
        return result;
    }
    
    public static int getArrayIndex(double[] arr,double mu, double eps) {
        
        int k=0;
        for(int i=0;i<arr.length;i++){
            
//            if(Math.abs(arr[i]-mu)<=eps){
            if(arr[i]>=mu){
                k=i;
                break;
            }
        }
        return k;
    }
    
    public static void saveResults(double[][] X, String title) throws Exception {
        System.out.println("\n--- Writing results to file: " + title);
        FileWriter fw;
        try {
            fw = new FileWriter(new File(title));
            BufferedWriter bw = new BufferedWriter (fw);
            String tempStr = "";
            for (int i = 0; i < X.length; i++) {
                for (int j = 0; j < X[i].length; j++) {
                    bw.write(String.valueOf(X[i][j]) + ",");
                }
                bw.newLine();
            }
            bw.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    public static void run2(double dose_vol, int ITS) throws Exception,FileNotFoundException {
        
        String anim= "|/-\\";
        
        int TT = 14400;
        int N = TT*4;
        int t0 = 0;
        double Dt =(TT-t0)/(float) N;
        double X0 = dose_vol+35;
        double tD;
        
        double[] dose_times = linspace(0,120*60,121);
        int totIts = ITS*dose_times.length;
        double[][] X1 = new double[totIts][N];
        Date date = new Date();
        String start_time = date.toString();
        System.out.println("--- " + start_time);
        System.out.println("--- Total iterations: " + totIts);
        
        for (int i = 0; i<dose_times.length; i++) {
            tD = dose_times[i];
            for (int k = 0; k<ITS; k++) {
                
                int index = ITS*i + k;
                double pct_complete = Math.round((float) (index)/totIts*100*10.0)/10.0;
                
                String data = "\r" + "--- Calculating Results " + anim.charAt(index % anim.length()) + " "  + index + "/" + totIts + " [" + pct_complete + "%]" ;
                System.out.write(data.getBytes());
                
                double[] X = gastric_emptying(t0, TT, Dt, X0, N, tD);
                X1[index] = X;
            }
        }
        saveResults(X1,"results_emptying.txt");
    }
    
    public static void run(double dose, double dose_vol, int ITS) throws Exception,FileNotFoundException {
        
        String anim= "|/-\\";
        
        int NCOMPS = 15;
        int TT = 28800;
        int N = TT*4;
        int t0 = 0;
        double Dt =(TT-t0)/(float) N;
        double[] X0 = {dose_vol+35,
            10, 10, 10,
            10, 10, 10, 10,
            5, 5, 5,
            0, 0, 0,
            0
        };
        double[] Y0 = {dose,
            0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0,
            0, 0, 0,
            0
        };
        double tD;
        
        double[] dose_times = linspace(0,120*60,20);
        String[] phase_lbls = {"Phase_I.txt",
            "Phase_II.txt",
            "Phase_III.txt",
            "Phase_III_Late.txt"};
        int totIts = ITS*dose_times.length;
        double[][][] Z1 = new double[ITS*dose_times.length][N][NCOMPS];
        double[][][] Z2 = new double[ITS*dose_times.length][N][NCOMPS];
        double[][] X1 = null;
        double[][] X2 = null;
        
        Date date = new Date();
        String start_time = date.toString();
        System.out.println("--- " + start_time);
        System.out.println("--- Total iterations: " + totIts);
        
        for (int i = 0; i<dose_times.length; i++) {
            tD = dose_times[i];
            for (int k = 0; k<ITS; k++) {
                
                int index = ITS*i + k;
                double pct_complete = Math.round((float) (index)/totIts*100*10.0)/10.0;
                
                String data = "\r" + "--- Calculating Results " + anim.charAt(index % anim.length()) + " "  + index + "/" + totIts + " [" + pct_complete + "%]" ;
                System.out.write(data.getBytes());
                
                Pair pair = poisson_proc(t0, TT, Dt, X0, Y0, N, tD,  NCOMPS);
                Z1[index] = pair.getArray1();
                if (X1==null) { X1 = Z1[index]; }
                else { X1 = append(X1,Z1[index]); }
                Z2[index] = pair.getArray2();
                if (X2==null) { X2 = Z2[index]; }
                else { X2 = append(X2,Z2[index]); }
                
            }
        }
        saveResults(X1,"results_vols.txt");
        saveResults(X2,"results_PR.txt");
    }
    
    public static void main(String[] args) throws Exception,FileNotFoundException {
        int ITS = 0;
        double dose = 0;
        double dose_vol = 0;
        
        switch (args[0]) {
            case "optimize" :
                dose = Double.parseDouble((args[1]));
                dose_vol = Double.parseDouble((args[2]));
                ITS = Integer.parseInt(args[3]);
                run(dose,dose_vol,ITS);

            case "emptying" :
                dose_vol = Double.parseDouble((args[1]));
                ITS = Integer.parseInt(args[2]);
                run2(dose_vol,ITS);
                
            case "lambda" :
                System.out.println("\n--- Writing results to file: lmbd.txt");
                FileWriter fw;
                int TT = 9000;
                int N = TT;
                int t0 = 0;
                
                double[] t = linspace(t0,TT,N);
                double[] lmbd = new double[N];
                for (int i = 0; i < t.length-1; i++) {
                    lmbd[i] = lmbd_stomach(200,t[i]);
                }
                
                try {
                    fw = new FileWriter(new File("lmbda.txt"));
                    BufferedWriter bw = new BufferedWriter (fw);
                    String tempStr = "";
                    for (int j = 0; j < lmbd.length; j++) {
                        bw.write(String.valueOf(lmbd[j]) + ",");
                    }
                    bw.close();
                } catch (IOException ioe) {
                    ioe.printStackTrace();
                }
        }
    }
}

class Pair extends gi_fluids {
    
    private double[][] array1;
    private double[][] array2;
    public Pair(double[][] array1, double[][] array2)
    {
        this.array1 = array1;
        this.array2 = array2;
        
    }
    public double[][] getArray1() { return array1; }
    public double[][] getArray2() { return array2; }
}