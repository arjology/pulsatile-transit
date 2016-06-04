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
import org.apache.commons.math3.distribution.*;
import org.apache.commons.math3.stat.StatUtils;
import java.text.DecimalFormat;
import java.text.NumberFormat;

class gi_emptying {
    
    public static double[] original_params = {40.0,0.045/60,1.5,59*60.0};
    
    static double alpha = 0.2;
    static double beta = 1000;
    static double gamma = 0.25;
    static double VS = 10.0;
    static double s = 0.01;
    static double p1 = 0.11;
    static double p2 = 8.0*Math.pow(10,-7);
    static double p3 = 8.4;
    static double phi = 1.0/3600.0;
    
//    {44.48223503966336,
//    0.04351428717154171,
//    3.4364748697778102,
//    38.66479268738968};
    
//    {14.368359417242269,
//     0.027822196163989732,
//    3.2206471539143577,
//    91.82153248487211};

    public static double lmbd_stomach(double[] params, double vol, double t) {
        
        double a = params[0];
        double b = params[1]/60;
        double c = params[2];
        double tau = params[3]*60;
        
        double Time = t%(120*60) - tau;
        
        //double lmbd = c*Math.exp(-alpha*(vol-VS)/(float) VS)*(1/(1+a*Math.exp(-b*Time)))  + 0.00001;
        //double lmbd = ((c/2)/(1 + 100*Math.exp(-.5*(vol - VS)/VS)))*(1/(1+a*Math.exp(-b*Time)))  + 0.000001;
        double vol_scaling = ((c/2)/(1 + 35*Math.exp(-.5*(vol - VS)/VS)));
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

    public static double[] thalf(double[] params, double Vo, int tD1, int tD2, int N) {
        
        Random randomno = new Random();
        
        double alpha = 0.2;
        double beta = 1000.0;
        double gamma = 0.25;
        
        double packet_vol = 8.0, t_s = 0, VI = 0;
        double[] t_array = new double[(tD2-tD1)*N];
        
        int ind = 0;
        for (int tD = tD1; tD<tD2; tD++) {
            for (int k = 0; k<N; k++) {
                VI = Vo;
                t_s = 0;
                while (VI > 0.5*Vo) {
                    double pkt_fwd = 0;
                    double lambda_stomach = lmbd_stomach(params,VI,t_s+tD*60);
                    double t_s_next = Math.log(1-randomno.nextDouble())/(-lambda_stomach);
                    double var = 2/(1+Math.exp(-1.5*(VI-10)));
                    double Z = randomno.nextGaussian()*var;
                    pkt_fwd = ((Z + (float) packet_vol)/(1+beta*Math.exp(-gamma*(VI + (Z+packet_vol)))));
                    t_s += t_s_next;
                    VI -= pkt_fwd;
                }
                t_array[ind] = t_s;
                ind++;
            }
        }
        return t_array;
    }

    public static void mean_halftimes(double[] params) {
        
        int N = 10000;
        DecimalFormat f = new DecimalFormat("###.0000");
        
        double Vo = 200;
        double[] t200_pI = thalf(params, Vo, 0, 60, N);
        double[] t200_pII = thalf(params, Vo, 60, 105, N);
        double[] t200_pIII = thalf(params, Vo, 105, 115, N);
        
        Vo = 50;
        double[] t50_pI = thalf(params, Vo, 0, 60, N);
        double[] t50_pII = thalf(params, Vo, 60, 105, N);
        double[] t50_pIII = thalf(params, Vo, 105, 115, N);
        
        String[] res50 = {f.format(StatUtils.geometricMean(t50_pI)/60),
            f.format(StatUtils.geometricMean(t50_pII)/60),
            f.format(StatUtils.geometricMean(t50_pIII)/60)};
        String[] res200 = {f.format(StatUtils.geometricMean(t200_pI)/60),
            f.format(StatUtils.geometricMean(t200_pII)/60),
            f.format(StatUtils.geometricMean(t200_pIII)/60)};

        
        System.out.println("\n--- Emptying Half-Times:");
        System.out.println("\n--- 50 mL");
        System.out.println("Phase I\t\tPhase II\t\tPhase III");
        System.out.println(res50[0] + "\t\t" + res50[1] + "\t\t" + res50[2]);
        System.out.println("\n--- 200 mL");
        System.out.println("Phase I\t\tPhase II\tPhase III");
        System.out.println(res200[0] + "\t\t" + res200[1] + "\t\t" + res200[2]);
   
    }

    public static void original_mean_halftimes() {
        
        double[] params = {40,
            0.05,
            37,
            100};
        
        int N = 10000;
        DecimalFormat f = new DecimalFormat("###.0000");
        
        double Vo = 200;
        double[] t200_pI = thalf(params, Vo, 0, 60, N);
        double[] t200_pII = thalf(params, Vo, 60, 105, N);
        double[] t200_pIII = thalf(params, Vo, 105, 115, N);
        
        Vo = 50;
        double[] t50_pI = thalf(params, Vo, 0, 60, N);
        double[] t50_pII = thalf(params, Vo, 60, 105, N);
        double[] t50_pIII = thalf(params, Vo, 105, 115, N);
        
        String[] res50 = {f.format(StatUtils.geometricMean(t50_pI)/60),
            f.format(StatUtils.geometricMean(t50_pII)/60),
            f.format(StatUtils.geometricMean(t50_pIII)/60)};
        String[] res200 = {f.format(StatUtils.geometricMean(t200_pI)/60),
            f.format(StatUtils.geometricMean(t200_pII)/60),
            f.format(StatUtils.geometricMean(t200_pIII)/60)};

        
        System.out.println("\n--- Original Emptying Half-Times:");
        System.out.println("\n--- 50 mL");
        System.out.println("Phase I\t\tPhase II\t\tPhase III");
        System.out.println(res50[0] + "\t\t" + res50[1] + "\t\t" + res50[2]);
        System.out.println("\n--- 200 mL");
        System.out.println("Phase I\tPhase II\tPhase III");
        System.out.println(res200[0] + "\t\t" + res200[1] + "\t\t" + res200[2]);
    }
    
    public static double process(double[] params) {
        
        int N = 10000;
        
        double Vo = 200;
        double[] t200_pI = thalf(params, Vo, 0, 60, N);
        double[] t200_pII = thalf(params, Vo, 60, 105, N);
        double[] t200_pIII = thalf(params, Vo, 105, 115, N);
        
        Vo = 50;
        double[] t50_pI = thalf(params, Vo, 0, 60, N);
        double[] t50_pII = thalf(params, Vo, 60, 105, N);
        double[] t50_pIII = thalf(params, Vo, 105, 115, N);
        
        double err_pI200 = Math.pow(22.8*60-StatUtils.geometricMean(t200_pI),2);
        double err_pII200 = Math.pow(12.2*60-StatUtils.geometricMean(t200_pII),2);
        double err_pIII200 = Math.pow(4.92*60-StatUtils.geometricMean(t200_pIII),2);
        double err_200 = Math.pow(err_pI200+err_pII200+err_pIII200,0.5);
        double err_pI50 = Math.pow(60.6*60-StatUtils.geometricMean(t50_pI),2);
        double err_pII50 = Math.pow(17.4*60-StatUtils.geometricMean(t50_pII),2);
        double err_pIII50 = Math.pow(9.0*60-StatUtils.geometricMean(t50_pIII),2);
        double err_50 = Math.pow(err_pI50+err_pII50+err_pIII50,0.5);
        
        double tot_err = Math.pow(err_200+err_50,.5);
        
        return tot_err;
    }
    
    private static void optimize() throws Exception {
        
        int N = 10000;
        int i = 0;
        double new_err = 0;
        int opt_ind = 0;
        double opt_err = process(original_params);
        
        UniformRealDistribution dist_a = new UniformRealDistribution(1,50);
        UniformRealDistribution dist_b = new UniformRealDistribution();
        UniformRealDistribution dist_c = new UniformRealDistribution(1,50);
        UniformRealDistribution dist_tau = new UniformRealDistribution(0,120);
        
        double[] vals_a = new double[N];
        double[] vals_b = new double[N];
        double[] vals_c = new double[N];
        double[] vals_tau = new double[N];

        for (i = 0; i < N; i++) {
            vals_a[i] = dist_a.sample();
            vals_b[i] = dist_b.sample();
            vals_c[i] = dist_c.sample();
            vals_tau[i] = dist_tau.sample();
        }
        
        Date date = new Date();
        String start_time = date.toString();
        System.out.println("--- Starting ptimization [" + start_time + "]");
        
        for (i = 0; i < N; i++) {
            String anim= "|/-\\";

            double pct_complete = Math.round((float) (i)/N*100*10.0)/10.0;
            int index = i+1;
            String data = "\r" + "--- Calculating Results " + anim.charAt(index % anim.length()) + " "  + index + "/" + N + " [" + pct_complete + "%] Error: " + opt_err ;
            System.out.write(data.getBytes());
            
            double[] params = {vals_a[i],vals_b[i], vals_c[i], vals_tau[i]};

            new_err = process(params);
            if (new_err<opt_err) {
                opt_err = new_err;
                opt_ind = i;
            }
        }
        System.out.println("\n--- Optimization complete");
        System.out.println("--- Error value: " + opt_err + "\n");
        System.out.println("a\t\t\tb\t\t\tc\t\t\ttau");
        System.out.println(vals_a[opt_ind] + "\t" +
                           vals_b[opt_ind] + "\t" +
                           vals_c[opt_ind] + "\t" +
                           vals_tau[opt_ind]);
        
        double[] opt_params = {vals_a[opt_ind],vals_b[opt_ind],vals_c[opt_ind],vals_tau[opt_ind]};
        mean_halftimes(opt_params);
        
        original_mean_halftimes();
        
    }
    
    public static void saveResults(double[] X, String title) throws Exception {
        System.out.println("\n--- Writing results to file: " + title);
        FileWriter fw;
        try {
            fw = new FileWriter(new File(title));
            BufferedWriter bw = new BufferedWriter (fw);
            String tempStr = "";
            for (int i = 0; i < X.length; i++) {
                bw.write(String.valueOf(X[i]) + ",");
            }
            bw.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    public static void run() throws Exception,FileNotFoundException {
        
        int N = 5000;
        
        System.out.println("\n--- " + N + " iterations");
        System.out.println("\n--- 200 mL volume...");
        double Vo = 200;
        System.out.println("\t* Phase I");
        double[] t200_pI = thalf(original_params, Vo, 0, 60, N);
        System.out.println("\t* Phase II");
        double[] t200_pII = thalf(original_params, Vo, 60, 105, N);
        System.out.println("\t* Phase III");
        double[] t200_pIII = thalf(original_params, Vo, 105, 115, N);
        
        System.out.println("\n--- 50 mL volume...");
        Vo = 50;
        System.out.println("\t* Phase I");
        double[] t50_pI = thalf(original_params, Vo, 0, 60, N);
        System.out.println("\t* Phase II");
        double[] t50_pII = thalf(original_params, Vo, 60, 105, N);
        System.out.println("\t* Phase III");
        double[] t50_pIII = thalf(original_params, Vo, 105, 115, N);
        
        System.out.println("\n--- Saving Results...");
        saveResults(t200_pI,"results_t200_pI.txt");
        saveResults(t200_pII,"results_t200_pII.txt");
        saveResults(t200_pIII,"results_t200_pIII.txt");
        saveResults(t50_pI,"results_t50_pI.txt");
        saveResults(t50_pII,"results_t50_pII.txt");
        saveResults(t50_pIII,"results_t50_pIII.txt");
        
    }
    
    public static void main(String[] args) throws Exception,FileNotFoundException {
        
        switch (args[0]) {
            case "optimize" :
                optimize();
                
            default :
                run();
        }
    }
}
