/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rhwlab.chipseq.cluster;

import java.util.Random;
import org.apache.commons.math3.stat.StatUtils;

/**
 *
 * @author gevirl
 */

// running calculation of  mean and variance 
public class Welford {
    double M2;
    int count;
    double mean;
    
    public Welford(){
        
    }
    
    public void update(double x){
        ++count;   
        double delta = x - mean;
        mean = mean + delta/count;
        double delta2 = x - mean;
        M2 = M2 + delta*delta2;
    }
    
    public double getMean(){
        return mean;
    }
    
    public double getVariance(){
        return M2/count;
    }
    
    static public void main(String[] args) {
        Random rnd = new Random();
        double[] x = new double[50000];
        for (int i=0 ; i<x.length ; ++i){
            x[i] = rnd.nextGaussian();
        }
        
        Welford welford = new Welford();
        for (int i=0 ; i<x.length ; ++i) {
            welford.update(x[i]);
        }
        System.out.printf("Welford Mean = %f\t var=%f\n", welford.getMean(),welford.getVariance());
        double mu = StatUtils.mean(x);
        double s = StatUtils.variance(x, mu);
        System.out.printf("Apache Mean = %f\t var=%f\n", mu,s);
        
    }
}
