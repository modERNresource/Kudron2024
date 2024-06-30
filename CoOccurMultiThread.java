/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rhwlab.chipseq.cluster;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import org.rhwlab.bedformat.BedBase;

/**
 *
 * @author gevirl
 */
public class CoOccurMultiThread implements Callable<Welford[]> {

    String[] clusters;
    String[] tfs;
    short[] peakTF;
    int[] peakCluster;
    TreeMap<String, Integer> peakCounts = new TreeMap<>();
    
    ArrayList<Callable<int[]>> workers = new ArrayList<>();
    int nWork;
    int nIter;
    Welford[] w ;

    public CoOccurMultiThread(File peakBedFile, int nIterPerWorker, int nWorkers) throws Exception {
        this.nWork = nWorkers;
        this.nIter = nIterPerWorker;

        TreeMap<String, Integer> clusterMap = new TreeMap<>();
        TreeMap<String, Short> tfMap = new TreeMap<>();
        ArrayList<Integer> peakClusterList = new ArrayList<>();
        ArrayList<Short> peakTFList = new ArrayList<>();

        BufferedReader reader = new BufferedReader(new FileReader(peakBedFile));
        String line = reader.readLine();
        while (line != null) {
            BedBase bed = new BedBase(line);
            String clusterName = bed.getValue(9);
            if (!clusterName.startsWith("Sing")) {
                Integer cIndex = clusterMap.get(clusterName);
                if (cIndex == null) {
                    cIndex = clusterMap.size();
                    clusterMap.put(clusterName, cIndex);
                }
                String tfName = bed.getName().split("_")[0];
                Short tfIndex = tfMap.get(tfName);
                if (tfIndex == null) {
                    tfIndex = (short) tfMap.size();
                    tfMap.put(tfName, tfIndex);
                }
                peakClusterList.add(cIndex);
                peakTFList.add(tfIndex);

                Integer count = peakCounts.get(tfName);
                if (count == null) {
                    peakCounts.put(tfName, 1);
                } else {
                    peakCounts.put(tfName, count + 1);
                }
            }
            line = reader.readLine();
        }

        peakTF = new short[peakTFList.size()];
        int i = 0;
        for (Short s : peakTFList) {
            peakTF[i] = s;
            ++i;
        }

        peakCluster = new int[peakClusterList.size()];
        i = 0;
        for (Integer c : peakClusterList) {
            peakCluster[i] = c;
            ++i;
        }
        reader.close();

        tfs = new String[tfMap.size()];
        for (Map.Entry<String, Short> e : tfMap.entrySet()) {
            tfs[e.getValue()] = e.getKey();
        }

        clusters = new String[clusterMap.size()];
        for (Map.Entry<String, Integer> e : clusterMap.entrySet()) {
            clusters[e.getValue()] = e.getKey();
        }
        
        // make the workers
        for (int nw =0 ; nw<nWork ; ++nw){
            workers.add(new CoOccurWorker(this));
        }   
        
        
        int nPairs = (tfs.length * (tfs.length - 1)) / 2;
        w = new Welford[nPairs];
        for (int n = 0; n < nPairs; ++n) {
            w[n] = new Welford();
        }        
    }

    @Override
    public Welford[] call() throws Exception {
        ExecutorService service = Executors.newWorkStealingPool();
        for (int iter = 0; iter < nIter; ++iter) {
            System.out.println(iter*nWork);
            List<Future<int[]>> futures = service.invokeAll(workers);
            for (Future<int[]> future : futures){
                int[] values = future.get();
//System.out.printf("values.length = %d   w.length = %d\n", values.length,w.length);
                for (int n=0 ; n<values.length ; ++n){
                    w[n].update( (double)values[n]);
                }
            }
        }
        return w;
    }
    
    public void report(PrintStream stream){
        CoOccurWorker worker = (CoOccurWorker)this.workers.get(0);
        int[] c = worker.notRandomCounts();
        
        int p = 0;
        for (int tf1 = 0; tf1 < tfs.length - 1; ++tf1) {
            for (int tf2 = tf1 + 1; tf2 < tfs.length; ++tf2) {
                double mu = w[p].getMean();
                double sd = Math.sqrt(w[p].getVariance());
                double z = ((double)c[p] - mu) / sd;
                stream.printf("%s\t%d\t%s\t%d\t%f\t%f\t%f\t%f\n", tfs[tf1],peakCounts.get(tfs[tf1]),
                        tfs[tf2],peakCounts.get(tfs[tf2]), mu, sd, (double)c[p], z);
                ++p;
            }
        }        
    }
    
    static public void main(String[] args)throws Exception {
        int nWorkers = Integer.valueOf(args[0]);
        int nIter = Integer.valueOf(args[1]);
        File peakFile = new File(args[2]);
        File outFile = new File(peakFile.getPath()+String.format("shuf.multi%d_%d.pairs",nWorkers,nIter));
        PrintStream stream = new PrintStream(outFile);
        CoOccurMultiThread co = new CoOccurMultiThread(peakFile, nIter, nWorkers);
        co.call();
        co.report(stream);
        stream.close();

    }
}
