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
import java.util.Map.Entry;
import java.util.TreeMap;

/**
 *
 * @author gevirl
 */
public class CoOccurFile {
    ArrayList<String[]> list = new ArrayList<>();
    TreeMap<String,Integer> tfIndex = new TreeMap<>();
    
    public CoOccurFile(File file)throws Exception {
        int i=0;
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line = reader.readLine();
        while (line != null) {
            String[] tokens = line.split("\t");
            list.add(tokens);    
            
            Integer index = tfIndex.get(tokens[0]);
            if (index == null){
                index = i;
                ++i;
                tfIndex.put(tokens[0], index);
            }
            index = tfIndex.get(tokens[2]);
            if (index == null){
                index = i;
                ++i;
                tfIndex.put(tokens[2], index);
            }            
            line = reader.readLine();
        }
        reader .close();
    }
    
    public void reportAsMatrix(PrintStream stream){
        double[][] mat = new double[tfIndex.size()][tfIndex.size()];
        for (String[] tokens : list){
            mat[tfIndex.get(tokens[0])][tfIndex.get(tokens[2])]= Double.valueOf(tokens[7]);
            mat[tfIndex.get(tokens[2])][tfIndex.get(tokens[0])]= Double.valueOf(tokens[7]);
        }
        TreeMap<Integer,String> reverseMap = new TreeMap<>();
        for (Entry<String,Integer> e : tfIndex.entrySet()){
            reverseMap.put(e.getValue(),e.getKey());
        }
        // prin the column headings
        boolean first = true;
        for (String tf : reverseMap.values()){
            if (first){
                stream.printf("%s", tf);
            } else {
                stream.printf("\t%s", tf);
            }
            first = false;
        }
        stream.println();
        
        for (int r=0 ; r<mat.length ; ++r){
            stream.printf("%f", mat[r][0]);
            for (int c=1 ; c<mat[r].length ; ++c){
                stream.printf("\t%f", mat[r][c]);
            }
            stream.println();
        }
        int iausdf=0;
    }
    static public void main(String[] args)throws Exception {
        File file = new File("/net/waterston/vol9/ChipSeqPipeline/test/AllWormPeaks.TF.noDups.clustered.bed.pairs");
        PrintStream stream = new PrintStream(file.getPath().replace(".pairs", ".matrix.pairs"));
        CoOccurFile co = new CoOccurFile(file);
        co.reportAsMatrix(stream);
        int iuasdf=0;
    }
}
