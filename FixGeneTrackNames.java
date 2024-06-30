/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rhwlab.chipseq.cluster;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.TreeSet;

/**
 *
 * @author gevirl
 */
public class FixGeneTrackNames {
    static public void main(String[] args)throws Exception {
        TreeSet<String> labels_0 = new TreeSet<>();
        TreeSet<String> labels_11 = new TreeSet<>();
        File file = new File("/net/waterston/vol9/WS285/c_elegans.PRJNA13758.WS285.annotations.wormbase.genePred");
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line = reader.readLine();
        while (line != null){
            String[] tokens = line.split("\t");
            labels_0.add(tokens[0].split(":")[0]);
            labels_11.add(tokens[11].split(":")[0]);
            line = reader.readLine();
        }
        reader.close();
    }
}
