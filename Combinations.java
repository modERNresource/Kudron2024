/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rhwlab.chipseq.cluster;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

/**
 *
 * @author gevirl
 */
public class Combinations {
    List<Object[]> combinations = new ArrayList<>();
    
    public Combinations(Collection items,int n){
        
        if (n == 1){
            for (Object item : items){
                Object[] obj = new Object[1];
                obj[0] = item;
                combinations.add(obj);
            }
        } else {
            List itemList = new ArrayList(items);
            for (int i=0 ; i<items.size()-1 ; ++i){
                
                List nextList = new ArrayList<>();
                for (int j=i+1 ; j<items.size() ; ++j){
                    nextList.add(itemList.get(j));
                }
                Combinations nextComb = new Combinations(nextList,n-1);
                for (Object[] objs : nextComb.combinations){
                    Object[] c = new Object[objs.length+1] ;
                    c[0] = itemList.get(i);
                    for (int j=0 ; j<objs.length ; ++j){
                        c[j+1] = objs[j];
                    }
                    combinations.add(c);
                }
            }
        }        
    }
    static public int[][] IntegerCombinatons(int n,int k){
        int size = 1;
        for (int i=0 ; i<k ; ++i){
            size = size*(n-i)/(i+1);
        }
        int[][] ret = new int[size][];
        for (int i=0 ; i<size ; ++i){
            ret[i] = new int[k];
        }
        IntegerCombinatons(0,n-1,k,ret);
        return ret;
    }
    static public void IntegerCombinatons(int start,int end,int k,int[][] comb){
        
    }
    

    static public void main(String[] args) {
        
        int[][] c = IntegerCombinatons(7,2);
        
        List<String> set = new ArrayList<>();
        set.add("One");
        set.add("Two");
        set.add("Three");
        set.add("Four");
        set.add("Five");
        set.add("Six");
        set.add("Seven");
        Combinations comb = new Combinations(set,3);
        int dhd=0;
    }
}
