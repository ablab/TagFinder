package ru.spbau.bioinf.palign;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;

public class EMassAdapter {

    private static final ProcessBuilder pb = new ProcessBuilder("emass.exe");
    private static Process process;
    private static  OutputStream processOutput;

    private static InputStream processInput;

    private static BufferedReader reader;
    private static BufferedWriter writer;

    static {
        process = null;
        try {
            process = pb.start();
            processOutput = process.getOutputStream();
            processInput = process.getInputStream();            
            reader = new BufferedReader(new InputStreamReader(processInput));
            writer = new BufferedWriter(new OutputStreamWriter(processOutput));           
        } catch (IOException e) {
            e.printStackTrace();  
        }
    }

    public static void main(String[] args) {
        ArrayList<double[]> d = getPeaks(new int[]{5,4,8,2,1});
        for (double[] doubles : d) {
            System.out.println(doubles[0] + " " + doubles[1]);
        }
    }

    public static ArrayList<double[]> getPeaks(int[] atomCount) {
        ArrayList<double[]> ans = new ArrayList<double[]>();        
        String formula = "C" + atomCount[0] + "H" + atomCount[1] + "N" + atomCount[2] + "O" + atomCount[3] + "S" + atomCount[4] + "\n";
        System.out.println("formula = " + formula);
        try {
            writer.write(formula);
            writer.flush();

            String s = null;
            while ((s = reader.readLine())!= null) {
                if (s.startsWith("formula")) {
                    continue;
                }
                if (s.startsWith("end")) {
                    break;
                }
                String[] tokens = s.split(" ");
                ans.add(new double[]{Double.parseDouble(tokens[0]), Double.parseDouble(tokens[1])});
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return ans;
    }
}
