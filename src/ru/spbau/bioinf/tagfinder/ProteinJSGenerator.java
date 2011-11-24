package ru.spbau.bioinf.tagfinder;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;
import ru.spbau.bioinf.tagfinder.util.ReaderUtil;

public class ProteinJSGenerator {
    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        List<Protein> proteins = conf.getProteins();
        PrintWriter out = ReaderUtil.createOutputFile(new File("proteins.js"));
        out.println("var proteins = [");
        for (Protein protein : proteins) {
            out.println("'" + protein.getSequence()+ "',");
        }
        out.println("];");
        out.close();

    }
}
