package ru.spbau.bioinf.tagfinder;


public class AllDataGenerator {
    public static void main(String[] args) throws Exception {
        //ValidTags.main(args);
        System.out.println("\\documentclass{article}[12pt]\n" +
                "\n" +
                "\\usepackage{amsmath}\n" +
                "\\usepackage{amsthm}\n" +
                "\\usepackage{booktabs}\n" +
                "\\usepackage{epsfig}\n" +
                "\\usepackage{graphicx}\n" +
                "\\usepackage{lscape}\n" +
                "\\usepackage{multirow}\n" +
                "\\usepackage{natbib}\n" +
                "\\usepackage{setspace}\n" +
                "\\usepackage{epstopdf}\n" +
                "\n" +
                "\\def\\STbar{{\\overline{\\mathrm{ST}}}}\n" +
                "\\def\\STtilde{{\\widetilde{\\mathrm{ST}}}}\n" +
                "\n" +
                "\\setlength{\\topmargin}{0in}\n" +
                "\\setlength{\\headheight}{12pt}\n" +
                "\\setlength{\\headsep}{0.3in}\n" +
                "\\setlength{\\textheight}{8.7in}\n" +
                "\\setlength{\\oddsidemargin}{0in}\n" +
                "\\setlength{\\evensidemargin}{0in}\n" +
                "\\setlength{\\textwidth}{6.5in}\n" +
                "\n" +
                "\\title{Peptide Sequence Tags for Top-Down Spectra}\n" +
                "\\author{}\n" +
                "%\\date{}\n" +
                "\n" +
                "\\begin{document}\n" +
                "\n" +
                "\\maketitle\n" +
                "\n" +
                "\\doublespacing\n" +
                "\\begin{abstract}\n" +
                "\n" +
                "\n" +
                "\\end{abstract}" +
                "\n" +
                "%TEXT" +
                "\n\n" +
                "% MONO-TAGS");
        CalculateRelation.tableMono();
        TexTableGenerator.tableMonoCorrect();

        System.out.println("\n" +
                "% SPECTRUM GRAPHS\n");

        TexTableGenerator.tableCorrectAndProper();

        System.out.println("\n");

        CalculateRelation.tableCorrectVsProper();

        System.out.println();

        CalculateRelation.tableLongestCorrect();

        System.out.println();


        KdTableGenerator.printTablesProper(1);
        KdTableGenerator.printTablesCorrect(1);

        System.out.println();



        TexTableGenerator.tablesCorrectErrAndAdv();

        System.out.println("\n\\begin{landscape}\n");
        System.out.println("\\begin{table}[ht]\\footnotesize\n" +
                "\\vspace{3mm}\n" +
                "{\\centering\n" +
                "\\begin{center}\n" +
                "\\begin{tabular}{|c|c|}\n" +
                "  \\hline\n" +
                "  a pair of correct aa-s & an incorrect aa \\\\\n" +
                "  \\hline\n" +
                "  AG or GA & K \\\\\n" +
                "  \\hline\n" +
                "  AG or GA & Q \\\\\n" +
                "  \\hline\n" +
                "  AD or DA & \\multirow{3}{*}{W} \\\\\n" +
                "  EG or GE & \\\\\n" +
                "  SV or VS & \\\\\n" +
                "  \\hline\n" +
                "  GV or VG & R \\\\\n" +
                "  \\hline\n" +
                "  GG & N \\\\\n" +
                "  \\hline\n" +
                "\\end{tabular}\n" +
                "\\end{center}\n" +
                "\\par}\n" +
                "\\centering\n" +
                "\\caption{Pairs of consecutive amino acids that can be mistaken for a single one.}\n" +
                "\\vspace{3mm}\n" +
                "\\label{table:errors-vs}\n" +
                "\\end{table}\n");
        System.out.println("\\end{landscape}");

        System.out.println("\n\n% TOP-SCORING TAGS\n");


        IntencityTableGenerator.tableTopscoreAllAndAvg();

        System.out.println("\n\n\n% UNIDENTIFIED SPECTRA\n");

        System.out.println("\\begin{landscape}\n");

        UnmatchedStatistics.tableUnident();
        System.out.println();
        System.out.println();
        GenerateMatchesTable.tableMatches();



        System.out.println("\n\n\n%NEIGHBORS TABLE\n\n\n");
        System.out.println("\\end{landscape}");

        ///CalculateRelation.main(args);
        //IntencityTableGenerator.main(args);
        //UnmatchedStatistics.main(args);

        //KdTableGenerator.main(args);
        System.out.println("\n" +
                "\\end{document}");
    }
}
