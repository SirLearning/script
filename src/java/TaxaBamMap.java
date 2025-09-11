package pgl.LAW.tmp.vmap4;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.ParseException;
import pgl.AppAbstract;
import pgl.infra.utils.IOUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

public class TaxaBamMap extends AppAbstract {
    private String bamFiles = "";
    private String depthS;
    private String bamS;
    private String outFile;
    private String taxaRunFile;

    public TaxaBamMap(String[] args) {
        creatAppOptions();
        retrieveAppParameters(args);
        createTaxaBamMap();
    }

    public static void main(String[] args) {
        String[] wapArgs = new String[] {
                "-d", "/data/home/dazheng/vmap4/01ctDepth",
                "-b", "/data/home/dazheng/vmap4/00data/bam_ABD",
                "-o", "/data/home/dazheng/vmap4/WAP.taxaBamMap.txt"
        };
        String[] subArgs = new String[] {
                "-d", "/data/home/dazheng/vmap4/01ctDepth",
                "-b", "/data/home/dazheng/vmap4/02subBams",
                "-o", "/data/home/dazheng/vmap4/test.taxaBamMap.txt"
        };
        // new TaxaBamMap(wapArgs);
        // new TaxaBamMap(subArgs);
        new TaxaBamMap(args);
    }

    private static List<File> getPureBamFiles(List<File> input) {
        List<File> output = new ArrayList<>();
        for (File file : input) {
            if (file.getName().endsWith("bam")) output.add(file);
        }
        return output;
    }

    private String[] getBamFilesFromFile() throws IOException {
        BufferedReader reader = IOUtils.getTextReader(bamFiles);
        List<String> bamList = new ArrayList<>();
        String line;
        while ((line = reader.readLine()) != null && !line.isEmpty()) {
            bamList.add(line);
        }
        return bamList.toArray(new String[0]);
    }

    private static File findBamBasedOnDepth(List<File> bams, File depth) {
        String factor = depth.getName().split("\\.")[0];
        File returnBam = null;
        for (File bam: bams) {
            if (bam.getName().startsWith(factor)) returnBam = bam;
        }
        return returnBam;
    }

    @Override
    public void creatAppOptions() {
        options.addOption("f", true, "Bam files list file");
        options.addOption("d", true, "Input depth files data site");
        options.addOption("b", true, "Input bam files data site");
        options.addOption("o", true, "Output taxaBamMap.txt");
        options.addOption("t", true, "Output file of taxa and sequencing run number");
    }

    @Override
    public void retrieveAppParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            this.bamFiles = line.getOptionValue("f");
            this.depthS = line.getOptionValue("d");
            this.bamS = line.getOptionValue("b");
            this.outFile = line.getOptionValue("o");
            this.taxaRunFile = line.getOptionValue("t");
        } catch (ParseException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void printInstructionAndUsage() {

    }

    private void createTaxaBamMap() {
        String cmd = "samtools view -H ";
        List<File> depthList = IOUtils.getFileListInDirContains(depthS, "summary");
        List<File> bamList = List.of();
        try (
                BufferedWriter bw = new BufferedWriter(new FileWriter(outFile, true));
                BufferedWriter bw1 = new BufferedWriter(new FileWriter(taxaRunFile, true))
        ) {   // append content to the file before
            if (bamFiles == null || bamFiles.isEmpty()) {
                bamList = getPureBamFiles(IOUtils.getFileListInDirContains(bamS, "bam"));
            } else {
                String[] bamFiles = getBamFilesFromFile();
                for (String bamFile : bamFiles) {
                    bamList.add(new File(bamS, bamFile));
                }
            }
            File depth;
            File bam;
            ExecutorService pool = Executors.newFixedThreadPool(10);
            List<Future<BamHeader>> futureList = new ArrayList<>();
            for (File file : depthList) {
                depth = file;
                bam = findBamBasedOnDepth(bamList, depth);
                BufferedReader br = IOUtils.getTextReader(depth.getAbsolutePath());
                TaxonRead tr = new TaxonRead(cmd + bam.getAbsolutePath());
                Future<BamHeader> f = pool.submit(tr);
                BamHeader tmpBamHeader;
                futureList.add(f);
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.startsWith("total")) {
                        String[] elements = line.split("\t");
                        //System.out.println(elements[3]);
                        //System.out.println(f.get().taxa);
                        //System.out.println(f.get().taxa + "\t" + elements[3] + "\t" + bam.getAbsolutePath());
                        tmpBamHeader = f.get();
                        bw.write(tmpBamHeader.runNumber + "\t" + elements[3] + "\t" + bam.getAbsolutePath() + "\n");
                        bw1.write(tmpBamHeader.taxa + "\t" + tmpBamHeader.runNumber + "\n");
                    }
                }
                br.close();
                //pool.shutdownNow();
            }
            pool.close();
        } catch (IOException | ExecutionException | InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    static class TaxonRead implements Callable<BamHeader> {
        String command;
        public TaxonRead (String command) {
            this.command = command;
        }

        @Override
        public BamHeader call() throws Exception {
            BamHeader bh;
            String taxon = null;
            String run = null;
            try {
                Runtime rt = Runtime.getRuntime();
                Process p = rt.exec(command);
                System.out.println(command);
                String temp;
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                //System.out.println(br.readLine());
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("@RG")) {
                        String[] elements = temp.split("\t");
                        for (String element : elements) {
                            if (element.startsWith("SM")) {
                                taxon = element.split(":")[1];
                                break;
                            }
                        }
                    }
                    if (temp.startsWith("@PG")) {
                        String[] elements = temp.split("\t");
                        for (String element : elements) {
                            if (element.startsWith("CL") && element.endsWith(".gz")) {
                                run = temp.split("/")[temp.split("/").length - 1].split("_")[0];
                                break;
                            }
                        }
                    }
                }
                bh = new BamHeader(taxon, run);
                p.waitFor();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
            return bh;
        }
    }
}
