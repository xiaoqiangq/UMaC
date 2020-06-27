package com.libs.commandline;

import com.libs.core.ObjectSave;
import com.libs.core.Segment;
import com.libs.core.Extract;
import com.libs.core.Extract2;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.io.File;
import java.io.IOException;
import java.util.Properties;
import java.util.logging.Logger;

/**
 * Created by wei_qiang on 11/12/18.
 */
public class UMaC {

    private static final Logger logger = Logger.getLogger(org.apache.commons.cli.CommandLine.class.getName());

    public static void main(String[] args) throws IOException, ClassNotFoundException {
        String usage = "\nCommandLine\n";
        usage = usage + "\nUsage: java -jar *.jar <COMMAND>[OPTIONS] \n" +
                "java -Xmx4g -jar *.jar extract -I *.bam -O *.txt # -Ref * -Tumor * -TFx 0.05 -Cutoff 2  -Window *.format \n"+
                "java -Xmx4g -jar *.jar saveRef -I *LD.gz -O * # -Min 50 -Max 1000 \n" +
                "java -Xmx4g -jar *.jar saveTumor -I *LD.gz -O * # -Min 50 -Max 1000 \n" ;//+
        //"java -Xmx4g -jar *.jar segment -I *rate.gz -Window *.format -O *.wig \n" ;

        String cmd;
        if (args.length >0){
            if (args[0].equals("extract") || args[0].equals("saveRef") || args[0].equals("saveTumor") || args[0].equals("segment") || args[0].equals("extract2")) {
                cmd = args[0];
            }
            else {
                logger.severe("command is not recognized!\n" + usage);
                return;
            }
        }
        else {
            System.out.println(usage);
            return;
        }
        /* if checking is OK, will be run the run*/

        runCmd(cmd, args);
    }

    /** Run the cmd */
    private static void runCmd(String cmd, String[] args) throws IOException, ClassNotFoundException {

        /* get the start time of running*/
        //long start = System.currentTimeMillis();
        CommandLine commandLine = null;
        CommandLineParser parser = new DefaultParser();

        /* get the options using createOptions methods*/
        Options options = createOptions(cmd);
        try {
            if (options != null)
                commandLine = parser.parse(options, args);
        }
        catch (ParseException parseException) {
            logger.severe("Invalid command line parameters!");
        }

        /* checking whether the option is right*/
        Properties properties;
        if (cmd.equals("extract")) {
            if (isValidated(commandLine, "extract"))
            {
                properties = getProperties(commandLine, "extract");
                new Extract(properties).run();

                //new IchorCNAMethylation(properties).run();
            }
            else {
                printHelp(options, "extract");
            }
        } else if (cmd.equals("saveRef") || cmd.equals("saveTumor")) {
            if (isValidated(commandLine, cmd))
            {
                properties = getProperties(commandLine, cmd);
                new ObjectSave(properties, cmd).run();//
            }
            else {
                printHelp(options, cmd);
            }

        }else if (cmd.equals("segment")) {
            if (isValidated(commandLine, "segment"))
            {
                properties = getProperties(commandLine, "segment");
                new Segment(properties).run();//
            }
            else {
                printHelp(options, "segment");
            }

        }else if (cmd.equals("extract2"))
        {
            if (isValidated(commandLine, "extract2"))
            {
                properties = getProperties(commandLine, "extract2");
                new Extract2(properties).run();
            }
            else {
                printHelp(options, "extract2");
            }
        }
    }
    private static Properties getProperties(CommandLine line, String cmd)
    {
        Properties properties = new Properties();
        if (cmd.equals("extract")) {
            properties.put("input", line.getOptionValue("input"));
            properties.put("output", line.getOptionValue("output"));

            if (line.hasOption("tumorRate"))
            {
                properties.put("tumorRate", line.getOptionValue("tumorRate"));
            }

            if (line.hasOption("refMethyRate"))
            {
                properties.put("refMethyRate", line.getOptionValue("refMethyRate"));
            }

            if (line.hasOption("tumorMethyRate"))
            {
                properties.put("tumorMethyRate", line.getOptionValue("tumorMethyRate"));
            }

            if (line.hasOption("window"))
            {
                properties.put("window", line.getOptionValue("window"));
            }

            if (line.hasOption("cutoff"))
            {
                properties.put("cutoff", line.getOptionValue("cutoff"));
            }

        }
        else if (cmd.equals("saveRef") || cmd.equals("saveTumor"))
        {
            properties.put("input", line.getOptionValue("input"));
            properties.put("output", line.getOptionValue("output"));

            if (line.hasOption("min"))
            {
                properties.put("min", line.getOptionValue("min"));
            }

            if (line.hasOption("max"))
            {
                properties.put("max", line.getOptionValue("max"));
            }
        }
        else if (cmd.equals("segment"))
        {
            properties.put("input", line.getOptionValue("input"));
            properties.put("output", line.getOptionValue("output"));

            if (line.hasOption("window"))
            {
                properties.put("window", line.getOptionValue("window"));
            }
        }else if (cmd.equals("extract2")) {
            properties.put("input", line.getOptionValue("input"));
            properties.put("output", line.getOptionValue("output"));

            if (line.hasOption("tumorRate"))
            {
                properties.put("tumorRate", line.getOptionValue("tumorRate"));
            }

            if (line.hasOption("refMethyRate"))
            {
                properties.put("refMethyRate", line.getOptionValue("refMethyRate"));
            }

            if (line.hasOption("tumorMethyRate"))
            {
                properties.put("tumorMethyRate", line.getOptionValue("tumorMethyRate"));
            }

            if (line.hasOption("cutoff"))
            {
                properties.put("cutoff", line.getOptionValue("cutoff"));
            }

            if (line.hasOption("interval"))
            {
                properties.put("interval", line.getOptionValue("interval"));
            }
        }

        return properties;
    }

    /**create default options*/
    private static Options createOptions(String cmd) {
        Options options = new Options();
        if (cmd.equals("extract")) {

            Option input = Option.builder("I").argName("FILE").hasArg()
                    .desc("input a bam file (required)")
                    .longOpt("input")
                    .build();
            options.addOption(input);

            Option output = Option.builder("O").argName("FILE").hasArg()
                    .desc("output a bed file (required)")
                    .longOpt("output")
                    .build();
            options.addOption(output);

            Option tumorRate = Option.builder("TFx").argName("double").hasArg()
                    .desc("initial tumor fraction, default 0.05")
                    .longOpt("tumorRate")
                    .build();
            options.addOption(tumorRate);

            Option refMethyRate = Option.builder("Ref").argName("FILE").hasArg()
                    .desc("ref CpG sites' methylation rate, default ./resource/Normal")
                    .longOpt("refMethyRate")
                    .build();
            options.addOption(refMethyRate);

            Option tumorMethyRate = Option.builder("Tumor").argName("FILE").hasArg()
                    .desc("tumor CpG sites' methylation rate, default ./resource/Tumor")
                    .longOpt("tumorMethyRate")
                    .build();
            options.addOption(tumorMethyRate);

            Option window = Option.builder("Window").argName("File").hasArg()
                    .desc("bins across the whole genome, default ./resource/format.en")
                    .longOpt("window")
                    .build();
            options.addOption(window);

            Option cutoff = Option.builder("Cutoff").argName("int").hasArg()
                    .desc("cut off the number of reads' CpG sites, default 2")
                    .longOpt("cutoff")
                    .build();
            options.addOption(cutoff);

            return options;
        }
        else if (cmd.equals("saveRef") || cmd.equals("saveTumor")) {

            Option input = Option.builder("I").argName("FILE").hasArg()
                    .desc("input CpG.gz file (required)")
                    .longOpt("input")
                    .build();
            options.addOption(input);

            Option output = Option.builder("O").argName("FILE").hasArg()
                    .desc("output object model files (required)")
                    .longOpt("output")
                    .build();
            options.addOption(output);

            Option min = Option.builder("Min").argName("int").hasArg()
                    .desc("depth min for each CpG sites, default 30")
                    .longOpt("min")
                    .build();
            options.addOption(min);

            Option max = Option.builder("Max").argName("int").hasArg()
                    .desc("depth max for each CpG sites, default 1000")
                    .longOpt("max")
                    .build();
            options.addOption(max);

            return options;
        }
        else if (cmd.equals("segment")) {

            Option input = Option.builder("I").argName("FILE").hasArg()
                    .desc("input CpG.gz file (required)")
                    .longOpt("input")
                    .build();
            options.addOption(input);

            Option output = Option.builder("O").argName("FILE").hasArg()
                    .desc("output object files (required)")
                    .longOpt("output")
                    .build();
            options.addOption(output);

            Option window = Option.builder("Window").argName("File").hasArg()
                    .desc("bins across the whole genome, default ./resource/format.en")
                    .longOpt("window")
                    .build();
            options.addOption(window);

            Option cutoff = Option.builder("Cutoff").argName("int").hasArg()
                    .desc("cut off the number of reads' CpG sites")
                    .longOpt("cutoff")
                    .build();
            options.addOption(cutoff);

            return options;
        }else if (cmd.equals("extract2")) {

            Option input = Option.builder("I").argName("FILE").hasArg()
                    .desc("input bam file (required)")
                    .longOpt("input")
                    .build();
            options.addOption(input);

            Option output = Option.builder("O").argName("FILE").hasArg()
                    .desc("output details for all reads (required)")
                    .longOpt("output")
                    .build();
            options.addOption(output);

            Option tumorRate = Option.builder("TFx").argName("double").hasArg()
                    .desc("initial tumor fraction, default 0.05")
                    .longOpt("tumorRate")
                    .build();
            options.addOption(tumorRate);

            Option refMethyRate = Option.builder("Ref").argName("FILE").hasArg()
                    .desc("ref CpG sites' methylation rate")
                    .longOpt("refMethyRate")
                    .build();
            options.addOption(refMethyRate);

            Option tumorMethyRate = Option.builder("Tumor").argName("FILE").hasArg()
                    .desc("tumor CpG sites' methylation rate")
                    .longOpt("tumorMethyRate")
                    .build();
            options.addOption(tumorMethyRate);

            Option cutoff = Option.builder("Cutoff").argName("int").hasArg()
                    .desc("cut off the number of reads' CpG sites, default 0")
                    .longOpt("cutoff")
                    .build();
            options.addOption(cutoff);

            Option interval = Option.builder("L").argName("chr:start-end").hasArg()
                    .desc("interval")
                    .longOpt("interval")
                    .build();
            options.addOption(interval);

            return options;
        }

        return null;
    }

    /**validated the every option*/
    private static boolean isValidated(CommandLine line, String cmd) {
        boolean tag = true;

        if (cmd.equals("extract"))
        {
            if ((!line.hasOption("input")) || (!new File(line.getOptionValue("input")).isFile()))
            {
                logger.severe("The input file is not correctly specified!");
                tag = false;
            }

            if (!line.hasOption("output"))
            {
                logger.severe("The output file is not correctly specified!");
                tag = false;
            }

        }else if (cmd.equals("saveRef") || cmd.equals("saveTumor"))
        {
            if ((!line.hasOption("input")) || (!new File(line.getOptionValue("input")).isFile()))
            {
                logger.severe("The input file is not correctly specified!");
                tag = false;
            }

            if (!line.hasOption("output"))
            {
                logger.severe("The output file is not correctly specified!");
                tag = false;
            }
        }else if (cmd.equals("segment"))
        {
            if ((!line.hasOption("input")) || (!new File(line.getOptionValue("input")).isFile()))
            {
                logger.severe("The input file is not correctly specified!");
                tag = false;
            }

            if (!line.hasOption("output"))
            {
                logger.severe("The output file is not correctly specified!");
                tag = false;
            }

            if (!line.hasOption("window"))
            {
                logger.severe("The window file is not correctly specified!");
                tag = false;
            }
        }else if(cmd.equals("extract2"))
        {
            if ((!line.hasOption("input")) || (!new File(line.getOptionValue("input")).isFile()))
            {
                logger.severe("The input file is not correctly specified!");
                tag = false;
            }

            if (!line.hasOption("output"))
            {
                logger.severe("The output file is not correctly specified!");
                tag = false;
            }
        }

        return tag;
    }

    /**when something wrong, the print help will be displayed */
    private static void printHelp(Options options, String command) {
        System.out.println();
        String cmdLineSyntax = "java -jar *.jar " + command + " [OPTIONS]\n";
        HelpFormatter formatter = new HelpFormatter();
        //formatter.setOptionComparator(new OptionComarator());
        formatter.printHelp(cmdLineSyntax, options);
    }

    /* display the options by opts_order
     private static class OptionComarator<T extends Option> implements Comparator<T>
     {
     private static final String OPTS_ORDER = "rpnbktco";

     public int compare(T o1, T o2) {
     return OPTS_ORDER.indexOf(o1.getLongOpt().charAt(0)) - OPTS_ORDER.indexOf(o2.getLongOpt().charAt(0));
     }
     }
     */
}
