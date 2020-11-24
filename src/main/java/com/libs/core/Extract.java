package com.libs.core;

import com.libs.file.DeClassLoader;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;

/**
 * Created by wei_qiang on 11/12/18.
 */

public class Extract
{
    private final String output;
    private final File input;
    private double tumorRate = 0.05;

    private final String refMethyRateFile;
    private final String tumorMethyRateFile;

    private final File windows;
    private int cutoff = 2;

    private int referenceIndex = 0; //"chr1" or "1"

    public Extract (Properties properties)
    {
        this.input = new File(properties.getProperty("input"));
        this.output = properties.getProperty("output");

        IOUtil.assertFileIsReadable(input);

        if ( !(properties.getProperty("refMethyRate")== null)) {
            this.refMethyRateFile = properties.getProperty("refMethyRate");
        }else
        {
            String fileName = new File(System.getProperty("java.class.path")).getParent();
            refMethyRateFile= fileName+"/resource/Normal";
        }

        if ( !(properties.getProperty("tumorMethyRate")== null)) {
            this.tumorMethyRateFile = properties.getProperty("tumorMethyRate");
        }else
        {
            String fileName = new File(System.getProperty("java.class.path")).getParent();
            this.tumorMethyRateFile= fileName+"/resource/Tumor";
        }

        if ( !(properties.getProperty("tumorRate")== null)) {
            this.tumorRate = Double.parseDouble(properties.getProperty("tumorRate"));
        }

        if ( !(properties.getProperty("cutoff")== null)) {
            this.cutoff = Integer.parseInt(properties.getProperty("cutoff"));
        }

        if ( !(properties.getProperty("window")== null)) {
            this.windows = new File(properties.getProperty("window"));
        }else
        {
            String fileName = new File(System.getProperty("java.class.path")).getParent();
            this.windows= new File(fileName+"/resource/format.en");
        }
    }

    public void run() throws IOException, ClassNotFoundException
    {

// get all bins from window files
        List<List<String>> segments = readSegment(windows);
// analysis bam files

        String regex = "\\s+";

/*
 * estimate tumorRate TFx
 */

        ObjectInputStream istumor;
        HashMap<String,float[]> hashTumorMethyRate = null;

        if(tumorMethyRateFile.length() != 0)
        {
            //istumor = new ObjectInputStream(new FileInputStream(tumorMethyRateFile+".chr1"));
            DeClassLoader loaderTumor = new DeClassLoader();
            ByteArrayInputStream tumorByte = new ByteArrayInputStream(loaderTumor.file2Bytes(new File (tumorMethyRateFile+".chr1.en")));

            istumor = new ObjectInputStream(tumorByte);

            hashTumorMethyRate = (HashMap<String,float[]>) istumor.readObject();

            istumor.close();
        }
        //ObjectInputStream isref = new ObjectInputStream(new FileInputStream(refMethyRateFile+".chr1"));

        DeClassLoader loaderRef = new DeClassLoader();
        ByteArrayInputStream refByte = new ByteArrayInputStream(loaderRef.file2Bytes(new File (refMethyRateFile+".chr1.en")));
        ObjectInputStream isref = new ObjectInputStream(refByte);

        HashMap<String,float[]> hashNormalMethyRate;

        hashNormalMethyRate = (HashMap<String,float[]>) isref.readObject();
        isref.close();

///

        SamReader reader = SamReaderFactory.makeDefault().open(input);

        List<SAMSequenceRecord> sequences = reader.getFileHeader().getSequenceDictionary().getSequences();

        for(SAMSequenceRecord ssr : sequences)
        {
            if(ssr.getSequenceName().equals("chr1"))
            {
                referenceIndex=1;
            }
        }
/////////////////////////////////
        SAMRecordIterator iterator_tmp;
        if(referenceIndex ==1)
        {
            iterator_tmp = reader.query("chr1",0,0,false);
        }else
        {
            iterator_tmp = reader.query("1",0,0,false);
        }

        int xmflag = 0;

        if(iterator_tmp.hasNext())
        {
            String[] samRecord_split = iterator_tmp.next().getSAMString().split(regex);

            for (String samRecord_split2: samRecord_split)
            {
                if(samRecord_split2.startsWith("XM:Z"))
                {
                    break;
                }
                xmflag++;
            }
        }
        iterator_tmp.close();

        ////////////////////////////////////////////////////////////////
        int itera = 20;

        for (int i = 0; i < itera; i++)
        {
            double tumorReads = 0;
            double normalReads = 0;

            HashMap<String,SAMRecord>  hashPairEndSamRecord = new HashMap<String, SAMRecord>();

            SAMRecordIterator iterator;

            if(referenceIndex ==1)
            {
                iterator = reader.query("chr1",0,0,false);
            }else
            {
                iterator = reader.query("1",0,0,false);
            }

            List<Double> tumorNorPro;
            int flag;

            int startReadS;
            int startReadE;
            int endReadS;
            int endReadE;

            SAMRecord startRead;

            while (iterator.hasNext())
            {
                SAMRecord samRecord = iterator.next();
                flag = samRecord.getFlags();
                String readName = samRecord.getReadName();

                //readName = readName.split("_")[0];
                //readName = readName.split("_")[2];

                //String[] readNames = readName.split(":");
                //readName = readNames[readNames.length-4]+":"+readNames[readNames.length-3]+":"+readNames[readNames.length-2]+":"+readNames[readNames.length-1];

                String str;

                String chrom = samRecord.getReferenceName();
                    chrom = chrom.replace("chr","");

                startReadS = samRecord.getAlignmentStart();
                startReadE = samRecord.getAlignmentEnd();
                endReadS = startReadS;
                endReadE = startReadE;

                List list;

                if(flag == 99 || flag == 163 || flag == 147 || flag == 83)
                {
                    if (hashPairEndSamRecord.containsKey(readName))
                    {
                        startRead = hashPairEndSamRecord.get(readName);

                        startReadS = startRead.getAlignmentStart();
                        startReadE = startRead.getAlignmentEnd();

                        String[]  samRecord_split_start = startRead.getSAMString().split(regex);
                        String[] samRecord_split_end = samRecord.getSAMString().split(regex);

                        if ( startReadE - endReadS >=0) // overlap in pair end
                        {
                            if ( startReadE - endReadE >= 0 ) //-1
                            {
                                str = samRecord_split_start[xmflag].substring(5);
                            }
                            else if (startReadE - endReadS+1 < samRecord_split_end[xmflag].substring(5).length())
                            {
                                str = samRecord_split_start[xmflag].substring(5)+samRecord_split_end[xmflag].substring(5).substring(startReadE - endReadS+1);

                            }else
                            {
                                str = samRecord_split_start[xmflag].substring(5);
                            }

                        }else
                        {
                            String gaps = "";
                            int n = endReadS- startReadE-1;

                            if (n > 0)
                            {
                                char[] chars = new char[n];
                                Arrays.fill(chars, 'N');
                                gaps = new String(chars);
                            }

                            str = samRecord_split_start[xmflag].substring(5)+gaps+samRecord_split_end[xmflag].substring(5);
                        }

                        list = rateMethylation(str);

                        tumorNorPro = bayes(list,chrom,startReadS,hashNormalMethyRate,hashTumorMethyRate);

                        if (tumorNorPro.get(2).intValue() >= cutoff)
                        {
                            tumorReads += tumorNorPro.get(0);
                            normalReads += tumorNorPro.get(1);
                        }

                        hashPairEndSamRecord.remove(readName);
                    }
                    else
                    {
                        hashPairEndSamRecord.put(readName, samRecord);
                    }

                }
                else if (flag == 0 || flag == 16)
                {
                    startRead = samRecord;

                    startReadS = startRead.getAlignmentStart();
                    //startReadE = startRead.getAlignmentEnd();

                    String[]  samRecord_split_start = startRead.getSAMString().split(regex);

                    str = samRecord_split_start[xmflag].substring(5);

                    list = rateMethylation(str);

                    tumorNorPro = bayes(list,chrom,startReadS,hashNormalMethyRate,hashTumorMethyRate);

                    if (tumorNorPro.get(2).intValue() >= cutoff)
                    {
                        tumorReads += tumorNorPro.get(0);
                        normalReads += tumorNorPro.get(1);
                    }

                }
                //else {} // >1024
            }
            iterator.close();

            if ( (int) (tumorReads/(tumorReads+normalReads)*1000) > tumorRate*1000) // 0.001
            {
                tumorRate = tumorReads/(tumorReads+normalReads);
            }else {
                break;
            }
        }


        ////
        //DecimalFormat df = new DecimalFormat("0.000");

        //BlockCompressedOutputStream rateFile = new BlockCompressedOutputStream(output);
        //Writer writer = new OutputStreamWriter(rateFile, "UTF-8");
        //String line = "##tumorRate="+tumorRate+"\n";
        //line += "#ReadName"+"\t"+"Chrom"+"\t"+"StartReadS"+"\t"+"EndReadE"+"\tStrand\t"+"TumorPro"+"\t"+"NormalPro" +"\t"+"Length"+"\t"+"MethlyateRate"+"\t"+"allCpGSites"+"\t"+"TumorCpGSites"+"\t"+"NormalCpGSites"+"\t"+"MethylationStatus"+"\t"+"EndReadS";
        //writer.write(line+"\n");

        FileWriter fw = new FileWriter(output);
        BufferedWriter bw = new BufferedWriter(fw);

        //String[] chroms={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};
        String[] chroms={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"};

        int index = 0;

        for (String chrom1 : chroms)
        {
            index++;

            if (tumorMethyRateFile.length() != 0) {
                //istumor = new ObjectInputStream(new FileInputStream(tumorMethyRateFile+".chr"+chroms[i]));

                DeClassLoader loaderTumor = new DeClassLoader();
                ByteArrayInputStream tumorByte = new ByteArrayInputStream(loaderTumor.file2Bytes(new File(tumorMethyRateFile+".chr" + chrom1+".en")));

                istumor = new ObjectInputStream(tumorByte);

                hashTumorMethyRate = (HashMap<String, float[]>) istumor.readObject();
            }
            //ObjectInputStream isref = new ObjectInputStream(new FileInputStream(refMethyRateFile+".chr"+chroms[i]));

            loaderRef = new DeClassLoader();
            refByte = new ByteArrayInputStream(loaderRef.file2Bytes(new File(refMethyRateFile+ ".chr" + chrom1+".en")));
            isref = new ObjectInputStream(refByte);

            hashNormalMethyRate = (HashMap<String, float[]>) isref.readObject();
            isref.close();

            for(String segmenti : segments.get(index-1))
            {
                String[] str_split = segmenti.split(":");

                double tumorReads = 0.0;
                double normalReads = 0.0;

                HashMap<String, SAMRecord> hashPairEndSamRecord = new HashMap<String, SAMRecord>();

                SAMRecordIterator iterator;

                if(referenceIndex ==1)
                {
                    iterator = reader.query("chr"+chrom1, Integer.parseInt(str_split[1]), Integer.parseInt(str_split[2]), false);
                }else
                {
                    iterator = reader.query(chrom1, Integer.parseInt(str_split[1]), Integer.parseInt(str_split[2]), false);
                }

                List<Double> tumorNorPro;
                int flag;

                //double methlyateRate;
                //double allCpGSites;
                int startReadS;
                int startReadE;
                int endReadS;
                int endReadE;

                SAMRecord startRead;

                while (iterator.hasNext())
                {
                    SAMRecord samRecord = iterator.next();
                    flag = samRecord.getFlags();
                    String readName = samRecord.getReadName();

                    //readName = readName.split("_")[0];
                    //readName = readName.split("_")[2];

                    //String[] readNames = readName.split(":");
                    //readName = readNames[readNames.length - 4] + ":" + readNames[readNames.length - 3] + ":" + readNames[readNames.length - 2] + ":" + readNames[readNames.length - 1];

                    String str;

                    String chrom = samRecord.getReferenceName();
                        chrom = chrom.replace("chr","");

                    startReadS = samRecord.getAlignmentStart();
                    startReadE = samRecord.getAlignmentEnd();
                    endReadS = startReadS;
                    endReadE = startReadE;

                    List list;

                    if (flag == 99 || flag == 163 || flag == 147 || flag == 83) {
                        if (hashPairEndSamRecord.containsKey(readName)) {
                            startRead = hashPairEndSamRecord.get(readName);

                            startReadS = startRead.getAlignmentStart();
                            startReadE = startRead.getAlignmentEnd();

                            String[] samRecord_split_start = startRead.getSAMString().split(regex);
                            String[] samRecord_split_end = samRecord.getSAMString().split(regex);

                            if (startReadE - endReadS >= 0) // overlap in pair end
                            {
                                if (startReadE - endReadE >= 0) //-1
                                {
                                    str = samRecord_split_start[xmflag].substring(5);
                                } else if (startReadE - endReadS + 1 < samRecord_split_end[xmflag].substring(5).length()) {
                                    str = samRecord_split_start[xmflag].substring(5) + samRecord_split_end[xmflag].substring(5).substring(startReadE - endReadS + 1);

                                } else {
                                    str = samRecord_split_start[xmflag].substring(5);
                                }

                            } else {
                                String gaps = "";
                                int n = endReadS - startReadE - 1;

                                if (n > 0) {
                                    char[] chars = new char[n];
                                    Arrays.fill(chars, 'N');
                                    gaps = new String(chars);
                                }

                                str = samRecord_split_start[xmflag].substring(5) + gaps + samRecord_split_end[xmflag].substring(5);
                            }

                            list = rateMethylation(str);
                            //methlyateRate = (Double) list.get(1);
                            //allCpGSites = (Double) list.get(0);

                            tumorNorPro = bayes(list, chrom, startReadS, hashNormalMethyRate, hashTumorMethyRate);

                            if (tumorNorPro.get(2).intValue() >= cutoff) {
                                double tumortmp = tumorNorPro.get(0);
                                double normaltmp = tumorNorPro.get(1);

                                tumorReads += tumortmp;
                                normalReads += normaltmp;

                            }

                            hashPairEndSamRecord.remove(readName);
                        } else {
                            hashPairEndSamRecord.put(readName, samRecord);
                        }

                    } else if (flag == 0 || flag == 16) {

                        startRead = samRecord;
                        startReadS = startRead.getAlignmentStart();
                        //startReadE = startRead.getAlignmentEnd();

                        String[] samRecord_split_start = startRead.getSAMString().split(regex);

                        str = samRecord_split_start[xmflag].substring(5);

                        list = rateMethylation(str);
                        //methlyateRate = (Double) list.get(1);
                        //allCpGSites = (Double) list.get(0);

                        tumorNorPro = bayes(list, chrom, startReadS, hashNormalMethyRate, hashTumorMethyRate);

                        if (tumorNorPro.get(2).intValue() >= cutoff) {
                            double tumortmp = tumorNorPro.get(0);
                            double normaltmp = tumorNorPro.get(1);

                            tumorReads += tumortmp;
                            normalReads += normaltmp;

                        }

                    }
                    //else { } // >1024

                }
                iterator.close();


                if (tumorReads + normalReads == 0) {
                    bw.write(str_split[0] + "\t" + str_split[1] + "\t" + str_split[2] + "\t" + str_split[3] + "\t" + str_split[4] + "\t" + 0 + "\t" + 0 + "\t" + 0 + "\n");
                } else {
                    int d1 = (int) tumorReads;
                    int d2 = (int) normalReads;
                    double ratio = tumorReads / (tumorReads + normalReads);

                    bw.write(str_split[0] + "\t" + str_split[1] + "\t" + str_split[2] + "\t" + str_split[3] + "\t" + str_split[4] + "\t" + d1 + "\t" + d2 + "\t" + ratio + "\n");
                }
            }
        }

        bw.close();
        reader.close();

        //writer.close();
        //rateFile.close();
        //istumor.close();

        //new MakeTabixCompressedIndex(output).run();
    }


    public List<Double> bayes(List list, String chrom, int pos, HashMap<String,float[]> hashNormalMethyRate, HashMap<String,float[]> hashTumorMethyRate) throws IOException
    {
        //double allCpGSites = (Double) list.get(0);
        //double methlyateRate = (Double) list.get(1);
        //List<Integer> methyNum = (List<Integer>) list.get(2);
        //List<Integer> unmethyNum  = (List<Integer>) list.get(3);
        List<Integer> methyAndunmethy = (List<Integer>) list.get(4);
        List<Integer> methyAndunmethyInt  = (List<Integer>) list.get(5);

        double tumorData = 1.0;
        double normalData = 1.0;
        double total;

        double refCpGSites = 0;
        double tumorCpGSites = 0;


        for (int var =0; var < methyAndunmethyInt.size()-1; var++)
        {
            String tmp = methyAndunmethyInt.get(var).toString() + methyAndunmethyInt.get(var + 1).toString();

            int pos2 = pos + methyAndunmethy.get(var);
            //int posNext = pos + methyAndunmethy.get(var + 1);

            if (hashNormalMethyRate.containsKey(chrom + ":" + pos2) && hashTumorMethyRate.containsKey(chrom + ":" + pos2))
            {
                float[] normalMethyRate = hashNormalMethyRate.get(chrom+":"+pos2);

                if (tmp.equals("11"))
                {
                    if(normalMethyRate[0] > 0)
                    {
                        normalData *= normalMethyRate[0];
                    }else
                    {
                        normalData *= 0.05;
                    }

                } else if (tmp.equals("10"))
                {
                    if(normalMethyRate[1] > 0)
                    {
                        normalData *= normalMethyRate[1];
                    }else
                    {
                        normalData *= 0.05;
                    }

                } else if (tmp.equals("01"))
                {
                    if(normalMethyRate[2] > 0)
                    {
                        normalData *= normalMethyRate[2];
                    }else
                    {
                        normalData *= 0.05;
                    }

                } else if (tmp.equals("00"))
                {
                    if(normalMethyRate[3] > 0)
                    {
                        normalData *= normalMethyRate[3];
                    }else
                    {
                        normalData *= 0.05;
                    }
                }
                refCpGSites+=1.0;
                //}

                //if ()
                //{
                float[] tumorMethyRate = hashTumorMethyRate.get(chrom+":"+pos2);

                if (tmp.equals("11"))
                {
                    if(tumorMethyRate[0] > 0)
                    {
                        tumorData *= tumorMethyRate[0];
                    }else
                    {
                        tumorData *= 0.05;
                    }

                } else if (tmp.equals("10"))
                {
                    if(tumorMethyRate[1] > 0)
                    {
                        tumorData *= tumorMethyRate[1];
                    }else
                    {
                        tumorData *= 0.05;
                    }

                } else if (tmp.equals("01"))
                {
                    if(tumorMethyRate[2] > 0)
                    {
                        tumorData *= tumorMethyRate[2];
                    }else
                    {
                        tumorData *= 0.05;
                    }

                } else if (tmp.equals("00"))
                {
                    if(tumorMethyRate[3] > 0)
                    {
                        tumorData *= tumorMethyRate[3];
                    }else
                    {
                        tumorData *= 0.05;
                    }
                }

                tumorCpGSites+=1.0;
            }
        }

        tumorData *= tumorRate;
        normalData *= 1 - tumorRate;

        total = tumorData+normalData;

        tumorData = tumorData/total;
        normalData = normalData/total;

        List<Double> double1 = new ArrayList<Double>();

        double1.add(tumorData);
        double1.add(normalData);
        double1.add(refCpGSites);
        double1.add(tumorCpGSites);

        return double1;
    }

    private List rateMethylation(String str)
    {
        int methylatedNum = 0;
        int unmethylatedNum = 0;
        double allCpGSites;
        double methlyateRate= 10;

        List list = new ArrayList();
        List<Integer> methyNum = new ArrayList<Integer>();
        List<Integer> unmethyNum = new ArrayList<Integer>();
        List<Integer> methyAndunmethy = new ArrayList<Integer>();
        List<Integer> methyAndunmethyInt = new ArrayList<Integer>();

        for(int i=0;i<str.length();i++)
        {
            char ch = str.charAt(i);

            if (ch == 'Z')
            {
                methylatedNum++;
                methyNum.add(i);
                methyAndunmethy.add(i);
                methyAndunmethyInt.add(1);

            }else if (ch == 'z')
            {
                unmethylatedNum++;
                unmethyNum.add(i);
                methyAndunmethy.add(i);
                methyAndunmethyInt.add(0);
            }
        }

        allCpGSites = methylatedNum + unmethylatedNum;

        if (allCpGSites >0)
        {
            methlyateRate = methylatedNum/allCpGSites;
        }

        list.add(allCpGSites);
        list.add(methlyateRate);
        list.add(methyNum);
        list.add(unmethyNum);
        list.add(methyAndunmethy);
        list.add(methyAndunmethyInt);

        return list;
    }


    private List<List<String>> readSegment(File windows) throws IOException
    {
        DeClassLoader loaderWindows = new DeClassLoader();
        ByteArrayInputStream windowsByte = new ByteArrayInputStream(loaderWindows.file2Bytes(windows));

        BufferedReader bfReader = new BufferedReader(new InputStreamReader(windowsByte));

        String window;

        List<List<String>> segments = new ArrayList<List<String>>();

        List<String> segment = null;

        String[] str_split;
        String chromTemp = "";

        while((window = bfReader.readLine()) != null)
        {
            if (window.startsWith("#") || window.startsWith("fixedStep")) {
                continue;
            }
            str_split = window.split("\t");

            String chrom = str_split[0];
            int start = Integer.parseInt(str_split[1]);
            int end = Integer.parseInt(str_split[2]);

            if (chrom.equals("X"))
            {
                break;
            }

            if (segment == null)
            {
                segment = new ArrayList<String>();
                segment.add(chrom+":"+start+":"+end+":"+str_split[3]+":"+str_split[4]);
                chromTemp = chrom;

            }else if (chromTemp.equals(chrom))
            {
                segment.add(chrom+":"+start+":"+end+":"+str_split[3]+":"+str_split[4]);
            }else
            {
                segments.add(segment);

                segment = new ArrayList<String>();
                segment.add(chrom+":"+start+":"+end+":"+str_split[3]+":"+str_split[4]);
                chromTemp = chrom;
            }
        }

        segments.add(segment);

        bfReader.close();
        windowsByte.close();


        return segments;
    }

}
