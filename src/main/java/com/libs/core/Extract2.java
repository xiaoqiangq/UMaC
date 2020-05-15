package com.libs.core;

import com.libs.file.DeClassLoader;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.IOUtil;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;

/**
 * Created by wei_qiang on 11/12/18.
 */

public class Extract2
{
    private final String output;
    private final File input;
    private double tumorRate = 0.05;
    private boolean tumorRateTrainFlag = true;

    private String interval = "";
    private String refMethyRateFile = "";
    private String tumorMethyRateFile = "";

    private int cutoff = 0;

    private int referenceIndex = 0; //"chr1" or "1"


    public Extract2(Properties properties)
    {
        this.input = new File(properties.getProperty("input"));
        this.output = properties.getProperty("output");

        IOUtil.assertFileIsReadable(input);

        if ( !(properties.getProperty("interval")== null)) {
            this.interval = properties.getProperty("interval");
        }

        if ( !(properties.getProperty("refMethyRate")== null)) {
            this.refMethyRateFile = properties.getProperty("refMethyRate");
        }

        if ( !(properties.getProperty("tumorMethyRate")== null)) {
            this.tumorMethyRateFile = properties.getProperty("tumorMethyRate");
        }

        if ( !(properties.getProperty("tumorRate")== null)) {
            this.tumorRate = Double.parseDouble(properties.getProperty("tumorRate"));
            this.tumorRateTrainFlag = false;
        }

        if ( !(properties.getProperty("cutoff")== null)) {
            this.cutoff = Integer.parseInt(properties.getProperty("cutoff"));
        }
    }

    public void run() throws IOException
    {
        String regex = "\\s+";
        /*
         * estimate tumorRate
         * **/

        ObjectInputStream istumor = null;
        if(tumorMethyRateFile.length() != 0)
        {
            //istumor = new ObjectInputStream(new FileInputStream(tumorMethyRateFile+".chr1"));
            DeClassLoader loaderTumor = new DeClassLoader();
            ByteArrayInputStream tumorByte = new ByteArrayInputStream(loaderTumor.file2Bytes(new File (tumorMethyRateFile+".chr1.en")));

            istumor = new ObjectInputStream(tumorByte);
        }
        //ObjectInputStream isref = new ObjectInputStream(new FileInputStream(refMethyRateFile+".chr1"));

        DeClassLoader loaderRef = new DeClassLoader();
        ByteArrayInputStream refByte = new ByteArrayInputStream(loaderRef.file2Bytes(new File (refMethyRateFile+".chr1.en")));
        ObjectInputStream isref = new ObjectInputStream(refByte);

        HashMap<String,float[]> hashNormalMethyRate = null;
        HashMap<String,float[]> hashTumorMethyRate = null;

        try {
            hashNormalMethyRate = (HashMap<String,float[]>) isref.readObject();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
        isref.close();

        if(tumorMethyRateFile.length() != 0)
        {
            try {
                hashTumorMethyRate = (HashMap<String,float[]>) istumor.readObject();
            } catch (ClassNotFoundException e) {
                e.printStackTrace();
            }
            istumor.close();
        }


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

////////////////////
        ///////////////////////////
        if(tumorRateTrainFlag)
        {
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
                    //  //readName = readName.split("_")[2];

                    // String[] readNames = readName.split(":");
                    // readName = readNames[readNames.length-4]+":"+readNames[readNames.length-3]+":"+readNames[readNames.length-2]+":"+readNames[readNames.length-1];

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

                            if (tumorNorPro.get(2).intValue() >= 2)
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

                        if (tumorNorPro.get(2).intValue() >= 2)
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
                    //System.out.println(tumorRate);
                }else {
                    break;
                }
            }
        }


        ///
        //final SAMFileWriter writerBam = new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(reader.getFileHeader(), false, new File(output+".bam"));

        ////
        DecimalFormat df = new DecimalFormat("0.000");

        BlockCompressedOutputStream rateFile = new BlockCompressedOutputStream(output);
        Writer writer = new OutputStreamWriter(rateFile, "UTF-8");
        String line = "##tumorRate="+tumorRate+"\n";
        line += "#ReadName"+"\t"+"Chrom"+"\t"+"StartReadS"+"\t"+"EndReadE"+"\tStrand\t"+"TumorPro"+"\t"+"NormalPro" +"\t"+"Length"+"\t"+"MethlyateRate"+"\t"+"allCpGSites"+"\t"+"TumorCpGSites"+"\t"
                +"NormalCpGSites"+"\t"+"MethylationStatus"+"\t"+"EndReadS"+"\t"+"logLike1"+"\t"+"logLike2";
        writer.write(line+"\n");

        //String[] chroms={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};
        String[] chroms={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"};

        if (interval.equals(""))
        {
            for (String chrom1 : chroms)
            {
                istumor = null;
                if (tumorMethyRateFile.length() != 0) {
                    //istumor = new ObjectInputStream(new FileInputStream(tumorMethyRateFile+".chr"+chroms[i]));

                    DeClassLoader loaderTumor = new DeClassLoader();
                    ByteArrayInputStream tumorByte = new ByteArrayInputStream(loaderTumor.file2Bytes(new File(tumorMethyRateFile+".chr" + chrom1+".en")));

                    istumor = new ObjectInputStream(tumorByte);

                }
                //ObjectInputStream isref = new ObjectInputStream(new FileInputStream(refMethyRateFile+".chr"+chroms[i]));

                loaderRef = new DeClassLoader();
                refByte = new ByteArrayInputStream(loaderRef.file2Bytes(new File(refMethyRateFile+ ".chr" + chrom1+".en")));
                isref = new ObjectInputStream(refByte);

                hashNormalMethyRate = null;
                hashTumorMethyRate = null;

                try {
                    hashNormalMethyRate = (HashMap<String, float[]>) isref.readObject();
                } catch (ClassNotFoundException e) {
                    e.printStackTrace();
                }
                isref.close();

                if (tumorMethyRateFile.length() != 0) {
                    try {
                        hashTumorMethyRate = (HashMap<String, float[]>) istumor.readObject();
                    } catch (ClassNotFoundException e) {
                        e.printStackTrace();
                    }
                    istumor.close();
                }

                HashMap<String, SAMRecord> hashPairEndSamRecord = new HashMap<String, SAMRecord>();

                SAMRecordIterator iterator;

                if(referenceIndex ==1)
                {
                    iterator = reader.query("chr"+chrom1,0,0,false);
                }else
                {
                    iterator = reader.query(chrom1,0,0,false);
                }

                List<Double> tumorNorPro;
                int flag;

                double methlyateRate;
                double allCpGSites;
                int startReadS;
                int startReadE;
                int endReadS;
                int endReadE;

                SAMRecord startRead;

                while (iterator.hasNext()) {
                    SAMRecord samRecord = iterator.next();
                    flag = samRecord.getFlags();
                    String readName = samRecord.getReadName();

                    //readName = readName.split("_")[0];
                    //  //readName = readName.split("_")[2];

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

                    if (flag == 99 || flag == 163 || flag == 147 || flag == 83)
                    {
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
                            methlyateRate = (Double) list.get(1);
                            allCpGSites = (Double) list.get(0);

                            tumorNorPro = bayes(list, chrom, startReadS, hashNormalMethyRate, hashTumorMethyRate);

                            if (flag == 147 || flag == 99) {
                                line = readName + "\t" + chrom + "\t" + startReadS + "\t" + endReadE + "\t+\t" + tumorNorPro.get(0) + "\t" + tumorNorPro.get(1) + "\t" + str.length()
                                        + "\t" + df.format(methlyateRate) + "\t" + (int) allCpGSites + "\t" + tumorNorPro.get(2).intValue() + "\t" + tumorNorPro.get(3).intValue() + "\t" + str + "\t" + endReadS
                                        +"\t" +tumorNorPro.get(4) +"\t" +tumorNorPro.get(5); //+"\t"+hashPairEndSamRecord.size();
                            } else
                            {
                                line = readName + "\t" + chrom + "\t" + startReadS + "\t" + endReadE + "\t-\t" + tumorNorPro.get(0) + "\t" + tumorNorPro.get(1) + "\t" + str.length()
                                        + "\t" + df.format(methlyateRate) + "\t" + (int) allCpGSites + "\t" + tumorNorPro.get(2).intValue() + "\t" + tumorNorPro.get(3).intValue() + "\t" + str + "\t" + endReadS
                                        +"\t" +tumorNorPro.get(4) +"\t" +tumorNorPro.get(5); //+"\t"+hashPairEndSamRecord.size();
                            }

                            if (tumorNorPro.get(2).intValue() >= cutoff)
                            {
                                writer.write(line + "\n");

                                //writerBam.addAlignment(hashPairEndSamRecord.get(readName));
                                //writerBam.addAlignment(samRecord);
                            }

                            hashPairEndSamRecord.remove(readName);
                        } else {
                            hashPairEndSamRecord.put(readName, samRecord);
                        }

                    } else if (flag == 0 || flag == 16)
                    {

                        startRead = samRecord;
                        startReadS = startRead.getAlignmentStart();
                        startReadE = startRead.getAlignmentEnd();

                        String[] samRecord_split_start = startRead.getSAMString().split(regex);

                        str = samRecord_split_start[xmflag].substring(5);

                        list = rateMethylation(str);
                        methlyateRate = (Double) list.get(1);
                        allCpGSites = (Double) list.get(0);

                        tumorNorPro = bayes(list, chrom, startReadS, hashNormalMethyRate, hashTumorMethyRate);

                        if (flag == 0) {
                            line = readName + "\t" + chrom + "\t" + startReadS + "\t" + startReadE + "\t+\t" + tumorNorPro.get(0) + "\t" + tumorNorPro.get(1) + "\t" + str.length()
                                    + "\t" + df.format(methlyateRate) + "\t" + (int) allCpGSites + "\t" + tumorNorPro.get(2).intValue() + "\t" + tumorNorPro.get(3).intValue() + "\t" + str + "\t" + startReadE
                                    +"\t" +tumorNorPro.get(4) +"\t" +tumorNorPro.get(5); //+"\t"+hashPairEndSamRecord.size();
                        } else
                        {
                            line = readName + "\t" + chrom + "\t" + startReadS + "\t" + startReadE + "\t-\t" + tumorNorPro.get(0) + "\t" + tumorNorPro.get(1) + "\t" + str.length()
                                    + "\t" + df.format(methlyateRate) + "\t" + (int) allCpGSites + "\t" + tumorNorPro.get(2).intValue() + "\t" + tumorNorPro.get(3).intValue() + "\t" + str + "\t" + startReadE
                                    +"\t" +tumorNorPro.get(4) +"\t" +tumorNorPro.get(5); //+"\t"+hashPairEndSamRecord.size();
                        }

                        if (tumorNorPro.get(2).intValue() >= cutoff)
                        {
                            writer.write(line + "\n");
                            //writerBam.addAlignment(samRecord);
                        }
                    }
                    //else { } // >1024
                }
                iterator.close();

            }

        }else
        {
            String[] str_split = interval.split("[:\\-]");

            String chrom1 = str_split[0];
            int start = Integer.parseInt(str_split[1]);
            int end = Integer.parseInt(str_split[2]);

            istumor = null;
            if (tumorMethyRateFile.length() != 0) {
                //istumor = new ObjectInputStream(new FileInputStream(tumorMethyRateFile+".chr"+chroms[i]));

                DeClassLoader loaderTumor = new DeClassLoader();
                ByteArrayInputStream tumorByte = new ByteArrayInputStream(loaderTumor.file2Bytes(new File(tumorMethyRateFile+".chr" + chrom1+".en")));

                istumor = new ObjectInputStream(tumorByte);

            }
            //ObjectInputStream isref = new ObjectInputStream(new FileInputStream(refMethyRateFile+".chr"+chroms[i]));

            loaderRef = new DeClassLoader();
            refByte = new ByteArrayInputStream(loaderRef.file2Bytes(new File(refMethyRateFile+ ".chr" + chrom1+".en")));
            isref = new ObjectInputStream(refByte);

            hashNormalMethyRate = null;
            hashTumorMethyRate = null;

            try {
                hashNormalMethyRate = (HashMap<String, float[]>) isref.readObject();
            } catch (ClassNotFoundException e) {
                e.printStackTrace();
            }
            isref.close();

            if (tumorMethyRateFile.length() != 0) {
                try {
                    hashTumorMethyRate = (HashMap<String, float[]>) istumor.readObject();
                } catch (ClassNotFoundException e) {
                    e.printStackTrace();
                }
                istumor.close();
            }

            HashMap<String, SAMRecord> hashPairEndSamRecord = new HashMap<String, SAMRecord>();

            SAMRecordIterator iterator;

            //if(referenceIndex ==1)
            //{
             //   iterator = reader.query("chr"+chrom1,start,end,false);
            //}else
            //{
                iterator = reader.query(chrom1,start,end,false);
            //}

            List<Double> tumorNorPro;
            int flag;

            double methlyateRate;
            double allCpGSites;
            int startReadS;
            int startReadE;
            int endReadS;
            int endReadE;

            SAMRecord startRead;

            while (iterator.hasNext()) {
                SAMRecord samRecord = iterator.next();
                flag = samRecord.getFlags();
                String readName = samRecord.getReadName();

                //readName = readName.split("_")[0];
                // //readName = readName.split("_")[2];

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

                if (flag == 99 || flag == 163 || flag == 147 || flag == 83)
                {
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
                        methlyateRate = (Double) list.get(1);
                        allCpGSites = (Double) list.get(0);

                        tumorNorPro = bayes(list, chrom, startReadS, hashNormalMethyRate, hashTumorMethyRate);

                        if (flag == 147 || flag == 99) {
                            line = readName + "\t" + chrom + "\t" + startReadS + "\t" + endReadE + "\t+\t" + tumorNorPro.get(0) + "\t" + tumorNorPro.get(1) + "\t" + str.length()
                                    + "\t" + df.format(methlyateRate) + "\t" + (int) allCpGSites + "\t" + tumorNorPro.get(2).intValue() + "\t" + tumorNorPro.get(3).intValue() + "\t" + str + "\t" + endReadS
                                    +"\t" +tumorNorPro.get(4) +"\t" +tumorNorPro.get(5); //+"\t"+hashPairEndSamRecord.size();
                        } else
                        {
                            line = readName + "\t" + chrom + "\t" + startReadS + "\t" + endReadE + "\t-\t" + tumorNorPro.get(0) + "\t" + tumorNorPro.get(1) + "\t" + str.length()
                                    + "\t" + df.format(methlyateRate) + "\t" + (int) allCpGSites + "\t" + tumorNorPro.get(2).intValue() + "\t" + tumorNorPro.get(3).intValue() + "\t" + str + "\t" + endReadS
                                    +"\t" +tumorNorPro.get(4) +"\t" +tumorNorPro.get(5); //+"\t"+hashPairEndSamRecord.size();
                        }

                        if (tumorNorPro.get(2).intValue() >= cutoff)
                        {
                            writer.write(line + "\n");
                            //writerBam.addAlignment(hashPairEndSamRecord.get(readName));
                            //writerBam.addAlignment(samRecord);
                        }

                        hashPairEndSamRecord.remove(readName);
                    } else {
                        hashPairEndSamRecord.put(readName, samRecord);
                    }

                } else if (flag == 0 || flag == 16)
                {

                    startRead = samRecord;
                    startReadS = startRead.getAlignmentStart();
                    startReadE = startRead.getAlignmentEnd();

                    String[] samRecord_split_start = startRead.getSAMString().split(regex);

                    str = samRecord_split_start[xmflag].substring(5);

                    list = rateMethylation(str);
                    methlyateRate = (Double) list.get(1);
                    allCpGSites = (Double) list.get(0);

                    tumorNorPro = bayes(list, chrom, startReadS, hashNormalMethyRate, hashTumorMethyRate);

                    if (flag == 0) {
                        line = readName + "\t" + chrom + "\t" + startReadS + "\t" + startReadE + "\t+\t" + tumorNorPro.get(0) + "\t" + tumorNorPro.get(1) + "\t" + str.length()
                                + "\t" + df.format(methlyateRate) + "\t" + (int) allCpGSites + "\t" + tumorNorPro.get(2).intValue() + "\t" + tumorNorPro.get(3).intValue() + "\t" + str + "\t" + startReadE
                                +"\t" +tumorNorPro.get(4) +"\t" +tumorNorPro.get(5) ; //+"\t"+hashPairEndSamRecord.size();
                    } else
                    {
                        line = readName + "\t" + chrom + "\t" + startReadS + "\t" + startReadE + "\t-\t" + tumorNorPro.get(0) + "\t" + tumorNorPro.get(1) + "\t" + str.length()
                                + "\t" + df.format(methlyateRate) + "\t" + (int) allCpGSites + "\t" + tumorNorPro.get(2).intValue() + "\t" + tumorNorPro.get(3).intValue() + "\t" + str + "\t" + startReadE
                                +"\t" +tumorNorPro.get(4) +"\t" +tumorNorPro.get(5); //+"\t"+hashPairEndSamRecord.size();
                    }

                    if (tumorNorPro.get(2).intValue() >= cutoff)
                    {
                        writer.write(line + "\n");
                        //writerBam.addAlignment(samRecord);
                    }
                }
                //else { } // >1024
            }
            iterator.close();

        }

        reader.close();
        writer.close();
        //writerBam.close();

        //rateFile.close();
        //istumor.close();

        new MakeTabixCompressedIndex(output).run();

    }


    public List<Double> bayes(List list, String chrom, int pos, HashMap<String,float[]> hashNormalMethyRate, HashMap<String,float[]> hashTumorMethyRate) throws IOException
    {
        //double allCpGSites = (Double) list.get(0);
        //double methlyateRate = (Double) list.get(1);
        //List<Integer> methyNum = (List<Integer>) list.get(2);
       // List<Integer> unmethyNum  = (List<Integer>) list.get(3);
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

        double tumorlog = Math.log(tumorData) / Math.log(10);
        double normallog = Math.log(normalData) / Math.log(10);

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

        double1.add(tumorlog);
        double1.add(normallog);

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

}
