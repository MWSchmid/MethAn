#include <QtCore/QCoreApplication>
#include <QtCore>
#include "../../Rcount/source_V2/p502_SOURCE/dataStructure/databaseitem.h"
#include "../../Rcount/source_V2/p502_SOURCE/dataStructure/database.h"
#include "../../Rcount/source_V2/p502_SOURCE/dataStructure/mappingtreeitem.h"
#include <iostream>
#include <getopt.h>
#include <QTime>

//! very small class for a nucleotide to hold the information together
// note: we create just once, use setInit to write the info from the file, send it through annotation mapping, then append it to the vector
class nucleotide {
public:
    // set at beginning
    QString chrom;
    uint position;
    uint methylated;
    uint unmethylated;
    float coverage;
    bool ismethylated;
    float percentMethylated;
    QString strand;
    QString context;
    QStringList otherData;
    // set in getNucMapping
    QList<databaseItem*> bestmappingfeatures; //holds the pointers to the data vectors

    nucleotide(): strand(".") { }

    void setInit(QString &chrom, QString &position, QString &context, QString &coverage, QString &percentMethylated) {
        this->coverage = coverage.toFloat();
        this->percentMethylated = percentMethylated.toFloat();
        this->chrom = chrom;
        this->strand = ".";
        this->position = position.toUInt();
        this->methylated = floor(this->coverage*(this->percentMethylated/100.0));
        this->unmethylated = this->coverage-this->methylated;
        this->context = context;
        this->otherData.clear();
        this->otherData = otherData;
        this->bestmappingfeatures.clear();
        this->ismethylated = this->percentMethylated > 50.0;
    }

    bool hasSufReads(uint &minReads) {
        if (this->coverage >= minReads) { return(true); }
        else { return(false); }
    }
};

//! very small class for a region to hold the information together
// note: we create just once, use setInit to write the new info and then append it to the vector
class region {
public:
    // set at beginning
    QString chrom;
    uint start;
    uint end;
    QMap<QString, uint> contexts;
    uint nucs;
    // set in getRegionMapping
    QList<nucleotide*> nucleotides;

    region() {
        this->contexts.insert("CG", 0);
        this->contexts.insert("CHG", 0);
        this->contexts.insert("CHH", 0);
    }

    void setInit(nucleotide *nuc) {
        this->chrom = nuc->chrom;
        this->start = nuc->position;
        this->end = nuc->position;
        this->contexts.clear();
        this->nucleotides.clear();
        this->nucleotides << nuc;
        QMap<QString, uint>::iterator iter = this->contexts.begin();
        while (iter != this->contexts.end()) { iter.value() = 0; }
        this->nucs = 1;
        ++this->contexts[nuc->context];
    }
};

//! this function maps a nucleotide to genomic features
void getNucMapping(nucleotide &nuc, database &anno)
{
    //! TODO: needs to go via pointers - otherwise there is tooooooooo much data
    //! query speed could be optimized
    // get the full recursive mapping but using priorites
    QVector<databaseItem*> mapping = anno.bestRmapPosition(nuc.chrom, nuc.position);

    // now fill in the mappings and the best mapping variables of nuc
    if (mapping.isEmpty()) {
        nuc.bestmappingfeatures.clear();
    }
    else {
        foreach (databaseItem* element, mapping) {
            //std::cerr << element->data(0).toString().toStdString() << std::endl << std::flush;
            nuc.bestmappingfeatures << element;
        }
    }
}

//! this function loads the nucleotides, and performs the context and the annotation mappings
bool loadNucleotides(QList<nucleotide> &nucleotides, QString &nucfile, database &anno, uint &minReads)
{
    // variables
    bool rval = false;
    nucleotide nuc; // the current nucleotide
    int nuccounter = 0;


    // read the file and do all the things
    QFile file(nucfile);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        std::cerr << "Error: Cannot read file " << qPrintable(nucfile)
                  << ": " << qPrintable(file.errorString())
                  << std::endl;
        rval = false;
    }
    else {
        QTextStream in(&file);
        in.setCodec("UTF-8");

        QString curline;
        QString message;
        while(!in.atEnd()) {
            curline = in.readLine().trimmed();
            QStringList fields = curline.split('\t');
            // set the original nucleotide
            nuc.setInit(fields[0], fields[1], fields[2], fields[3], fields[4]);
            if (fields.length() > 5) {
                for (int i = 5; i < fields.size(); ++i) {
                    nuc.otherData << fields.at(i);
                }
            }
            // skip if not enough reads
            if ( !nuc.hasSufReads(minReads) ) { continue; }
            // do the mapping
            getNucMapping(nuc, anno);
            // append the nucleotide
            nucleotides.append(nuc);
            // just some time info
            ++nuccounter;
            //std::cerr << curline.toStdString() << std::endl << std::flush;
            //if (nuccounter == 100) {break;}
            if ( (nuccounter % 5000000) == 0 ) {
                message = QObject::tr("Processed %1 nucleotides that were not skipped").arg(QString::number(nuccounter));
                anno.print_time(message);
            }
        }
        message = QObject::tr("Processed %1 nucleotides that were not skipped").arg(QString::number(nuccounter));
        anno.print_time(message);


        file.close();
        if (file.error() != QFile::NoError) {
            std::cerr << "Error: Cannot read file " << qPrintable(nucfile)
                      << ": " << qPrintable(file.errorString())
                      << std::endl;
            rval = false;
        }
        else { rval = true; }

    }

    return(rval);
}
//! this function loads the nucleotides, and performs the context and the annotation mappings
bool loadRandomNucleotides(QList<nucleotide> &nucleotides, QString &randomProbs, QString &nucfile, database &anno, uint &minReads)
{
    QTime time = QTime::currentTime();
    qsrand((uint)time.msec());
    QStringList randomProbsList = randomProbs.split(",");
    QMap<QString,float> nucProbs;
    nucProbs["CG"] = randomProbsList.at(0).toFloat()*10000000;
    nucProbs["CHG"] = randomProbsList.at(1).toFloat()*10000000;
    nucProbs["CHH"] = randomProbsList.at(2).toFloat()*10000000;

    // variables
    bool rval = false;
    nucleotide nuc; // the current nucleotide
    int nuccounter = 0;

    // read the file and do all the things
    QFile file(nucfile);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        std::cerr << "Error: Cannot read file " << qPrintable(nucfile)
                  << ": " << qPrintable(file.errorString())
                  << std::endl;
        rval = false;
    }
    else {
        QTextStream in(&file);
        in.setCodec("UTF-8");

        QString curline;
        QString message;
        while(!in.atEnd()) {
            curline = in.readLine().trimmed();
            QStringList fields = curline.split('\t');
            // check if "significant"
            if ((qrand() % 10000000) > nucProbs.value(fields[2])) { continue; }
            // set the original nucleotide
            nuc.setInit(fields[0], fields[1], fields[2], fields[3], fields[4]);
            if (fields.length() > 5) {
                for (int i = 5; i < fields.size(); ++i) {
                    nuc.otherData << fields.at(i);
                }
            }
            // skip if not enough reads
            if ( !nuc.hasSufReads(minReads) ) { continue; }
            // do the mapping
            getNucMapping(nuc, anno);
            // append the nucleotide
            nucleotides.append(nuc);
            // just some time info
            ++nuccounter;
            if ( (nuccounter % 5000000) == 0 ) {
                message = QObject::tr("Processed %1 nucleotides that were not skipped").arg(QString::number(nuccounter));
                anno.print_time(message);
            }
        }
        message = QObject::tr("Processed %1 nucleotides that were not skipped").arg(QString::number(nuccounter));
        anno.print_time(message);


        file.close();
        if (file.error() != QFile::NoError) {
            std::cerr << "Error: Cannot read file " << qPrintable(nucfile)
                      << ": " << qPrintable(file.errorString())
                      << std::endl;
            rval = false;
        }
        else { rval = true; }

    }

    return(rval);
}

bool loadRandomNucleotidesNoMapping(QList<nucleotide> &nucleotides, QString &randomProbs, QString &nucfile, database &anno, uint &minReads)
{
    QTime time = QTime::currentTime();
    qsrand((uint)time.msec());
    QStringList randomProbsList = randomProbs.split(",");
    QMap<QString,float> nucProbs;
    nucProbs["CG"] = randomProbsList.at(0).toFloat()*10000000;
    nucProbs["CHG"] = randomProbsList.at(1).toFloat()*10000000;
    nucProbs["CHH"] = randomProbsList.at(2).toFloat()*10000000;

    // variables
    bool rval = false;
    nucleotide nuc; // the current nucleotide
    int nuccounter = 0;

    // read the file and do all the things
    QFile file(nucfile);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        std::cerr << "Error: Cannot read file " << qPrintable(nucfile)
                  << ": " << qPrintable(file.errorString())
                  << std::endl;
        rval = false;
    }
    else {
        QTextStream in(&file);
        in.setCodec("UTF-8");

        QString curline;
        QString message;
        while(!in.atEnd()) {
            curline = in.readLine().trimmed();
            QStringList fields = curline.split('\t');
            // check if "significant"
            if ((qrand() % 10000000) > nucProbs.value(fields[2])) { continue; }
            // set the original nucleotide
            nuc.setInit(fields[0], fields[1], fields[2], fields[3], fields[4]);
            if (fields.length() > 5) {
                for (int i = 5; i < fields.size(); ++i) {
                    nuc.otherData << fields.at(i);
                }
            }
            // skip if not enough reads
            if ( !nuc.hasSufReads(minReads) ) { continue; }
            // append the nucleotide
            nucleotides.append(nuc);
            // just some time info
            ++nuccounter;
            if ( (nuccounter % 5000000) == 0 ) {
                message = QObject::tr("Processed %1 nucleotides that were not skipped").arg(QString::number(nuccounter));
                anno.print_time(message);
            }
        }
        message = QObject::tr("Processed %1 nucleotides that were not skipped").arg(QString::number(nuccounter));
        anno.print_time(message);


        file.close();
        if (file.error() != QFile::NoError) {
            std::cerr << "Error: Cannot read file " << qPrintable(nucfile)
                      << ": " << qPrintable(file.errorString())
                      << std::endl;
            rval = false;
        }
        else { rval = true; }

    }

    return(rval);
}

//! if requested, the function obtains the regions, and performs the mapping
bool loadRegions(QList<nucleotide> &nucleotides, const uint &minNuc, const uint &maxDist, QList<region> &regions)
{
    // variables
    bool rval = false;
    region reg; // the current region


    // init a first region
    if (!nucleotides.isEmpty()) {
        reg.setInit(&nucleotides[0]);
    }

    // go through the nucleotides
    for (int i = 1; i < nucleotides.count(); ++i) {
        // check if the nuc belongs to previous region
        if (nucleotides.at(i).chrom != reg.chrom) {
            if (reg.nucs >= minNuc) {
                //getRegionMapping(reg, anno);
                regions.append(reg);
            }
            reg.setInit(&nucleotides[i]);
        }
        else {
            if ((nucleotides.at(i).position - reg.end) <= maxDist) {
                reg.end = nucleotides.at(i).position;
                ++reg.contexts[nucleotides.at(i).context];
                ++reg.nucs;
                reg.nucleotides << &nucleotides[i];
            }
            else {
                if (reg.nucs >= minNuc) {
                    //getRegionMapping(reg, anno);
                    regions.append(reg);
                }
                reg.setInit(&nucleotides[i]);
            }
        }
    }
    // append the last region if suitable
    if (reg.nucs >= minNuc) {
        //getRegionMapping(reg, anno);
        regions.append(reg);
    }

    if ( !regions.isEmpty() ) {rval = true;}
    return(rval);
}

//! this function writes the nucleotides to a file
bool writeNucleotides(QList<nucleotide> &nucleotides, QString &fileName)
{
    bool rval = false;
    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        std::cerr << "Error: Cannot read file " << qPrintable(fileName)
                  << ": " << qPrintable(file.errorString())
                  << std::endl;
        rval = false;
    }
    else {
        QTextStream out(&file);
        out.setCodec("UTF-8");
        int i;

        // File writing
        foreach (nucleotide nuc, nucleotides) {
            out <<
                   nuc.chrom << '\t' <<
                   nuc.strand << '\t' <<
                   nuc.position << '\t' <<
                   nuc.methylated << '\t' <<
                   nuc.unmethylated << '\t' <<
                   nuc.context << '\t';
            if (nuc.bestmappingfeatures.isEmpty()) {
                out << "none,intergenic";
            }
            else {
                out << nuc.bestmappingfeatures.at(0)->topParent()->data(0).toString() << ',' << nuc.bestmappingfeatures.at(0)->data(5).toString(); //the locusname and the feature
                for (i = 1; i < nuc.bestmappingfeatures.count(); ++i) {
                    if (nuc.bestmappingfeatures.at(i)->topParent()->data(0).toString() != nuc.bestmappingfeatures.at(i-1)->topParent()->data(0).toString()) { // write only one variant
                        out << '|' << nuc.bestmappingfeatures.at(i)->topParent()->data(0).toString() << ',' << nuc.bestmappingfeatures.at(i)->data(5).toString(); //the locusname and the feature
                    }
                }
            }
            if (!nuc.otherData.isEmpty()) {
                foreach (QString entry, nuc.otherData) {
                    out << '\t' << entry;
                }
            }
            out << '\n';
        }
        out.flush();

        file.close();
        if (file.error() != QFile::NoError) {
            std::cerr << "Error: Cannot write file " << qPrintable(fileName)
                      << ": " << qPrintable(file.errorString())
                      << std::endl;
            rval = false;
        }
        else { rval = true; }
    }
    return(rval);
}



//! this function write the regions to a file
bool writeRegions(QList<region> &regions, QString &fileName)
{
    bool rval = false;
    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        std::cerr << "Error: Cannot read file " << qPrintable(fileName)
                  << ": " << qPrintable(file.errorString())
                  << std::endl;
        rval = false;
    }
    else {
        QTextStream out(&file);
        out.setCodec("UTF-8");
        int i;

        // File writing
        foreach (region reg, regions) {
            out <<
                   reg.chrom << '\t' <<
                   reg.nucs << '\t' <<
                   reg.start << '\t' <<
                   reg.end << '\t' <<
                   reg.contexts["CG"] << '\t' <<
                   reg.contexts["CHG"] << '\t' <<
                   reg.contexts["CHH"] << '\t';
            // only write chrom and pos of nucleotide
            out << reg.nucleotides.at(0)->chrom << ':' << reg.nucleotides.at(0)->position;
            for (i = 1; i < reg.nucleotides.count(); ++i) {
                out << ',' << reg.nucleotides.at(i)->chrom << reg.nucleotides.at(i)->position;
            }
            out << '\n';
        }
        out.flush();

        file.close();
        if (file.error() != QFile::NoError) {
            std::cerr << "Error: Cannot write file " << qPrintable(fileName)
                      << ": " << qPrintable(file.errorString())
                      << std::endl;
            rval = false;
        }
        else { rval = true; }
    }
    return(rval);
}

//! calculate some basic stats
void calculate_stats(QList<nucleotide> &nucleotides)
{
    QMap<QPair<QString, QString>, QVector<float> > stats; // the pair contains: context and a bestmappingfeature
    QPair<QString, QString> curkey;
    QVector<float> emptyvec;
    float weight;
    float METHweight;
    float UNMETHweight;
    emptyvec.fill(0,4); // contains number of: methC unmethC methReads unmethReads
    foreach (nucleotide nuc, nucleotides) {
        if (nuc.bestmappingfeatures.isEmpty()) {
            curkey = qMakePair(nuc.context, QString("intergenic")); // 5 is the featurecol
            if (!stats.contains(curkey)) { stats.insert(curkey, emptyvec); }
            if (nuc.ismethylated) { ++stats[curkey][0]; }
            else { ++stats[curkey][1]; }
            stats[curkey][2] += nuc.methylated;
            stats[curkey][3] += nuc.unmethylated;
        }
        else {
            weight = 1/static_cast<float>(nuc.bestmappingfeatures.count());
            METHweight = nuc.methylated/static_cast<float>(nuc.bestmappingfeatures.count());
            UNMETHweight = nuc.unmethylated/static_cast<float>(nuc.bestmappingfeatures.count());
            foreach (databaseItem* element, nuc.bestmappingfeatures) {
                curkey = qMakePair(nuc.context, element->data(5).toString()); // 5 is the featurecol
                if (!stats.contains(curkey)) { stats.insert(curkey, emptyvec); }
                if (nuc.ismethylated) { stats[curkey][0] += weight; }
                else { stats[curkey][1] += weight; }
                stats[curkey][2] += METHweight;
                stats[curkey][3] += UNMETHweight;
            }
        }
    }
    std::cout << "context" << '\t' <<
                 "feature" << '\t' <<
                 "methylatedCs" << '\t' <<
                 "unmethylatedCs" << '\t' <<
                 "methylatedCoverage" << '\t' <<
                 "unmethylatedCoverage" << '\t' <<
                 std::endl;
    QList< QPair<QString, QString> > allkeys = stats.keys();
    foreach (curkey, allkeys) {
        std::cout << curkey.first.toStdString() << '\t' <<
                     curkey.second.toStdString() << '\t' <<
                     stats[curkey][0] << '\t' <<
                     stats[curkey][1] << '\t' <<
                     stats[curkey][2] << '\t' <<
                     stats[curkey][3] << std::endl;
    }

}


//! the function that is called if something with the arguments is wrong
void usage() {
    std::cerr << std::endl <<
                 "Usage is ./programname [options [-R regions]] [-S nucleotides] [-A annotation] [-O results]" << std::endl << std::endl << std::endl <<
                 "Options" << std::endl << std::endl << '\t' <<
                 "-m <uint> minimal number of reads/coverage (total)." << std::endl << '\t' << '\t' <<
                 "default: 0" << std::endl << std::endl << '\t' <<
                 "-n <uint> minimal number of nucleotides within a region. Must be above 3." << std::endl << '\t' << '\t' <<
                 "default: 3" << std::endl << std::endl << '\t' <<
                 "-d <uint> maximal distance between two nucleotides/regions to form one region. Must be positive." << std::endl << '\t' << '\t' <<
                 "default: 250" << std::endl << std::endl << '\t' <<
                 "-Z <CG,CHG,CHH-probabilities> if option is given, the nucleotides are not chosen randomly (for each context a separate probability)." << std::endl << '\t' << '\t' <<
                 "default: off" << std::endl << std::endl << std::endl <<
                 "-R <file> if option is given, regions are obtained using the parameters n and d. Report written in the file" << std::endl << '\t' << '\t' <<
                 "default: no regions calculated" << std::endl << std::endl << std::endl <<
                 "Notes" << std::endl << std::endl << '\t' <<
                 "the annotation file can be obtained via p502 format wizard (xml)." << std::endl << std::endl << '\t' <<
                 "the nucleotides file contains a tab-separated line with the fields:" << std::endl << '\t' << '\t' <<
                 "chromosome" << std::endl << '\t' << '\t' <<
                 "position" << std::endl << '\t' << '\t' <<
                 "context" << std::endl << '\t' << '\t' <<
                 "coverage" << std::endl << '\t' << '\t' <<
                 "percent methylated" << std::endl << std::endl << '\t' <<
                 "IMPORTANT: the file should be sorted - at least for the regions" << std::endl << std::endl << '\t' <<
                 "IMPORTANT: remove any kind of header before" << std::endl << std::endl << '\t' <<
                 "IMPORTANT: position counting needs to start with 0 (but they can be discontinuous)" << std::endl << std::endl << '\t' <<
                 "IMPORTANT: if Z is set, A can be set to SKIP - this will do only loading and regions" << std::endl << std::endl << '\t' <<
                 "IMPORTANT: O can be set to SKIP - this will omit writing mapped nucleotides (mainly useful for the random sets)" << std::endl << std::endl << '\t' <<
                 "NOTE: methylation state and coverage will refer to the average of the ancestral samples" << std::endl << std::endl << '\t' <<
                 "the statistics are calculated via pairs. If a position has two equally important mappings, they are both counted proportionally" << std::endl << std::endl << std::flush;
    exit(8);
}

int c;
extern char *optarg;
extern int optind, optopt, opterr;

int main(int argc, char *argv[])
{
    //QCoreApplication a(argc, argv);

    // check arguments
    if (argc < 4) { usage(); }

    // variables from command line
    QString nucfile = "";
    QString annofile = "";
    QString resultfile = "";
    QString regionfile = "";
    bool useRandom = false;
    QString randomProbs = "";
    uint minNuc = 3; // minimal number of nucleotides per DMR
    uint maxDist = 250; // maximal distance from one to the next nucleotides within a DMR
    uint minReads = 0;
    while ((c = getopt(argc, argv, "t:m:n:d:DZ:R:S:A:O:")) != -1) {
        switch(c) {
        case 'n':
            minNuc = atoi(optarg);
            break;
        case 'm':
            minReads = atoi(optarg);
            break;
        case 'd':
            maxDist = atoi(optarg);
            break;
        case 'Z':
            useRandom = true;
            randomProbs = optarg;
            break;
        case 'R':
            regionfile = optarg;
            break;
        case 'S':
            nucfile = optarg;
            break;
        case 'A':
            annofile = optarg;
            break;
        case 'O':
            resultfile = optarg;
            break;
        case ':':
            std::cerr << "some stuff not specified" << std::endl << std::flush;
            break;
        case '?':
            std::cerr << "unknown argument" << std::endl << std::flush;
            usage();
        }
    }
    if (minNuc < 3) { usage(); }

    /*
    //! some test for the random stuff
    QTime fixedTime = QTime::fromString("235959", "hh'mm'ss");
    QTime time = QTime::currentTime();
    uint forSeed = qAbs(time.msecsTo(fixedTime));
    qsrand(forSeed);
    float hurz = 100000000*threshold;
    int successcounter = 0;
    for (int i = 0; i < 35208436; i++) {
        if ((qrand() % 100000000) > hurz) { continue; }
        successcounter++;
    }
    std::cerr << forSeed << "\t" << successcounter << std::endl << std::flush;
    exit(8);
    */

    // create the database
    QVector<QVariant> headers;
    headers << "Sname" << "Schrom" << "Sstrand" << "Ustart" << "Uend" << "Sfeature" << "SassembledFeature" << "Upriority";
    database anno(headers);
    anno.print_time("START");
    if (useRandom && (annofile == "SKIP")) {
        anno.print_time("skipped annotation loading");
    } else {
        if ( anno.readData(annofile) ) { anno.print_time("annotation loaded"); }
        else {
            std::cerr << "ERROR: could not initialize database" << std::endl << std::flush;
            exit(8);
        }
    }

    //anno.writeData("TestAnno.xml");

    // load the nucleotides -
    QList<nucleotide> nucleotides;
    nucleotides.reserve(50000000);
    if (useRandom) {
        if (annofile == "SKIP") {
            if ( loadRandomNucleotidesNoMapping(nucleotides, randomProbs, nucfile, anno, minReads) ) { anno.print_time("loaded nucleotides and skipped mapping"); }
            else { anno.print_time("ERROR: could not load the nucleotides"); exit(8); }
        } else {
            if ( loadRandomNucleotides(nucleotides, randomProbs, nucfile, anno, minReads) ) { anno.print_time("mapping obtained"); }
            else { anno.print_time("ERROR: could not load the nucleotides"); exit(8); }
        }
    } else {
        if ( loadNucleotides(nucleotides, nucfile, anno, minReads) ) { anno.print_time("mapping obtained"); }
        else { anno.print_time("ERROR: could not load the nucleotides"); exit(8); }
    }
    if (resultfile == "SKIP") {
        anno.print_time("skipped writing nucleotides");
    } else {
        if ( writeNucleotides(nucleotides, resultfile) ) { anno.print_time("nucleotides written"); }
        else { anno.print_time("ERROR: could not write the nucleotides"); exit(8);}
    }

    // get regions if requested
    if (!regionfile.isEmpty()) {
        QList<region> regions;
        regions.reserve(1000000);
        if ( loadRegions(nucleotides, minNuc, maxDist, regions) ) { anno.print_time("regions obtained"); }
        else { anno.print_time("ERROR: could not load the regions"); exit(8);}
        if ( writeRegions(regions, regionfile) ) { anno.print_time("regions written");}
        else { anno.print_time("ERROR: could not write the regions"); exit(8);}
    }

    // calculate some basic stats
    calculate_stats(nucleotides);

    // give the end
    anno.print_time("END");

    //return a.exec();
    return(0);
}





