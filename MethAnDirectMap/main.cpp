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

//! this function writes a single nucleotide to a stream
void writeNucleotide(nucleotide &nuc, QTextStream &out)
{
    out << nuc.chrom << '\t' <<
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
        for (int i = 1; i < nuc.bestmappingfeatures.count(); ++i) {
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

//! this function loads, maps and writes the nucleotides - also does the stats
bool readMapWriteNucleotides(QString &nucfile, database &anno, QString &resultfile, uint &minReads)
{
    // variables
    bool rval = true;
    nucleotide nuc; // the current nucleotide
    int nuccounter = 0;

    // for the stats
    QMap<QPair<QString, QString>, QVector<float> > stats; // the pair contains: context and a bestmappingfeature
    QPair<QString, QString> curkey;
    QVector<float> emptyvec;
    float weight;
    float METHweight;
    float UNMETHweight;
    emptyvec.fill(0,4); // contains number of: methC unmethC methReads unmethReads

    // open in and out
    QFile infile(nucfile);
    if (!infile.open(QFile::ReadOnly | QFile::Text)) {
        std::cerr << "Error: Cannot read file " << qPrintable(nucfile)
                  << ": " << qPrintable(infile.errorString())
                  << std::endl;
        rval = false;
        return(rval);
    }
    QFile outfile(resultfile);
    if (!outfile.open(QFile::WriteOnly | QFile::Text)) {
        std::cerr << "Error: Cannot write file " << qPrintable(resultfile)
                  << ": " << qPrintable(outfile.errorString())
                  << std::endl;
        rval = false;
        return(rval);
    }

    // in and outstream
    QTextStream in(&infile);
    in.setCodec("UTF-8");
    QTextStream out(&outfile);
    out.setCodec("UTF-8");

    // read, map, write
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
        // add to stats
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
        // write the nucleotide
        writeNucleotide(nuc, out);
        // just some time info
        ++nuccounter;
        if ( (nuccounter % 5000000) == 0 ) {
            message = QObject::tr("Processed %1 nucleotides that were not skipped").arg(QString::number(nuccounter));
            anno.print_time(message);
        }
    }
    out.flush(); // make sure all is written
    message = QObject::tr("Processed %1 nucleotides that were not skipped").arg(QString::number(nuccounter));
    anno.print_time(message);

    infile.close();
    if (infile.error() != QFile::NoError) {
        std::cerr << "Error: Cannot read file " << qPrintable(nucfile)
                  << ": " << qPrintable(infile.errorString())
                  << std::endl;
        rval = false;
        return(rval);
    }
    outfile.close();
    if (outfile.error() != QFile::NoError) {
        std::cerr << "Error: Cannot write file " << qPrintable(resultfile)
                  << ": " << qPrintable(outfile.errorString())
                  << std::endl;
        rval = false;
        return(rval);
    }

    // print the stats
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
    std::cout << std::flush;

    return(rval);
}



//! the function that is called if something with the arguments is wrong
void usage() {
    std::cerr << std::endl <<
                 "Usage is ./programname [-S nucleotides] [-A annotation] [-O results]" << std::endl << std::endl << std::endl <<
                 "Options" << std::endl << std::endl << '\t' <<
                 "-m <uint> minimal number of reads/coverage (total)." << std::endl << '\t' << '\t' <<
                 "default: 0" << std::endl << std::endl << '\t' <<
                 "Notes" << std::endl << std::endl << '\t' <<
                 "the annotation file can be obtained via p502 format wizard (xml)." << std::endl << std::endl << '\t' <<
                 "the nucleotides file contains a tab-separated line with the fields:" << std::endl << '\t' << '\t' <<
                 "chromosome" << std::endl << '\t' << '\t' <<
                 "position" << std::endl << '\t' << '\t' <<
                 "context" << std::endl << '\t' << '\t' <<
                 "coverage" << std::endl << '\t' << '\t' <<
                 "percent methylated" << std::endl << std::endl << '\t' <<
                 "IMPORTANT: remove any kind of header before" << std::endl << std::endl << '\t' <<
                 "IMPORTANT: position counting needs to start with 0 (but they can be discontinuous)" << std::endl << std::endl << '\t' <<
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
    if (argc < 3) { usage(); }

    // variables from command line
    QString nucfile = "";
    QString annofile = "";
    QString resultfile = "";
    uint minReads = 0;
    while ((c = getopt(argc, argv, "m:S:A:O:")) != -1) {
        switch(c) {
        case 'm':
            minReads = atoi(optarg);
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

    // create the database
    QVector<QVariant> headers;
    headers << "Sname" << "Schrom" << "Sstrand" << "Ustart" << "Uend" << "Sfeature" << "SassembledFeature" << "Upriority";
    database anno(headers);
    anno.print_time("START");
    if ( anno.readData(annofile) ) { anno.print_time("annotation loaded"); }
    else {
        std::cerr << "ERROR: could not initialize database" << std::endl << std::flush;
        exit(8);
    }

    // read-map-write nucleotides
    if ( readMapWriteNucleotides(nucfile, anno, resultfile, minReads) ) { anno.print_time("read, mapped and wrote nucleotides"); }
    else { anno.print_time("ERROR: could not read/map/write nucleotides"); exit(8); }

    // give the end
    anno.print_time("END");

    //return a.exec();
    return(0);
}





