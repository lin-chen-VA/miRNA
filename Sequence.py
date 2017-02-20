#!/usr/bin/python

import csv;
import re;
import urllib2;
import sys;
from bs4 import BeautifulSoup;

def read(fileName):
    ''' Read csv file
    Args:
        fileName (string), csv file name
    Return
        2d array
    '''
    with open(fileName, 'rU') as csvfile:
        reader = csv.reader(csvfile);
        l = list(reader);
    return l;

def getTable(l, startIndex, col1, col2):
    ''' Process 2d array
    Args:
        l (list), 2d array of miRNA records
        startIndex (int), index of starting row
        col1 (int), a specific column
        col2 (int), a specific column
    '''
    miRNA = [];
    for index, ele in enumerate(l):
        if index >= startIndex:
            miRNA.append([ele[col1], ele[col2]]);
    return miRNA;

def processMiRNATable(miRNATable, outputFile):
    ''' process miRNA table
    Args:
        miRNATable (list), each element contains transcript ID and Fold-Change
        outputFile (string), name of output file
    '''
    num = len(miRNATable);
    index = 0;
    with open(outputFile, 'wb') as f:
        line = "Transcript ID, Fold Change, Target Gene, RefSeq, Gene Symbol, TargetScan\n";
        f.write(line);
        for id, change in miRNATable:
            table = parseMiRBase('http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc='+id);
            for row in table:
                line = id+','+change;
                for ele in row:
                    line += ','+ele;
                line += '\n';
                f.write(line);
            index += 1;
            print 'Have done', '%10.1f%%'%(float(index)/num*100), id

def parseMiRBase(url):
    ''' parse the webpage at mirbase.org to get the target url at targetscan.org
    Args:
        url (string), url at mirbase.org
    Return:
        table (list), a table of target genes, representative transcript, and Gene Symbol
    '''
    table = [];
    response = urllib2.urlopen(url);
    soup = BeautifulSoup(response, 'html.parser');
    for r in soup.find_all('li'):
        if re.search(r'TARGETSCAN-VERT', str(r.contents[0]), re.M|re.I):
            for target in r.find_all('a'):
                #print target
                l = matching(target.get('href'), expressionTable);
                table += l;
    for index, row in enumerate(table):
        table[index].append(target.text);
    return table;

def matching(url, expressionTable):
    ''' get all matching representative transcripts and target genes
    Args:
        url (string), target url at targetscan.org
        expressionTable (dict), expression table
    Return:
        l (list), each element contains target gene, matching representative transcript, and Gene Symbol
    '''
    l = []
    inputs = urllib2.urlopen(url);
    soup = BeautifulSoup(inputs, 'html.parser');
    #table = soup.find('table', {'id': 'restable'});
    tds = soup.find_all('td');
    for index, td in enumerate(tds):
        if td.find('a'):
            s = td.find('a').text;
            if expressionTable.has_key(s.strip()):
                l.append([tds[index-1].find('a').text, tds[index].find('a').text, expressionTable[s.strip()]]);
    return l;

def list2dict(expressionTable):
    ''' convert 2d array to a dict
    Args:
        expressionTable (list), two column array, first column is Gene Symbol, second is RefSeq
    Return:
        dict, key is RefSeq, value is Gene Symbol
    '''
    d = {};
    for row in expressionTable:
        #if d.has_key(row[1].strip()):
            #d[row[1].strip()].append(row[0].strip());
        #else:
            #l = [];
            #l.append(row[0].strip());
        d[row[1].strip()] = row[0].strip();
    return d;

if __name__ == '__main__':
    l = read(sys.argv[1]);
    miRNATable = getTable(l, 1, 3, 4);
    del l;
    e = read(sys.argv[2]);
    expressionTable = getTable(e, 1, 3, 4);
    del e;
    expressionTable = list2dict(expressionTable);

    processMiRNATable(miRNATable, 'output.csv');

