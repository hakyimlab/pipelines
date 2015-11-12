#!/usr/bin/env python

SOURCE_DIR = 'input'
TARGET_DIR = 'output'

import gzip
import os
import sys
import sqlite3

def smart_open(source_file):
    if source_file.endswith('.txt'):
        return open(source_file)
    elif source_file.endswith('.gz'):
        return gzip.open(source_file)
    else:
        print "error: source file names should end in .txt or .gz; %s doesn't comply. exiting."%source_file
        sys.exit(1)

def smart_list(source_dir, including):
    source_files = [ x for x in os.listdir(source_dir) if including in x]
    if len(source_files) == 0:
        print "warning: no recognized source files (i.e., including %s) on %s"%(including, source_dir)
    return sorted(source_files)

def generate_weights_file():

    def data_rows_in(source_file):
        "Iterate over data rows in the source file, labeling fields and converting formats as required."
        def upconvert(x):
            for f in (int, float):
                try:
                    return f(x)
                except ValueError:
                    pass
            return x
        header = None
        for k, line in enumerate(smart_open(source_file)):
            if header is None:
                header = line.strip().split()
            else:
                yield dict(zip(header, map(upconvert, line.strip().split())))

    def source_files(source_dir=SOURCE_DIR):
        "List all relevant source files"
        for x in smart_list(source_dir, including='.allBetas.'):
            yield os.path.join(source_dir, x)

    class DB:
        "This encapsulates a single SQLite DB (for a given source file and alpha)."
        def __init__(self, source_file, alpha, target_dir=TARGET_DIR):
            tissue_name = os.path.basename(source_file).split('.')[0]
            db_filename = os.path.join(target_dir, '%s_%s.db'%(tissue_name, alpha))
            if not os.path.exists(target_dir):
                os.mkdir(target_dir)

            if os.path.exists(db_filename):
                os.unlink(db_filename)

            self.connection = sqlite3.connect(db_filename)

            self("CREATE TABLE weights (rsid TEXT, gene TEXT, weight DOUBLE, ref_allele CHARACTER, eff_allele CHARACTER, pval DOUBLE, N INTEGER, cis INTEGER)")
            self("CREATE INDEX weights_rsid ON weights (rsid)")
            self("CREATE INDEX weights_gene ON weights (gene)")
            self("CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")


        def __call__(self, sql, args=None):
            c = self.connection.cursor()
            if args:
                c.execute(sql, args)
            else:
                c.execute(sql)

        def close(self):
            self.connection.commit()

        def insert_row(self, row):
            self("INSERT INTO weights VALUES(?, ?, ?, ?, ?, NULL, NULL, NULL)", (row['rsid'], row['gene'], row['beta'], row['ref'],row['alt']))
            "alt allele is the dosage/effect allele in GTEx data"

    class MetaDB:
        "This handles all the DBs for each source file (tissue type)"
        def __init__(self, source_file):
            self.source_file = source_file
            self.dbs = {} # alpha -> DB object

        def insert_row(self, row):
            alpha = row['alpha']
            if alpha not in self.dbs:
                self.dbs[alpha] = DB(self.source_file, alpha)
            self.dbs[alpha].insert_row(row)

        def close(self):
            for db in self.dbs.values():
                db.close()



    for source_file in source_files():
        print "Processing %s..."%source_file
        meta_db = MetaDB(source_file=source_file)
        for row in data_rows_in(source_file):
            meta_db.insert_row(row)
        meta_db.close()


def add_extra_data():
    def data_rows_in(source_file):
        "Iterate over data rows in the source file, labeling fields and converting formats as required."
        def upconvert(x):
            for f in (int, float):
                try:
                    return f(x)
                except ValueError:
                    pass
            return x

        header = None
        for k, line in enumerate(smart_open(source_file)):
            if header is None:
                if 'gene' not in line: ## this is then not a header
                    header = "gene    ensid   mean.cvm        var.cvm lambda.var      lambda.frac.diff        mean.lambda.iteration   lambda.min      n.snps  R2      alpha   pval".split()
                else:
                    header = line.strip().split() ## we expect: gene    alpha   cvm     lambda.iteration        lambda.min      n.snps  R2      pval    genename
            else:
                if line.strip(): ###and 'gene' not in line: # some files, but not all, have a header
                    ret = dict(zip(header, map(upconvert, line.strip().split())))
                    if 'alpha' in ret:
                        yield ret

    def source_files(source_dir=SOURCE_DIR):
        "List all relevant source files"
        for x in smart_list(source_dir, including='.allResults.txt'):
            yield os.path.join(source_dir, x)

    class DB:
        "This encapsulates a single SQLite DB (for a given source file and alpha)."
        def __init__(self, source_file, alpha, target_dir=TARGET_DIR):
            tissue_name = os.path.basename(source_file).split('.')[0]
            db_filename = os.path.join(target_dir, '%s_%s.db'%(tissue_name, alpha))
            print(db_filename) ## delete after debug finished
            assert(os.path.exists(db_filename))
            self.connection = sqlite3.connect(db_filename)

            self("DROP INDEX IF EXISTS extra_gene")
            self("DROP TABLE IF EXISTS extra")
            self("CREATE TABLE extra (gene TEXT, genename TEXT, R2 DOUBLE,  `n.snps` INTEGER)")
            self("CREATE INDEX extra_gene ON extra (gene)")


        def __call__(self, sql, args=None):
            c = self.connection.cursor()
            if args:
                c.execute(sql, args)
            else:
                c.execute(sql)

        def close(self):
            self.connection.commit()

        def insert_row(self, row):
            if 'genename' in row:
                genename = row['genename']
            else:
                genename = row['gene']
            self("INSERT INTO extra VALUES(?, ?, ?, ?)", (row['gene'], genename, row['R2'], row['n.snps']))

    class MetaDB:
        "This handles all the DBs for each source file (tissue type)"
        def __init__(self, source_file):
            self.source_file = source_file
            self.dbs = {} # alpha -> DB object

        def insert_row(self, row):
            alpha = row['alpha']
            if alpha not in self.dbs:
                self.dbs[alpha] = DB(self.source_file, alpha)
            self.dbs[alpha].insert_row(row)

        def close(self):
            for db in self.dbs.values():
                db.close()


    for source_file in source_files():
        print "Processing %s..."%source_file
        meta_db = MetaDB(source_file=source_file)
        for row in data_rows_in(source_file):
            meta_db.insert_row(row)
        meta_db.close()


if __name__ == '__main__':
    generate_weights_file()
    add_extra_data()
