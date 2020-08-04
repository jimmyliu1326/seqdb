# seqdb

## Description
The shell script contains various functions for managing and constructing local sequence database.

## Usage
```
seqdb [subcommand] [options]

Available subcommands:
    build_db    construct a local sequence database
    extract     retrieve a target sequence from database
    remove      remove an existing target from database
    insert      add a new sequence to database
    replace     replace sequence of an existing target with a new sequence
    unbuild_db  deconstruct a local sequence database
    search      identify significant hits from an existing database given a query sequence
    build_tree  construct a cgSNP tree from closest genomic neighbours in the database  
```

## Installation

1. Clone the repository

```
git clone https://github.com/jimmyliu1326/seqdb.git
```

2. Modify file permissions

```
chmod +x seqdb/seqdb
```

3. Add seqdb to PATH by including the following line in ~/.bashrc

```
export PATH="/path/to/seqdb:$PATH"
```

## Dependecies
The following dependencies must be installed and added to PATH:

* seqkit >= 0.13.0
* mash >= 2.2.2
* phame = 1.0.3
* mashtree >= 1.1.2
* picard >= 2.23.3
* kipper >= 1.0.0
