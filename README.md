# PeptoVar

PeptoVar - **Pept**ides **o**f **Var**iations

## Overview

PeptoVar is a program for personalization of protein coding genes and peptidomes generation.

 - Easy to use (see *Usage* section) 
 
 - Generate peptidomes:
   - personalized for selected sample
   - unique peptides for pair of samples
   - all possible peptides for variations in population
   - *... multi-processing version will be available soon*

- Efficiently remove synonymous variations

- Has optional parameters for variation filtration and transcript selection


## Installation / Download

#### Using Homebrew on Mac OS X or Linux (linuxbrew)

    *... will be available soon*
    
to upgrade already installed PeptoVar to the newest version:

    brew update
    brew upgrade PeptoVar

#### Manual install

* download latest stable PeptoVar build from [release page](https://github.com/DMalko/PeptoVar/releases/latest)
* unzip the archive
* add resulting folder to your ``PATH`` variable
  * or add symbolic link for ``PeptoVar`` script to your ``bin`` folder
  * or use peptoVar directly by specifying full path to the executable script

#### Requirements

* Linux or MacOS
* python >= 3.5
* pysam module >= 0.11.2.2
 
## Usage

#### Peptides for a sample

This example illustrates usage for sample S00001:

    peptoVar -samples S00001 -peptlen 9 -var nonsyn -gff ./testdata/test.gff -vcf ./testdata/test.vcf.gz


#### Unique peptides for a pair of samples
This example illustrates usage for the pair of samples S00001 and S00002:

    peptoVar -samples S00001 S00002 -peptlen 9 -var nonsyn -gff ./testdata/test.gff -vcf ./testdata/test.vcf.gz


#### Generation of peptides using all possible variations in the set
This example illustrates usage virtual sample with all variations:

    peptoVar -peptlen 9 -var nonsyn -gff ./testdata/test.gff -vcf ./testdata/test.vcf.gz
    

## Documentation

Detailed documentation can be found at https://...

If you haven't found the answer to your question in the docs, or have any suggestions concerning new features, feel free to create an issue here, on GitHub, or write an email to dmitry.malko at gmail.com .

## Build

Dependancy:

- pysam module

To build PeptoVar from source:

- Clone repository

  ```
  git clone https://github.com/DMalko/peptoVar/PeptoVar.git
  ```

- Refresh git submodules

  ```
  git submodule update --init --recursive
  ```
  
- That is all. To generate peptide you need genome annotation in GFF format and VCF file with index. (... from ENSEMBL).


## License
Copyright (c) 2017, Dmitry Malko
All Rights Reserved

PeptoVar is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).

## Cite

has not published yet

## Files referenced in original paper

Can be found [here](https://github.com/...).
