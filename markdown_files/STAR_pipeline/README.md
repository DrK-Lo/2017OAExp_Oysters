# STAR and RSEM Pipeline

### Overview

This pipeline takes advantage of a genome mapper STAR, which performs transcript alignment by mapping to a reference genome. Importantly, STAR is suited for the de novo discovery of splice junctions, which can be leveraged for identifying novel exons and isoforms. This pipeline couples the STAR mapper with RSEM for transcript quantification. This approach attempts to probabilistically estimate transcript abundance rather than simply count the reads.This may be beneficial for improving transcript count estimates, by probabilistically resolving reads which map to multiple genes (multimappers).  

## Table of Contents

1. [Brief Description and Literature on Required Tools and Scripts](#one)
2. [Step 1 - Creating STAR index](#two)
3. [Step 2 - Mapping with STAR](#three)
4. [Step 3 - Running RSEM](#four)

### Brief Description and Literature on Required Tools and Scripts <a name="one"></a>

**Mapping**

*STAR* - . 

* [Website]()  
* [Publication]()

**Transcript Quantification**

*RSEM* - .

* [Website]()
* [Publication]()

### Step 1 - Creating STAR index <a name="two"></a>

### Step 2 - Mapping with STAR <a name="three"></a>

### Step 3 - Running RSEM <a name="four"></a>
