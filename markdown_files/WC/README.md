# Adult Exposure 2017 - Water Chemistry

### Collection

Tank chemistry was measured three times a week  (M,W,F) for temperature, salinity, and pH. Every two weeks (additional water sampling around the start of the exposure) borosilicate bottles were used to sample the water chemistry to measure DIC and AT to calculate the complete carbonate chemistry. Bottle samples were poisined with mercuric choloride and refridgerated prior to sampling.

### Instruments

* `Temperature` : Glass thermometer
* `Salinity` :  YSI 3200 conductivity probe (precision = 0.1 ppt)
* `pH` : Accumet solid state pH electrode (precision = 1mV)  
* `DIC` and `AT` : VINDTA 3C coupled alkalinity gram titration and coulometric DIC analyzer system


### Analysis

`pH`
  * Calculated from raw mV values using a two-point standard curve to calculate the pH, using NBS standard for pH7 and pH10.

`Complete Carbonate Chemistry`
1) `AT` was initially determined by titration and `DIC using a coulormetric approach integrated into the VINDTA 3C.
2) `**FILL IN ADDITONAL DETIAL ABOUT THE SPREADSHEET CORRECTIONS**
3) Corrected alkalinity, `AT`, and dissolved inorganic carbon, `CT` or `DIC`, was used along with the tank temperature (`temp_out`), VINDTA temperature (`temp_in`), tank salinity, tank pressure, and VINDTA pressure were used in an excel plugin for `CO2Sys`.
  * `CO2Sys Parameters`
    * `Roy et al. 1993` :  Constants
    * `Dickson` : KHSO4
    * `Seawater` : Output pH scale
    * `Lee et al. 2010` : BI Value
    
### Output Values

| Parameter | Collection | Description | Units |
|:---------:|:----------:|:-----------:|:-----:|
| pH_measured | M,W,F  | measured pH | unitless |
| Salinity | M,W,F | salinity | PSU |
| CT | Biweekly | dissolved inorganic carbon | |
| AT | Biweekly | Alkalinity | |
| pH_out | Biweekly | calculated pH | unitless |

