# SubsistenceStrategyTradeoffs

To replicate workflow and analysis of Bird et al. "A first empirical analysis of population stability in North America using radiocarbon records." First, what the code does and does not do:

The code is intended to be run with a csv of radiocarbon dates that includes several headings (lab number, sampling unit (heading SBox), conventional radiocarbon date, one sigma standard deviation, Country, Site, SiteID (generated in excel by sorting site alphabetically) and source). Any data can be input, but the headings and type of data within fit certain criteria. Radiocarbon_data.csv provides the information used for Bird et al.'s publication with locational information removed, as required by many SHPO's in the United States. 

This code does NOT clean radiocarbon data. It does not ensure that duplicate lab numbers, invalid labs, impossible dates, undesireable standard deviations, non-archaeological dates etc. are removed. It does not verify that the locational information provided for the radiocarbon age is physically possible (e.g. under the ocean or on bedrock) or accurate (e.g. that a site from a certain state has a trinomial matching that state). It DOES calculate the average date for multiple dates from the same archaeological site in an effort to correct for sampling bias. It also normalizes the SPD's (divides all SPD values by the number of dates used to produce the SPD) to make the SPD results comparable.

You may start at step 4 to run the dataset used in the publication, or at step 1 to create your own dataset.


STEP 1:
To create your own dataset, first you must acquire and clean a radiocarbon dataset. You may use more than one dataset as well. Sample the dataset using your chosen sampling method (the goal of this code is to analyze multiple sampling units, not simply one). 

Step 2:
Create a .csv named "radiocarbon_data.csv" with the following column headings: "labnumber" "SAMPLING UNIT" "date" "sd" "Lat" "Long" "Country" "Site" "SiteID" "Source", where labnumber is an alpha-numeric (dashes allowed) unique identifier, SAMPLING UNIT is your chosen sampling unit name with different unit names within, date is a positive integer for years before present (1950), sd is ONE standard deviation for years BP as provided by the radiocarbon lab, Country is the country from which the data came, Site is a unique identifier for the site shared by all radiocarbon ages from the same site, and Source is the origin of each radiocarbon date (either databases or references). In excel or a similar program, sort first by lat, then long, then Site. Under the SiteID column, type "1" at the first non-heading row (presumably row 2). In the following row in the same column, type "=if(H2 = H3, H2, H2+1)" where H2 is the cell you just typed "1". This provides a unique numeric identifier for each site that will be better handled by the sampling bias code. Copy this column and replace with values. You can now re-sort the dataset without affecting the siteID.

STEP 3:
Open the R-script in RStudio. Type "Control-F" to open a search menu. In the "Find" box, type "Sbox". In the "Replace" box, type "SAMPLING UNIT" according to your column heading in your "radiocarbon_data.csv" file.

STEP 4 (Start here if you're using the original dataset):
Open the R-script in RStudio. Beginning at the top of the script, click the "source" button to run the whole script in order. You may click "run" for each line. Note that groups of lines with no spaces can be highlighted run together, but the entire script cannot be highlighted and run. Be sure to install all packages listed in the preamble of the script.  If you encounter any problems, please ensure that your column headings match those listed in Step 2.
